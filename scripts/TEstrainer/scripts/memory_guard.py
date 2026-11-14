#!/usr/bin/env python3
"""Limit the combined RSS of a TEstrainer subprocess and its process group."""

import argparse
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import time


MIB = 1024 * 1024


class MemoryLimitExceeded(Exception):
    pass


def available_memory(requested, proc_root=Path("/proc")):
    """Return a conservative startup budget in MiB."""
    available = requested * MIB
    meminfo = proc_root / "meminfo"
    if meminfo.exists():
        for line in meminfo.read_text().splitlines():
            if line.startswith("MemAvailable:"):
                available = min(available, int(line.split()[1]) * 1024)
                break

    cgroup_file = proc_root / "self/cgroup"
    mountinfo_file = proc_root / "self/mountinfo"
    if cgroup_file.exists() and mountinfo_file.exists():
        memberships = [line.split(":", 2) for line in cgroup_file.read_text().splitlines()]
        for line in mountinfo_file.read_text().splitlines():
            before, after = line.split(" - ", 1)
            mount = before.split()
            fs = after.split()
            if fs[0] == "cgroup2":
                membership = next((parts[2] for parts in memberships if parts[1] == ""), None)
                limit_name, usage_name = "memory.max", "memory.current"
            elif fs[0] == "cgroup" and "memory" in fs[2].split(","):
                membership = next((parts[2] for parts in memberships if "memory" in parts[1].split(",")), None)
                limit_name, usage_name = "memory.limit_in_bytes", "memory.usage_in_bytes"
            else:
                continue
            if membership is None:
                continue
            root = Path(mount[3])
            try:
                relative = Path(membership).relative_to(root)
            except ValueError:
                continue
            mountpoint = Path(mount[4])
            current = mountpoint / relative
            while True:
                limit_path = current / limit_name
                usage_path = current / usage_name
                if limit_path.is_file() and usage_path.is_file():
                    limit = limit_path.read_text().strip()
                    if limit != "max":
                        remaining = int(limit) - int(usage_path.read_text().strip())
                        available = min(available, max(0, remaining))
                if current == mountpoint:
                    break
                current = current.parent
    # Reserve one eighth for the scheduler, shared memory and unmonitored overhead.
    return available * 7 // 8 // MIB


def group_rss_kib(pgid):
    result = subprocess.run(
        ["ps", "-A", "-o", "pgid=", "-o", "rss="],
        check=True, stdout=subprocess.PIPE, text=True,
    )
    total = 0
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) == 2 and int(fields[0]) == pgid:
            total += int(fields[1])
    return total


def run_guarded(command, limit_mib, label, output, fallback_source=None,
                fallback_output=None, remove_on_fallback=(), allow_empty_output=False):
    temp_path = None
    output_file = None
    process = None
    peak_rss = 0
    start = time.monotonic()

    def terminate(signum, frame):
        raise InterruptedError(f"received signal {signum}")

    previous = signal.signal(signal.SIGTERM, terminate)
    try:
        if output:
            output_file = tempfile.NamedTemporaryFile(
                mode="wb", prefix=".testrainer-", dir=os.path.dirname(os.path.abspath(output)),
                delete=False,
            )
            temp_path = output_file.name
        process = subprocess.Popen(command, stdout=output_file, start_new_session=True)
        while process.poll() is None:
            rss = group_rss_kib(process.pid)
            peak_rss = max(peak_rss, rss)
            if rss > limit_mib * 1024:
                raise MemoryLimitExceeded(
                    f"{label} exceeded its {limit_mib} MiB process-group RSS budget "
                    f"({rss // 1024} MiB observed)"
                )
            time.sleep(0.2)
        if process.returncode:
            raise RuntimeError(f"{label} exited with status {process.returncode}")
        if output_file:
            output_file.close()
            output_file = None
            if not allow_empty_output and os.path.getsize(temp_path) == 0:
                raise RuntimeError(f"{label} produced an empty output file")
            os.replace(temp_path, output)
            temp_path = None
        if label.startswith(("primary-mafft:", "trimming:")):
            print(f"{label}: peak sampled process-group RSS {peak_rss // 1024} MiB; "
                  f"elapsed {time.monotonic() - start:.1f} s", file=sys.stderr)
    except (MemoryLimitExceeded, MemoryError, RuntimeError, ValueError, OSError,
            subprocess.SubprocessError, InterruptedError, KeyboardInterrupt) as error:
        if process is not None and process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
        if isinstance(error, MemoryLimitExceeded) and fallback_source:
            fallback_temp = None
            try:
                for path in remove_on_fallback:
                    Path(path).unlink(missing_ok=True)
                if os.path.getsize(fallback_source) == 0:
                    raise ValueError(f"fallback sequence is empty: {fallback_source}")
                with tempfile.NamedTemporaryFile(
                    mode="wb", prefix=".testrainer-",
                    dir=os.path.dirname(os.path.abspath(fallback_output)), delete=False,
                ) as fallback:
                    fallback_temp = fallback.name
                    with open(fallback_source, "rb") as source:
                        shutil.copyfileobj(source, fallback)
                os.replace(fallback_temp, fallback_output)
                fallback_temp = None
            except (OSError, ValueError) as fallback_error:
                print(f"ERROR: {error}; could not retain starting sequence: {fallback_error}",
                      file=sys.stderr)
                return 1
            finally:
                if fallback_temp and os.path.exists(fallback_temp):
                    os.unlink(fallback_temp)
            print(f"WARNING: {error}; retaining the sequence entering this round, "
                  "skipping further alignment rounds, and passing it to final filtering "
                  "and classification.", file=sys.stderr)
            return 0
        print(f"ERROR: {error}", file=sys.stderr)
        if isinstance(error, MemoryLimitExceeded) and label.startswith(("primary-mafft:", "trimming:")):
            print("Try -J 1 or a larger -B within the scheduler allocation; "
                  "alternatively start a fresh run with -A fast (which can change results).",
                  file=sys.stderr)
        return 1
    finally:
        signal.signal(signal.SIGTERM, previous)
        if output_file:
            output_file.close()
        if temp_path and os.path.exists(temp_path):
            os.unlink(temp_path)
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--budget", type=int, help="Print usable budget for the requested MiB")
    parser.add_argument("--limit", type=int, help="Process-group RSS limit in MiB")
    parser.add_argument("--label")
    parser.add_argument("--output")
    parser.add_argument("--allow-empty-output", action="store_true")
    parser.add_argument("--fallback-source", help="Sequence entering this round")
    parser.add_argument("--fallback-output", help="Completed-family output on a memory breach")
    parser.add_argument("--remove-on-fallback", action="append", default=[],
                        help="Remove a partial output after a memory breach")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.budget is not None:
        if args.budget <= 0:
            parser.error("--budget must be positive")
        print(available_memory(args.budget))
        return 0
    if not args.limit or args.limit < 1 or not args.label or not args.command:
        parser.error("--limit, --label and a command are required")
    if bool(args.fallback_source) != bool(args.fallback_output):
        parser.error("--fallback-source and --fallback-output must be supplied together")
    if args.remove_on_fallback and not args.fallback_source:
        parser.error("--remove-on-fallback requires a fallback sequence")
    if args.allow_empty_output and not args.output:
        parser.error("--allow-empty-output requires --output")
    command = args.command
    if command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("a command is required")
    return run_guarded(command, args.limit, args.label, args.output, args.fallback_source,
                       args.fallback_output, args.remove_on_fallback, args.allow_empty_output)


if __name__ == "__main__":
    sys.exit(main())
