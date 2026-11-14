import os
import importlib.util
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts" / "TEstrainer" / "scripts"
GUARD = SCRIPTS / "memory_guard.py"


class MemoryGuardTests(unittest.TestCase):
    def run_guard(self, *args):
        return subprocess.run(
            [sys.executable, str(GUARD), *args],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )

    def test_combined_parent_and_child_memory_is_limited(self):
        if not shutil.which("ps"):
            self.skipTest("ps is unavailable")
        child = "import time; data=bytearray(35*1024*1024); time.sleep(3)"
        parent = (
            "import subprocess,sys,time; data=bytearray(35*1024*1024); "
            f"p=subprocess.Popen([sys.executable,'-c',{child!r}]); p.wait()"
        )
        result = self.run_guard("--limit", "60", "--label", "family", "--",
                                sys.executable, "-c", parent)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("family exceeded its 60 MiB", result.stderr)

    def test_failed_or_empty_output_does_not_replace_previous_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "alignment.fa"
            output.write_text("previous\n")
            for command in ("import sys; print('partial'); sys.exit(2)", "pass"):
                result = self.run_guard("--limit", "128", "--label", "mafft",
                                        "--output", str(output), "--",
                                        sys.executable, "-c", command)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(output.read_text(), "previous\n")
                self.assertEqual(list(Path(directory).iterdir()), [output])

    def test_successful_output_is_published(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "alignment.fa"
            result = self.run_guard("--limit", "128", "--label", "mafft",
                                    "--output", str(output), "--",
                                    sys.executable, "-c", "print('alignment')")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(output.read_text(), "alignment\n")

    def test_memory_breach_retains_sequence_and_removes_partial_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "round-input.fa"
            source.write_text(">family#Unknown\nACGT\n")
            completed = root / "complete.fa"
            partial = root / "further.fa"
            partial.write_text(">partial\nAAA\n")
            alignment = root / "alignment.fa"
            result = self.run_guard(
                "--limit", "20", "--label", "trimming:family.fa",
                "--fallback-source", str(source), "--fallback-output", str(completed),
                "--remove-on-fallback", str(partial), "--output", str(alignment), "--",
                sys.executable, "-c",
                "import time; data=bytearray(40*1024*1024); print('partial', flush=True); time.sleep(3)",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(completed.read_text(), source.read_text())
            self.assertFalse(partial.exists())
            self.assertFalse(alignment.exists())
            self.assertIn("WARNING: trimming:family.fa exceeded its 20 MiB", result.stderr)
            self.assertIn("skipping further alignment rounds", result.stderr)
            self.assertEqual(sorted(path.name for path in root.iterdir()),
                             ["complete.fa", "round-input.fa"])

    def test_fallback_does_not_hide_other_errors_or_missing_source(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, completed = root / "source.fa", root / "complete.fa"
            source.write_text(">family\nACGT\n")
            result = self.run_guard(
                "--limit", "128", "--label", "family", "--fallback-source", str(source),
                "--fallback-output", str(completed), "--", sys.executable, "-c",
                "import sys; sys.exit(7)",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("exited with status 7", result.stderr)
            self.assertFalse(completed.exists())
            source.unlink()
            result = self.run_guard(
                "--limit", "20", "--label", "family", "--fallback-source", str(source),
                "--fallback-output", str(completed), "--", sys.executable, "-c",
                "import time; data=bytearray(40*1024*1024); time.sleep(3)",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("could not retain starting sequence", result.stderr)
            self.assertFalse(completed.exists())

    def test_empty_output_can_be_allowed_for_trf(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "trf"
            result = self.run_guard("--limit", "128", "--label", "trf",
                                    "--output", str(output), "--allow-empty-output",
                                    "--", sys.executable, "-c", "pass")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(output.read_bytes(), b"")

    def test_budget_accounts_for_cgroup_ancestor_usage(self):
        spec = importlib.util.spec_from_file_location("memory_guard", GUARD)
        guard = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(guard)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            proc = root / "proc"
            (proc / "self").mkdir(parents=True)
            (proc / "meminfo").write_text("MemAvailable: 10000000 kB\n")
            mount = root / "cgroup"
            (mount / "job").mkdir(parents=True)
            (proc / "self/cgroup").write_text("0::/allocated/job\n")
            (proc / "self/mountinfo").write_text(
                f"1 0 0:1 /allocated {mount} rw - cgroup2 cgroup rw\n"
            )
            (mount / "memory.max").write_text(str(3 * 1024**3))
            (mount / "memory.current").write_text(str(1024**3))
            (mount / "job/memory.max").write_text("max")
            (mount / "job/memory.current").write_text("0")
            self.assertEqual(guard.available_memory(8192, proc), 1792)


class ConcurrencyTests(unittest.TestCase):
    def run_wrapper(self, threads, usable_mib, override=None):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "genome.fa").write_text(">chr1\nACGT\n")
            (root / "library.fa").write_text(">family\nACGT\n")
            bin_dir = root / "bin"
            bin_dir.mkdir()
            fake_python = bin_dir / "python3"
            fake_python.write_text(
                "#!/bin/sh\n"
                "case \" $* \" in *' --budget '*) echo \"$TEST_USABLE_MIB\"; exit 0;; esac\n"
                "exit 99\n"
            )
            fake_python.chmod(0o755)
            wrapper = root / "wrapper.sh"
            source = (ROOT / "scripts/TEstrainer/TEstrainer_for_earlGrey.sh").read_text()
            source = source.replace(
                "STRAIN_SCRIPTS=/data/toby/EarlGrey/scripts/TEstrainer/scripts/",
                f"STRAIN_SCRIPTS={SCRIPTS}/",
            )
            wrapper.write_text(source)
            command = ["bash", str(wrapper), "-g", str(root / "genome.fa"),
                       "-l", str(root / "library.fa"), "-d", str(root / "output"),
                       "-t", str(threads), "-q"]
            if override is not None:
                command.extend(["-J", str(override)])
            env = os.environ.copy()
            env["PATH"] = str(bin_dir) + os.pathsep + env["PATH"]
            env["TEST_USABLE_MIB"] = str(usable_mib)
            return subprocess.run(command, env=env, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, text=True)

    def test_auto_scales_with_threads_and_available_memory(self):
        for threads, memory, jobs in ((2, 32768, 2), (16, 32768, 4), (16, 4096, 1)):
            with self.subTest(threads=threads, memory=memory):
                result = self.run_wrapper(threads, memory)
                self.assertEqual(result.returncode, 99, result.stderr)
                self.assertIn(f"{jobs} alignment jobs", result.stdout)

    def test_override_respects_cpu_and_memory_caps(self):
        for threads, memory, override, jobs in ((16, 32768, 7, 7), (2, 32768, 7, 2),
                                                (16, 2048, 7, 2)):
            with self.subTest(threads=threads, memory=memory, override=override):
                result = self.run_wrapper(threads, memory, override)
                self.assertEqual(result.returncode, 99, result.stderr)
                self.assertIn(f"{jobs} alignment jobs", result.stdout)

    def test_rejects_invalid_override(self):
        result = self.run_wrapper(16, 32768, 0)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("alignment jobs must be a positive integer or auto", result.stderr)


class EntrypointTests(unittest.TestCase):
    def test_strainer_failure_does_not_publish_stale_library(self):
        for entrypoint in ("earlGrey", "earlGreyLibConstruct"):
            for function in ("strainer", "strainerResume"):
                with self.subTest(entrypoint=entrypoint, function=function), \
                        tempfile.TemporaryDirectory() as directory:
                    root = Path(directory)
                    strainer_dir = root / "sample_strainer"
                    stale_dir = strainer_dir / "TS_sample-families.fa_old"
                    stale_dir.mkdir(parents=True)
                    (stale_dir / "sample-families.fa.strained").write_text(">old\nACGT\n")
                    script_dir = root / "scripts/TEstrainer"
                    script_dir.mkdir(parents=True)
                    fake = script_dir / "TEstrainer_for_earlGrey.sh"
                    fake.write_text("#!/bin/sh\nprintf '%s\\n' \"$@\" > \"$FAKE_ARGS\"\nexit 33\n")
                    fake.chmod(0o755)
                    source = (ROOT / entrypoint).read_text()
                    functions = source[source.index("strainer()\n"):source.index(
                        "\n# Subprocess ", source.index("strainerResume()"))]
                    env = os.environ.copy()
                    env.update(OUTDIR=str(root), SCRIPT_DIR=str(root / "scripts"),
                               species="sample", genome="genome.fa", ProcNum="16",
                               Flank="100", num="1", no_seq="20", min_seq="3",
                               strainMemory="64", strainJobs="auto", strainAlignment="localpair",
                               quietBar="no", FAKE_ARGS=str(root / "arguments"))
                    result = subprocess.run(["bash", "-c", functions + "\n" + function],
                                            env=env, stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE, text=True)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("ERROR: TEstrainer failed", result.stderr)
                    self.assertFalse((strainer_dir / "sample-families.fa.strained").exists())
                    arguments = (root / "arguments").read_text().splitlines()
                    self.assertEqual(arguments[arguments.index("-J") + 1], "auto")


class MrepsTests(unittest.TestCase):
    def test_only_ambiguous_base_failure_is_tolerated(self):
        cases = (
            (0, "", True, False),
            (17, "Error: Too many N's in the window", True, True),
            (17, "Different error", False, False),
            (42, "Error: Too many N's in the window", False, False),
        )
        for status, diagnostic, succeeds, warns in cases:
            with self.subTest(status=status, diagnostic=diagnostic), \
                    tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                sequence = root / "family.fasta"
                sequence.write_text(">family\nACGTNNNN\n")
                mreps = root / "mreps"
                mreps.write_text(
                    "#!/bin/sh\n"
                    "echo '1 -> 40 : 2 3 4 5 6 7 8'\n"
                    "printf '%s\\n' \"$TEST_MREPS_ERROR\" >&2\n"
                    'exit "$TEST_MREPS_STATUS"\n'
                )
                mreps.chmod(0o755)
                env = os.environ.copy()
                env.update(PATH=str(root) + os.pathsep + env["PATH"],
                           TEST_MREPS_STATUS=str(status), TEST_MREPS_ERROR=diagnostic)
                result = subprocess.run(
                    ["bash", str(SCRIPTS / "mreps_parser.sh"), "-i", str(sequence)],
                    env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                )
                self.assertEqual(result.returncode == 0, succeeds, result.stderr)
                self.assertEqual("WARNING: mreps did not finish" in result.stderr, warns)
                output = root / "family.fasta.mreps"
                self.assertEqual(output.exists(), succeeds)
                if succeeds:
                    self.assertTrue(output.read_text())


class FamilyMemoryFallbackTests(unittest.TestCase):
    def test_round_input_survives_each_family_stage_and_other_families_continue(self):
        if not shutil.which("parallel"):
            self.skipTest("GNU parallel is unavailable")
        try:
            import Bio  # noqa: F401
        except ImportError as error:
            self.skipTest(str(error))
        helper = """#!/usr/bin/env python3
import os
from pathlib import Path
import shutil
import sys
import time

name = Path(sys.argv[0]).name
args = sys.argv[1:]
def value(flag):
    return args[args.index(flag) + 1]
def log(message):
    with open(os.environ["EVENT_LOG"], "a") as out:
        out.write(message + "\\n")
def exceed(stage, path):
    if os.environ["FAIL_STAGE"] == stage and "/run_2/" in str(path) and Path(path).name == "bad.fasta":
        data = bytearray(40 * 1024 * 1024)
        time.sleep(3)

if name == "initial_mafft_setup.py":
    run = value("-r")
    family = value("-s")
    base = Path(value("-d")) / ("run_" + run)
    log("prepare:" + run + ":" + family)
    target = base / "to_align" / family
    target.write_text((base / "raw" / family).read_text())
    exceed("preparation", target)
elif name == "TEtrim.py":
    run = value("-n")
    family = Path(value("-i")).name
    base = Path(value("-d")) / ("run_" + run)
    log("trim:" + run + ":" + family)
    target = base / ("TEtrim_complete" if run == "3" else "TEtrim_further") / family
    target.write_text(">" + family.removesuffix(".fasta") + "#Unknown\\nACGTACGT\\n")
    exceed("trimming", base / "raw" / family)
elif name == "trf":
    exceed("trf", args[0])
elif name == "mafft":
    exceed("primary-mafft", args[-1])
    sys.stdout.write(Path(args[-1]).read_text())
elif name == "makeblastdb":
    Path(value("-in") + ".nsq").touch()
elif name == "sa-ssr":
    Path(args[-1]).touch()
elif name == "Rscript":
    log("final-filter")
    library = Path(value("-i"))
    dest = Path(value("-d")) / "trf"
    shutil.copyfile(library, dest / (library.name + ".nonsatellite"))
    (dest / (library.name + ".satellites")).touch()
elif name == "RepeatClassifier":
    log("classify")
    library = Path(value("-consensi"))
    shutil.copyfile(library, str(library) + ".classified")
"""
        for stage in ("trf", "preparation", "primary-mafft", "trimming"):
            with self.subTest(stage=stage), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                scripts, bin_dir = root / "scripts", root / "bin"
                scripts.mkdir()
                bin_dir.mkdir()
                for name in ("splitter.py", "trf_parser.py"):
                    shutil.copyfile(SCRIPTS / name, scripts / name)
                guard = (SCRIPTS / "memory_guard.py").read_text().replace(
                    "    command = args.command\n",
                    '    if args.fallback_source and "/run_2/" in args.fallback_source '
                    'and args.label == os.environ["FAIL_STAGE"] + ":bad.fasta":\n'
                    "        args.limit = 20\n"
                    "    command = args.command\n",
                )
                (scripts / "memory_guard.py").write_text(guard)
                for name in ("initial_mafft_setup.py", "TEtrim.py"):
                    (scripts / name).write_text(helper)
                (scripts / "indexer.py").write_text("pass\n")
                (scripts / "mreps_parser.sh").write_text('#!/bin/sh\n: > "$2.mreps"\n')
                for name in ("trf", "mafft", "makeblastdb", "sa-ssr", "Rscript", "RepeatClassifier"):
                    executable = bin_dir / name
                    executable.write_text(helper)
                    executable.chmod(0o755)
                (root / "genome.fa").write_text(">chr1\nACGTACGT\n")
                (root / "library.fa").write_text(">bad#Unknown\nAAAA\n>good#Unknown\nCCCC\n")
                wrapper = root / "wrapper.sh"
                source = (ROOT / "scripts/TEstrainer/TEstrainer_for_earlGrey.sh").read_text()
                wrapper.write_text(source.replace(
                    "STRAIN_SCRIPTS=/data/toby/EarlGrey/scripts/TEstrainer/scripts/",
                    f"STRAIN_SCRIPTS={scripts}/",
                ))
                env = os.environ.copy()
                env.update(PATH=str(bin_dir) + os.pathsep + env["PATH"],
                           EVENT_LOG=str(root / "events"), FAIL_STAGE=stage)
                result = subprocess.run(
                    ["bash", str(wrapper), "-g", str(root / "genome.fa"),
                     "-l", str(root / "library.fa"), "-d", str(root / "output"),
                     "-t", "2", "-r", "3", "-B", "4", "-J", "1", "-q"],
                    env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=120,
                )
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn(f"WARNING: {stage}:bad.fasta exceeded its 20 MiB", result.stderr)
                self.assertIn("skipping further alignment rounds", result.stderr)
                library = (root / "output/library.fa.strained").read_text()
                self.assertEqual(library.count(">bad#Unknown\nACGTACGT\n"), 1)
                self.assertEqual(library.count(">good#Unknown\nACGTACGT\n"), 1)
                events = (root / "events").read_text().splitlines()
                self.assertNotIn("prepare:3:bad.fasta", events)
                self.assertIn("prepare:3:good.fasta", events)
                self.assertIn("final-filter", events)
                self.assertIn("classify", events)


class TEtrimTests(unittest.TestCase):
    def setUp(self):
        try:
            import Bio  # noqa: F401
            import numpy  # noqa: F401
            import pandas  # noqa: F401
        except ImportError as error:
            self.skipTest(str(error))
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.run_dir = self.root / "run_1"
        for directory in ("mafft", "TEtrim_con", "TEtrim_bp", "TEtrim_unaln",
                          "TEtrim_blast", "TEtrim_mafft", "TEtrim_complete", "TEtrim_further"):
            (self.run_dir / directory).mkdir(parents=True)
        (self.run_dir / "mafft" / "family.fa").write_text(
            ">family\nACGTACGT\n>s1\nACG-ACGT\n>s2\nACG-ACGT\n>s3\nACG-ACGT\n"
        )
        self.bin = self.root / "bin"
        self.bin.mkdir()
        blast = self.bin / "blastn"
        blast.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "if '-out' in sys.argv:\n"
            " open(sys.argv[sys.argv.index('-out')+1], 'w').write('8\\t8\\t8\\t100\\n')\n"
            "else:\n"
            " print('family\\ts1\\t100\\nfamily\\ts2\\t100\\nfamily\\ts3\\t100')\n"
        )
        mafft = self.bin / "mafft"
        mafft.write_text(
            "#!/usr/bin/env python3\n"
            "import os,sys\n"
            "open(os.environ['FAKE_MAFFT_LOG'], 'w').write(' '.join(sys.argv[1:]))\n"
            "if os.environ.get('FAKE_MAFFT_FAIL'): sys.exit(1)\n"
            "sys.stdout.write(open(sys.argv[-1]).read())\n"
        )
        blast.chmod(0o755)
        mafft.chmod(0o755)

    def run_trim(self, flag, fail=False):
        env = os.environ.copy()
        env["PATH"] = str(self.bin) + os.pathsep + env["PATH"]
        env["FAKE_MAFFT_LOG"] = str(self.root / "mafft.log")
        if fail:
            env["FAKE_MAFFT_FAIL"] = "1"
        return subprocess.run(
            [sys.executable, str(SCRIPTS / "TEtrim.py"), "-i",
             str(self.run_dir / "mafft" / "family.fa"), "-d", str(self.root),
             "-n", "1", flag],
            env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )

    def test_trim_preserves_consensus_and_selects_alignment_mode(self):
        result = self.run_trim("--fast")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--6merpair --memsave", (self.root / "mafft.log").read_text())
        self.assertIn("ACGACGT", (self.run_dir / "TEtrim_complete" / "family.fa").read_text())

    def test_second_mafft_failure_is_not_silently_accepted(self):
        result = self.run_trim("--large", fail=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--large", (self.root / "mafft.log").read_text())
        self.assertFalse((self.run_dir / "TEtrim_mafft" / "mafft_family.fa").exists())
        self.assertFalse((self.run_dir / "TEtrim_complete" / "family.fa").exists())


class BlastPreparationTests(unittest.TestCase):
    def setUp(self):
        try:
            import Bio  # noqa: F401
            import pandas  # noqa: F401
            import pyranges  # noqa: F401
        except ImportError as error:
            self.skipTest(str(error))
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for directory in ("run_0/og", "run_1/raw", "run_1/initial_blast",
                          "run_1/self_search", "run_1/to_align", "run_1/TEtrim_complete"):
            (self.root / directory).mkdir(parents=True)
        for directory in ("run_0/og", "run_1/raw"):
            (self.root / directory / "family.fa").write_text(">family\n" + "ACGT" * 25 + "\n")
        (self.root / "run_1/raw/family.fa.trf").write_text("@family\n")
        (self.root / "genome.fa").write_text(">chr1\n" + "ACGT" * 2500 + "\n")
        self.table = self.root / "hits.tsv"
        line = "family\tchr1\t90\t100\t1\t100\t100\t{}\t{}\t10000\t1e-20\t{}\t100\n"
        with self.table.open("w") as handle:
            handle.write(line.format(1000, 1099, 1000))
            for _ in range(50000):
                handle.write(line.format(2000, 2099, 10))
            handle.write(line.format(5000, 5099, 1000))
        bin_dir = self.root / "bin"
        bin_dir.mkdir()
        blast = bin_dir / "blastn"
        blast.write_text(
            "#!/usr/bin/env python3\n"
            "import os,shutil,sys\n"
            "if os.environ.get('FAKE_BLAST_FAIL'): sys.exit(2)\n"
            "output=sys.argv[sys.argv.index('-out')+1]\n"
            "if '-db' in sys.argv: shutil.copyfile(os.environ['FAKE_BLAST_TABLE'], output)\n"
            "else: open(output,'w').write('a\\t90\\t100\\nb\\t90\\t100\\n')\n"
        )
        blast.chmod(0o755)
        self.env = os.environ.copy()
        self.env["PATH"] = str(bin_dir) + os.pathsep + self.env["PATH"]
        self.env["FAKE_BLAST_TABLE"] = str(self.table)
        self.env["OPENBLAS_NUM_THREADS"] = "1"

    def run_setup(self):
        return subprocess.run(
            [sys.executable, str(SCRIPTS / "initial_mafft_setup.py"),
             "-d", str(self.root), "-r", "1", "-s", "family.fa",
             "-g", str(self.root / "genome.fa"), "-f", "0", "-n", "2", "-D"],
            env=self.env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )

    def test_selects_top_hits_across_chunks_with_stable_ties(self):
        result = self.run_setup()
        self.assertEqual(result.returncode, 0, result.stderr)
        hits = (self.root / "run_1/self_search/family.fa_check_1").read_text()
        self.assertIn("chr1:999-1099", hits)
        self.assertIn("chr1:4999-5099", hits)

    def test_blast_failure_is_not_treated_as_no_hits(self):
        self.env["FAKE_BLAST_FAIL"] = "1"
        result = self.run_setup()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.root / "run_1/TEtrim_complete/family.fa").exists())


@unittest.skipUnless(os.environ.get("TEST_TESTRAINER_INTEGRATION"), "integration test is opt-in")
class PipelineIntegrationTests(unittest.TestCase):
    def test_fast_mode_runs_both_alignments_and_retains_failed_round(self):
        self.run_alignment_mode("fast")

    def test_localpair_mode_runs_both_alignments_and_retains_failed_round(self):
        self.run_alignment_mode("localpair")

    def run_alignment_mode(self, mode):
        required = ("mafft", "parallel", "blastn", "makeblastdb", "trf", "sa-ssr", "mreps", "Rscript")
        for program in required:
            if not shutil.which(program):
                self.skipTest(f"{program} is unavailable")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rng = random.Random(41)
            consensus = "".join(rng.choices("ACGT", k=800))
            genome = []
            for _ in range(8):
                genome.append("".join(rng.choices("ACGT", k=500)))
                genome.append("".join(c if rng.random() > .02 else rng.choice("ACGT") for c in consensus))
            (root / "genome.fa").write_text(">chr1\n" + "".join(genome) + "\n")
            (root / "library.fa").write_text(">family#Unknown\n" + consensus + "\n")
            bin_dir = root / "bin"
            bin_dir.mkdir()
            real_mafft = shutil.which("mafft")
            mafft = bin_dir / "mafft"
            mafft.write_text(
                "#!/usr/bin/env python3\n"
                "import os,sys\n"
                "with open(os.environ['FAKE_MAFFT_LOG'],'a') as log: log.write(' '.join(sys.argv[1:])+'\\n')\n"
                f"os.execv({real_mafft!r}, [{real_mafft!r}]+sys.argv[1:])\n"
            )
            mafft.chmod(0o755)
            classifier = bin_dir / "RepeatClassifier"
            classifier.write_text("#!/bin/sh\nexit 9\n")
            classifier.chmod(0o755)
            env = os.environ.copy()
            env["PATH"] = str(bin_dir) + os.pathsep + env["PATH"]
            env["FAKE_MAFFT_LOG"] = str(root / "mafft.log")
            env.pop("MAFFT_TMPDIR", None)
            wrapper = root / "wrapper.sh"
            source = (ROOT / "scripts/TEstrainer/TEstrainer_for_earlGrey.sh").read_text()
            source = source.replace(
                "STRAIN_SCRIPTS=/data/toby/EarlGrey/scripts/TEstrainer/scripts/",
                f"STRAIN_SCRIPTS={SCRIPTS}/",
            )
            wrapper.write_text(source)
            command = ["bash", str(wrapper), "-g", str(root / "genome.fa"),
                       "-l", str(root / "library.fa"), "-d", str(root / "output"),
                       "-t", "2", "-f", "100", "-r", "1", "-B", "4", "-J", "1"]
            result = subprocess.run(command + ["-A", mode, "-q"], env=env,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    text=True, timeout=120)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("RepeatClassifier exited with status 9", result.stderr)
            calls = (root / "mafft.log").read_text().splitlines()
            self.assertEqual(len(calls), 2, result.stdout + result.stderr)
            if mode == "fast":
                self.assertTrue(all("--6merpair" in call and "--memsave" in call for call in calls))
            else:
                self.assertTrue(all("--localpair" in call and "--large" in call for call in calls))
                run_dir = root / "output/run_1"
                for input_path, output_path, direction in (
                    (run_dir / "to_align/family.fasta", run_dir / "mafft/family.fasta", True),
                    (run_dir / "TEtrim_unaln/unaln_family.fasta", run_dir / "TEtrim_mafft/mafft_family.fasta", False),
                ):
                    comparison_command = [real_mafft, "--quiet", "--thread", "1", "--localpair"]
                    if direction:
                        comparison_command.append("--adjustdirectionaccurately")
                    comparison = subprocess.run(comparison_command + [str(input_path)], env=env,
                                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
                    self.assertEqual(comparison.stdout, output_path.read_bytes())
            self.assertTrue((root / "output/run_1/to_align/family.fasta").exists())
            self.assertTrue((root / "output/run_1/mafft/family.fasta").exists())
            self.assertFalse((root / "output/library.fa.strained").exists())
            wrong_mode = "fast" if mode == "localpair" else "localpair"
            mismatch = subprocess.run(command + ["-A", wrong_mode, "-R", "-q"], env=env,
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                      text=True, timeout=120)
            self.assertNotEqual(mismatch.returncode, 0)
            self.assertIn("cannot resume", mismatch.stderr)
            self.assertTrue((root / "output/run_1/.complete").exists())
            resumed = subprocess.run(command + ["-A", mode, "-R", "-q"], env=env,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                     text=True, timeout=120)
            self.assertNotEqual(resumed.returncode, 0)
            self.assertIn("RepeatClassifier exited with status 9", resumed.stderr)
            self.assertEqual(len((root / "mafft.log").read_text().splitlines()), 2)
            if mode == "fast":
                classifier.write_text("#!/bin/sh\nexit 0\n")
                missing = subprocess.run(command + ["-A", mode, "-R", "-q"], env=env,
                                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                         text=True, timeout=120)
                self.assertNotEqual(missing.returncode, 0)
                self.assertIn("did not produce a classified library", missing.stderr)
                self.assertFalse((root / "output/library.fa.strained").exists())


if __name__ == "__main__":
    unittest.main()
