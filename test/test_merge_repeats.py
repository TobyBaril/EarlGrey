import csv
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "mergeRepeats.R"


def _has_r_deps():
    if not shutil.which("Rscript"):
        return False
    check = subprocess.run(
        [
            "Rscript",
            "-e",
            'pkgs <- c("tidyverse","plyr","dplyr","magrittr","data.table"); '
            "quit(status = ifelse(all(sapply(pkgs, requireNamespace, quietly = TRUE)), 0, 1))",
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return check.returncode == 0


class MergeRepeatsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not _has_r_deps():
            raise unittest.SkipTest("Rscript or required R packages are unavailable")

    def run_merge_repeats(self, rows, lowend, genome_size=10000):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            in_gff = root / "in.rmerge.gff.filtered"
            in_gff.write_text("\n".join("\t".join(map(str, row)) for row in rows) + "\n")
            merged_bed = root / "merged.bed"
            revised_table = root / "merged.revisedTable"
            filt_bed = root / "filtered.bed"
            filt_sum = root / "filtered.summary"
            result = subprocess.run(
                [
                    "Rscript",
                    str(SCRIPT),
                    str(in_gff),
                    str(merged_bed),
                    str(genome_size),
                    str(revised_table),
                    str(filt_bed),
                    str(filt_sum),
                    lowend,
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with revised_table.open() as handle:
                unfiltered = {row["Family"]: row for row in csv.DictReader(handle, delimiter="\t")}
            with filt_sum.open() as handle:
                filtered = {row["Family"]: row for row in csv.DictReader(handle, delimiter="\t")}
            filt_bed_lines = [line for line in filt_bed.read_text().splitlines() if line.strip()]
            return unfiltered, filtered, filt_bed_lines

    def test_lowend_yes_excludes_short_repeats_from_filtered_summary(self):
        rows = [
            ("chr1", "RepeatMasker", "DNA/TcMar", 1, 200, 100, "+", ".", "ID=1;NAME=fam1"),
            ("chr1", "RepeatMasker", "DNA/TcMar", 301, 380, 50, "+", ".", "ID=2;NAME=fam2"),
            ("chr1", "RepeatMasker", "LTR/Gypsy", 501, 600, 50, "+", ".", "ID=3;NAME=fam3;FULLY_NESTED=TRUE"),
            ("chr1", "RepeatMasker", "LTR/Gypsy", 701, 950, 100, "+", ".", "ID=4;NAME=fam4;FULLY_NESTED=TRUE"),
        ]
        unfiltered, filtered, filt_bed_lines = self.run_merge_repeats(rows, lowend="yes")
        self.assertEqual(len(filt_bed_lines), 2)
        self.assertEqual(unfiltered["DNA"]["Number of Elements"], "2")
        self.assertEqual(unfiltered["DNA"]["Length Occupied"], "280")
        self.assertEqual(unfiltered["LTR-nested"]["Number of Elements"], "2")
        self.assertEqual(unfiltered["LTR-nested"]["Length Occupied"], "350")
        self.assertEqual(filtered["DNA"]["Number of Elements"], "1")
        self.assertEqual(filtered["DNA"]["Length Occupied"], "200")
        self.assertEqual(filtered["Total Interspersed Repeat"]["Number of Elements"], "1")
        self.assertEqual(filtered["Total Interspersed Repeat"]["Length Occupied"], "200")
        self.assertEqual(filtered["Unmasked"]["Length Occupied"], "9800")
        self.assertEqual(filtered["LTR-nested"]["Number of Elements"], "1")
        self.assertEqual(filtered["LTR-nested"]["Length Occupied"], "250")

    def test_lowend_no_retains_short_repeats_in_filtered_summary(self):
        rows = [
            ("chr1", "RepeatMasker", "DNA/TcMar", 1, 200, 100, "+", ".", "ID=1;NAME=fam1"),
            ("chr1", "RepeatMasker", "DNA/TcMar", 301, 380, 50, "+", ".", "ID=2;NAME=fam2"),
        ]
        unfiltered, filtered, filt_bed_lines = self.run_merge_repeats(rows, lowend="no")
        self.assertEqual(len(filt_bed_lines), 2)
        self.assertEqual(filtered["DNA"]["Number of Elements"], unfiltered["DNA"]["Number of Elements"])
        self.assertEqual(filtered["DNA"]["Length Occupied"], "280")


if __name__ == "__main__":
    unittest.main()
