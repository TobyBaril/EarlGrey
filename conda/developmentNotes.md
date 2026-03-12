# Updating Earl Grey to remove Python 3.9 dependency

Bjorn wants me to update Earl Grey to remove the Python 3.9 dependency. This is because Python 3.9 is no longer supported and we want to ensure that Earl Grey can run on newer versions of Python.

I will build Earl Grey with the meta dependency removed and test it to ensure that it still works correctly. I will also update the documentation to reflect the change in dependencies.

Build the conda package:

```bash
conda build .
```

Make a new conda environment and install the package:

```bash
conda create -n earlgrey_python_test -c conda-forge -c bioconda earlgrey --use-local
conda activate earlgrey_python_test
```

This still builds with 3.9. I will force it to build with 3.10 by adding `python >=3.10` to the `run` dependencies in `meta.yaml`. After making this change, I will rebuild the package and test it again. I also removed the dependency on `ncls` version 0.0.64 and replaced it with `ncls` without a specific version to allow for newer versions of `ncls` to be used.

```bash
conda build .
conda create -n earlgrey_python_test -c conda-forge -c bioconda earlgrey --use-local
conda activate earlgrey_python_test
```

Link the dfam libraries and configure RepeatMasker:

```bash
ln -sf ln -sf /data/toby/tools/earlgrey_databases/Libraries/famdb/* \
    /data/toby/miniforge3/envs/earlgrey_python_test/share/RepeatMasker/Libraries/famdb/

cd /data/toby/miniforge3/envs/earlgrey_python_test/share/RepeatMasker
perl configure
```

Run a test command to ensure that Earl Grey is working correctly:

```bash
cd /data/toby/testDIR/
earlGrey -g /data/toby/tools/test.fasta -s test_python312 -o . -t 16
```

This passed with a small fix to the `backSwap.py` script where I changed the regular expression for splitting the input file from `'\s+'` to `r'\s+'` to ensure that it correctly handles whitespace. After making this change, I will test the script again to ensure that it works correctly with the new dependencies.

```bash
cd /data/toby/testDIR/
earlGrey -g /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1A5.fa -s test_python312_2 -o . -t 32
```

This also passed successfully, confirming that Earl Grey is now compatible with Python 3.10 and does not have a dependency on Python 3.9. I will now update the documentation to reflect these changes and ensure that users are aware of the new dependencies when installing Earl Grey.