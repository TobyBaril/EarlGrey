# Development notes for the pangenome pipeline of Earlgrey

### February 2026
We decided to get rid of the EarlGrey set up steps, such a the DFAM configuration etc, in the pangenome pipeline. Before using the pangenome pipeline, a user would have to run the normal installation and set up script. So there is no need to repeat the steps here. The only thing left to do in the pangenome pipeline is to carry out the checks. 

The pipeline expects a config.yaml file called `config/config.yaml`. This file contains the necessary information to run the pipeline such as the species and their genome files. It can also contain some additional options, but for some defaults will be applied (see the validate_parameters function for more details).

I am testing with snakemake version 9.9.0. 

`snakemake --cores 1 --dag ` will run the parameter checks and create the dag. 

Now copy all the text which is not part of the check messages into a `temp.txt` file and use the following command to create the DAG plot. 

`cat temp.txt | dot -Tsvg > dag.svg` 

If you want a diagram that represents a summary of the DAG, independently of the different samples, you can use `--rulegraph` instead. 

`snakemake --cores 1` to run the pipeline.