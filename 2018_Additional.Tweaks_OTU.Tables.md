# Adjustments to Minimum OTU Size prior to Diversity Assessment and Differential Abundance Testing in R.

### January 8, 2018

Reforming OTU tables based on a minimum OTU size of 1, to allow for singletons and doubletons to be present for alpha diversity testing. Singletons are often required for more accurate diversity measurements and assessments of coverage, although minimum cleaning steps will be done in R using Phyloseq to mitigate as many erroneous reads as possible. Basis of diversity on Shannon or Simpson metrics is also less prone to error from sequencing or PCR bias (Jost, 2007; Marcon, 2015, several papers from Anne Chao).

Also reforming OTU tables in a much more stringent fashion to reduce noise, which may help refine differential abundance testing later, and to ensure that any OTUs undergoing DA testing are true. I am filtering at a level of 0.005% minimum OTU threshold (0.005% of total read count, assuming approximate total of 500,000 or 900,000 reads depending on dataset; this leads to minimum OTU size of 25 or 45 reads, respectively). Filtering threshold is based on quality assessments with Bokulich et al, 2011 paper. Following the production of OTU tables, .biom files will be imported into R to undergo normalization and DA testing with various methods.

#### Script for DSS Feces/Contents datasets, to obtain ALL OTUs, including singletons and doubletons.

```
#!/bin/bash
##

#SBATCH -o /share/magalab/Kat/DSS/pick.otus/all.otus.out
#SBATCH -e /share/magalab/Kat/DSS/pick.otus/all.otus.err
#SBATCH -n 4
#SBATCH -p gc
#SBATCH -t 8-00:00:00

echo "Setting up output and error directories, running 4 ntasks, partitioned to gc computer."
echo "This script is intended to assign OTUs and include singletons and doubletons for DSS Feces and DSS Contents datasets."

echo "Loading QIIME"

module load qiime/1.9.1

echo "Picking OTUs for feces and content samples using open reference method."

time pick_open_reference_otus.py -i /share/magalab/Kat/DSS/chimera.detection/feces.screened.nochimeras.good.fna -o /share/magalab/Kat/DSS/pick.otus/feces_otus_fixed_all/ -a -O 8 --min_otu_size 1 -r /share/magalab/Kat/DSS/qiime/gg_13_8_otus/rep_set/97_otus.fasta

time pick_open_reference_otus.py -i /share/magalab/Kat/DSS/chimera.detection/contents.screened.nochimeras.good.fna -o /share/magalab/Kat/DSS/pick.otus/contents_otus_fixed_all/ -a -O 8 --min_otu_size 1 -r /share/magalab/Kat/DSS/qiime/gg_13_8_otus/rep_set/97_otus.fasta

echo "Done picking OTUs."

echo "Summarizing the number of OTUs per sample"

time biom summarize-table -i /share/magalab/Kat/DSS/pick.otus/feces_otus_fixed_all/otu_table_mc1_w_tax_no_pynast_failures.biom -o /share/magalab/Kat/DSS/pick.otus/feces_otus_fixed_all/results_summary.txt

time biom summarize-table -i /share/magalab/Kat/DSS/pick.otus/contents_otus_fixed_all/otu_table_mc1_w_tax_no_pynast_failures.biom -o /share/magalab/Kat/DSS/pick.otus/contents_otus_fixed_all/results_summary.txt

echo "Done summarizing."
```

#### Script for DSS Feces/Contents datasets, to obtain a STRINGENT OTU table, using an approximately 0.005% minimum OTU size.

```
#!/bin/bash
##

#SBATCH -o /share/magalab/Kat/DSS/pick.otus/stringent.otus.out
#SBATCH -e /share/magalab/Kat/DSS/pick.otus/stringent.otus.err
#SBATCH -n 4
#SBATCH -p gc
#SBATCH -t 8-00:00:00

echo "Setting up output and error directories, running 4 ntasks, partitioned to gc computer."
echo "This script is intended to assign OTUs much more stringently, using a 0.005% of total reads threshold for minimum OTU size."
echo "For datasets with approximately 500,000 reads, the threshold is 25 reads to call an OTU. For the dataset with approximately 900,000 reads, the threshold is 45 reads."

echo "Loading QIIME"

module load qiime/1.9.1

echo "Picking OTUs for feces and content samples using open reference method."

time pick_open_reference_otus.py -i /share/magalab/Kat/DSS/chimera.detection/feces.screened.nochimeras.good.fna -o /share/magalab/Kat/DSS/pick.otus/feces_otus_fixed_stringent/ -a -O 8 --min_otu_size 45 -r /share/magalab/Kat/DSS/qiime/gg_13_8_otus/rep_set/97_otus.fasta

time pick_open_reference_otus.py -i /share/magalab/Kat/DSS/chimera.detection/contents.screened.nochimeras.good.fna -o /share/magalab/Kat/DSS/pick.otus/contents_otus_fixed_stringent/ -a -O 8 --min_otu_size 25 -r /share/magalab/Kat/DSS/qiime/gg_13_8_otus/rep_set/97_otus.fasta

echo "Done picking OTUs."

echo "Summarizing the number of OTUs per sample"

time biom summarize-table -i /share/magalab/Kat/DSS/pick.otus/feces_otus_fixed_stringent/otu_table_mc45_w_tax_no_pynast_failures.biom -o /share/magalab/Kat/DSS/pick.otus/feces_otus_fixed_stringent/results_summary.txt

time biom summarize-table -i /share/magalab/Kat/DSS/pick.otus/contents_otus_fixed_stringent/otu_table_mc25_w_tax_no_pynast_failures.biom -o /share/magalab/Kat/DSS/pick.otus/contents_otus_fixed_stringent/results_summary.txt

echo "Done summarizing."
```
