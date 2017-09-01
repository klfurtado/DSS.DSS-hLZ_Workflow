# DSS-hLZ Workflow
#### A clean, condensed report of all the scripts that worked for the DSS pilot study data, adapted here for DSS-hLZ study data

### 1. Assess Quality of Raw Reads with FASTQC

Ran FASTQC on forward reads only to try and find where barcodes might be placed. "." denotes current working directory.

```
module load fastqc/v0.11.2
fastqc DSS-hLZ_S1_L001_R1_001.fastq -o ./fastqc
```

Move HTML file outside the server (in this case, to my laptop) for viewing:

```
scp -r klfurtado@cabernet.genomecenter.ucdavis.edu:/share/magalab/Kat/DSS-hLZ/fastqc/DSS-hLZ_S1_L001_R1_001_fastqc.html //Users/kathleenfurtado/Documents/maga.lab/fastqc
```

### 2. Run Trimmomatic in Paired End (PE) mode to remove Illumina adapters and primers.

Script:
```
#!/bin/bash
##

#SBATCH -o /share/magalab/Kat/DSS-hLZ/trim/trimmomatic.pe.out
#SBATCH -e /share/magalab/Kat/DSS-hLZ/trim/trimmomatic.pe.err
#SBATCH -n 4
#SBATCH -p gc
#SBATCH -t 8-00:00:00
#SBATCH --mem-per-cpu=10000

echo "Loading Trimmomatic."

module load trimmomatic/0.33

echo "Running Trimmomatic in Paired End mode, using palindromic clipping to remove short contaminating adapter sequences."

trimmomatic PE -threads 4 /share/magalab/Kat/DSS-hLZ/DSS-hLZ_S1_L001_R1_001.fastq /share/magalab/Kat/DSS-hLZ/DSS-hLZ_S1_L001_R2_001.fastq ./dsshlz.trimmed.fwd.paired.fastq ./dsshlz.trimmed.fwd.unpaired.fastq ./dsshlz.trimmed.rev.paired.fastq ./dsshlz.trimmed.rev.unpaired.fastq ILLUMINACLIP:adapters.fa:1:30:10:2:keepBothReads TRAILING:25 SLIDINGWINDOW:5:25 MINLEN:50

echo "Check error and out files."
```

Output:
```
Module trimmomatic-0.33-static loaded.
TrimmomaticPE: Started with arguments: -threads 4 /share/magalab/Kat/DSS-hLZ/DSS-hLZ_S1_L001_R1_001.fastq /share/magalab/Kat/DSS-hLZ/DSS-hLZ_S1_L001_R2_001.fastq ./dsshlz.trimmed.fwd.paired.fastq ./dsshlz.trimmed.fwd.unpaired.fastq ./dsshlz.trimmed.rev.paired.fastq ./dsshlz.trimmed.rev.unpaired.fastq ILLUMINACLIP:adapters.fa:1:30:10:2:keepBothReads TRAILING:25 SLIDINGWINDOW:5:25 MINLEN:50
Using PrefixPair: 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
Using PrefixPair: 'CACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'CACTCTTTCCCTACACGACGCTCTTCCGATCT'
Using PrefixPair: 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC' and 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
Using PrefixPair: 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
Using Long Clipping Sequence: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT'
Using Long Clipping Sequence: 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
ILLUMINACLIP: Using 5 prefix pairs, 4 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Quality encoding detected as phred33
Input Read Pairs: 15729741 Both Surviving: 9972114 (63.40%) Forward Only Surviving: 4457725 (28.34%) Reverse Only Surviving: 318248 (2.02%) Dropped: 981654 (6.24%)
TrimmomaticPE: Completed successfully
```

### 3. Merge Forward and Reverse Reads with FLASH2.

Script:
```
#!/bin/bash
##

#SBATCH -o /share/magalab/Kat/DSS-hLZ/merge/merge.out
#SBATCH -e /share/magalab/Kat/DSS-hLZ/merge/merge.err
#SBATCH -n 4
#SBATCH -p gc
#SBATCH -t 8-00:00:00

echo "Setting up output and error directories, running 4 ntasks, partitioned to animal_sciences computer."

echo "Loading Flash2 module"
module load flash2/c41a82e

echo "Average read length 245 bp (trimmed), average fragment length should be between 290 and 300, and standard deviation of read lengths approximately 29."

echo "Options -r 245, -f 290, and -s 29 are used to calculate max overlap"

time flash2 /share/magalab/Kat/DSS-hLZ/trim/dsshlz.trimmed.fwd.paired.fastq /share/magalab/Kat/DSS-hLZ/trim/dsshlz.trimmed.rev.paired.fastq -o dsshlz.trimmed -r 245 -f 290 -s 29

echo "Merging complete"
```
**Note: script was actually partitioned to gc, not animal_sciences.*

Output:
```
[FLASH] Read combination statistics:
[FLASH]     Total pairs:       9972114
[FLASH]     Discarded pairs:   0
[FLASH]     Percent Discarded: 0.00%
[FLASH]     Combined pairs:    8691526
[FLASH]     Uncombined pairs:  1280588
[FLASH]     Percent combined:  87.16%
[FLASH]
[FLASH] Writing histogram files.
[FLASH]
[FLASH] FLASH v2.2.00 complete!
[FLASH] 1222.969 seconds elapsed
Merging complete
```

### 4. Attach any unmerged forward reads (any reads with potential barcodes) to the merged reads for demultiplexing:

```
cat /share/magalab/Kat/DSS-hLZ/trim/dsshlz.trimmed.fwd.unpaired.fastq /share/magalab/Kat/DSS-hLZ/merge/dsshlz.trimmed.extendedFrags.fastq /share/magalab/Kat/DSS-hLZ/merge/dsshlz.trimmed.notCombined_1.fastq > ./dsshlz.trimmed.merged.fwd.fastq
```

### 5. Validate mapping files before demultiplexing

Script:
```
#!/bin/bash
##

#SBATCH -o /share/magalab/Kat/DSS-hLZ/mapping.files/validate.out
#SBATCH -e /share/magalab/Kat/DSS-hLZ/mapping.files/validate.err
#SBATCH -n 4
#SBATCH -p gc
#SBATCH -t 8-00:00:00

echo "Setting up output and error directories, running 4 ntasks, partitioned to animal_sciences computer."

echo "Loading QIIME"

module load qiime/1.9.1

echo "Validating Mapping Files"

time validate_mapping_file.py -m /share/magalab/Kat/DSS-hLZ/mapping.files/MappingFiles_DSS-hLZ_Feces.txt -o /share/magalab/Kat/DSS-hLZ/mapping.files/validate.mapping

time validate_mapping_file.py -m /share/magalab/Kat/DSS-hLZ/mapping.files/MappingFiles_DSS-hLZ_Contents.txt -o /share/magalab/Kat/DSS-hLZ/mapping.files/validate.mapping

echo "Check log files"
```
**Note: script was actually partitioned to gc, not animal_sciences.*
Check .out file to see if there were any issues in the Mapping File.

### 6. Extract barcodes and Demultiplex

Script:
```
#!/bin/bash
##

#SBATCH -o /share/magalab/Kat/DSS-hLZ/demultiplex/demultiplex.out
#SBATCH -e /share/magalab/Kat/DSS-hLZ/demultiplex/demultiplex.err
#SBATCH -n 4
#SBATCH -p gc
#SBATCH -t 8-00:00:00
#SBATCH --mem-per-cpu=10000

echo "Setting up output and error directories, running 4 ntasks, partitioned to animal_sciences computer."

echo "Loading QIIME"

module load qiime/1.9.1

echo "Running extract_barcodes.py in order to use split_libraries_fastq.py."

time extract_barcodes.py -f /share/magalab/Kat/DSS-hLZ/merge/dsshlz.trimmed.merged.fwd.fastq -o /share/magalab/Kat/DSS-hLZ/demultiplex/contents_barcodes -c barcode_single_end -l 8 -m /share/magalab/Kat/DSS-hLZ/mapping.files/MappingFiles_DSS-hLZ_Contents.txt -a

time extract_barcodes.py -f /share/magalab/Kat/DSS-hLZ/merge/dsshlz.trimmed.merged.fwd.fastq -o /share/magalab/Kat/DSS-hLZ/demultiplex/feces_barcodes -c barcode_single_end -l 8 -m /share/magalab/Kat/DSS-hLZ/mapping.files/MappingFiles_DSS-hLZ_Feces.txt -a

echo "Demultiplexing with split_libraries_fastq.py."

time split_libraries_fastq.py -i /share/magalab/Kat/DSS-hLZ/merge/dsshlz.trimmed.merged.fwd.fastq -o /share/magalab/Kat/DSS-hLZ/demultiplex/splitlibfastq_contents_out/ -m /share/magalab/Kat/DSS-hLZ/mapping.files/MappingFiles_DSS-hLZ_Contents.txt -b /share/magalab/Kat/DSS-hLZ/demultiplex/contents_barcodes/barcodes.fastq -p 0.5 -r 5 --store_qual_scores --store_demultiplexed_fastq -q 24 --barcode_type 8 --phred_offset 33

time split_libraries_fastq.py -i /share/magalab/Kat/DSS-hLZ/merge/dsshlz.trimmed.merged.fwd.fastq -o /share/magalab/Kat/DSS-hLZ/demultiplex/splitlibfastq_feces_out/ -m /share/magalab/Kat/DSS-hLZ/mapping.files/MappingFiles_DSS-hLZ_Feces.txt -b /share/magalab/Kat/DSS-hLZ/demultiplex/feces_barcodes/barcodes.fastq -p 0.5 -r 5 --store_qual_scores --store_demultiplexed_fastq -q 24 --barcode_type 8 --phred_offset 33

echo "Check .err and .out files"
```

**Note: script was actually partitioned to gc, not animal_sciences.*

Output from demultiplexing:

```
klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/demultiplex/splitlibfastq_feces_out$ cat split_library_log.txt
Input file paths
Mapping filepath: /share/magalab/Kat/DSS-hLZ/mapping.files/MappingFiles_DSS-hLZ_Feces.txt (md5: c041c0bf6b3969cb46a775dd3229f35d)
Sequence read filepath: /share/magalab/Kat/DSS-hLZ/merge/dsshlz.trimmed.merged.fwd.fastq (md5: ee894ccb65d842fd7da1d70502cd158c)
Barcode read filepath: /share/magalab/Kat/DSS-hLZ/demultiplex/feces_barcodes/barcodes.fastq (md5: 64fa124004cdae03063b8c3596035855)

Quality filter results
Total number of input sequences: 14429839
Barcode not in mapping file: 13796159
Read too short after quality truncation: 0
Count of N characters exceeds limit: 24
Illumina quality digit = 0: 0
Barcode errors exceed max: 0

Result summary (after quality filtering)
Median sequence length: 301.00
468	25923
458	22423
464	19328
466	17803
467	15923
596	15389
471	15368
591	14375
474	12042
586	11905
481	11683
599	11305
583	10886
498	10810
580	10763
463	10693
607	10564
593	10371
473	10363
585	9825
490	9602
589	9537
588	9182
446	9181
598	9037
450	8960
485	8925
491	8885
454	8798
582	8781
483	8714
503	8447
477	8348
594	8275
452	8177
455	8076
489	8027
465	7881
502	7817
496	7676
595	7659
501	7654
443	7646
495	7593
494	7507
499	7457
470	7318
587	7269
456	7238
486	7236
447	7133
497	7080
479	6843
451	6739
581	5783
442	5722
493	5717
505	5619
453	5603
484	5463
441	5422
597	5399
590	5293
592	5277
605	3460
584	3448
600	3406
492	3383
500	3311
444	3276
445	2660
604	995
602	989
603	626
488	155
601	102
606	72
448	34
449	1

Total number seqs written	633656
```
```
klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/demultiplex/splitlibfastq_contents_out$ cat split_library_log.txt
Input file paths
Mapping filepath: /share/magalab/Kat/DSS-hLZ/mapping.files/MappingFiles_DSS-hLZ_Contents.txt (md5: b454fb612219beb187e43b987f786b7f)
Sequence read filepath: /share/magalab/Kat/DSS-hLZ/merge/dsshlz.trimmed.merged.fwd.fastq (md5: ee894ccb65d842fd7da1d70502cd158c)
Barcode read filepath: /share/magalab/Kat/DSS-hLZ/demultiplex/contents_barcodes/barcodes.fastq (md5: 64fa124004cdae03063b8c3596035855)

Quality filter results
Total number of input sequences: 14429839
Barcode not in mapping file: 13769538
Read too short after quality truncation: 0
Count of N characters exceeds limit: 27
Illumina quality digit = 0: 0
Barcode errors exceed max: 0

Result summary (after quality filtering)
Median sequence length: 301.00
554	14557
591	14375
544	12973
578	12880
543	12467
559	12221
552	12144
562	11949
586	11905
550	11610
522	11440
551	11274
527	11055
583	10886
580	10763
519	10648
535	10579
546	10507
542	10448
530	10411
537	10392
593	10371
558	10228
567	10035
560	10034
526	10034
579	9942
538	9935
585	9825
557	9582
589	9537
548	9337
525	9274
588	9182
533	9009
549	8882
532	8795
582	8781
528	8691
566	8621
529	8603
531	8442
534	8441
541	8413
523	8402
545	8316
594	8275
565	8175
553	8154
513	8136
540	7920
524	7765
595	7659
561	7591
564	7291
587	7269
556	6878
512	6635
547	6067
511	5824
581	5783
563	5597
577	5378
590	5293
592	5277
539	5195
521	4827
555	4785
520	4766
506	4745
515	4597
518	4005
508	3829
514	3797
584	3448
509	2471
507	2238
517	2038
516	1885
510	485

Total number seqs written	660274
```
Use output .fastq file to run Fastq_screen.

### 7. Run Fastq_screen to remove reads from PhiX or Pig Mitochondria

Script:
```
#!/bin/bash

#SBATCH -o /share/magalab/Kat/DSS/index.screen/dsshlz.fastq.screen.out
#SBATCH -e /share/magalab/Kat/DSS/index.screen/dsshlz.fastq.screen.err
#SBATCH -n 4
#SBATCH -p gc
#SBATCH -t 8-00:00:00

echo "Setting up output and error directories, running 4 ntasks."

echo "Loading bowtie2 v. 2.2.8"

module load bowtie2/2.2.8

echo "Running fastq_screen to check merged and trimmed reads for hits within PhiX and Sus Scrofa mitochondrial genome."

time fastq_screen_v0.9.5/fastq_screen --threads 8 --subset 0 --aligner bowtie2 --force --nohits /share/magalab/Kat/DSS-hLZ/demultiplex/splitlibfastq_feces_out/feces.fastq

time fastq_screen_v0.9.5/fastq_screen --threads 8 --subset 0 --aligner bowtie2 --force --nohits /share/magalab/Kat/DSS-hLZ/demultiplex/splitlibfastq_contents_out/contents.fastq

echo "Reads should be screened, check output."
echo "Run FastQC!"
```

Output:
```
klfurtado@cabernet:/share/magalab/Kat/DSS/index.screen$ cat feces_screen.txt
#Fastq_screen version: 0.9.5	#Aligner: bowtie2	#Processing all reads in FASTQ files
Genome	#Reads_processed	#Unmapped	%Unmapped	#One_hit_one_genome	%One_hit_one_genome	#Multiple_hits_one_genome	%Multiple_hits_one_genome	#One_hit_multiple_genomes	%One_hit_multiple_genomes	Multiple_hits_multiple_genomes	%Multiple_hits_multiple_genomes
pig.mt	633656	633656	100.00	0	0.00	0	0.00	0	0.00	0	0.00
phi.x.174	633656	633141	99.92	512	0.08	3	0.00	0	0.00	0	0.00

%Hit_no_genomes: 99.92
```
```
#Fastq_screen version: 0.9.5	#Aligner: bowtie2	#Processing all reads in FASTQ files
Genome	#Reads_processed	#Unmapped	%Unmapped	#One_hit_one_genome	%One_hit_one_genome	#Multiple_hits_one_genome	%Multiple_hits_one_genome	#One_hit_multiple_genomes	%One_hit_multiple_genomes	Multiple_hits_multiple_genomes	%Multiple_hits_multiple_genomes
pig.mt	660274	660258	100.00	16	0.00	0	0.00	0	0.00	0	0.00
phi.x.174	660274	659150	99.83	1063	0.16	61	0.01	0	0.00	0	0.00

%Hit_no_genomes: 99.83
```

### 8. Convert back to FASTA

Script:
```
klfurtado@cabernet:/share/magalab/Kat/DSS/index.screen$ cat dsshlz.convert2fasta.sh
#!/bin/bash
##

#SBATCH -o /share/magalab/Kat/DSS/index.screen/dsshlz.convert.out
#SBATCH -e /share/magalab/Kat/DSS/index.screen/dsshlz.convert.err
#SBATCH -n 4
#SBATCH -p gc
#SBATCH -t 8-00:00:00

echo "Setting up output and error directories, running 4 ntasks."

echo "Loading QIIME"

module load qiime/1.9.1

echo "Converting screened FASTQ files to FASTA+QUAL for chimera detection."

time convert_fastaqual_fastq.py -c fastq_to_fastaqual -f /share/magalab/Kat/DSS/index.screen/dsshlz.feces.screened.fastq -o /share/magalab/Kat/DSS-hLZ/screened/ -F

time convert_fastaqual_fastq.py -c fastq_to_fastaqual -f /share/magalab/Kat/DSS/index.screen/dsshlz.contents.screened.fastq -o /share/magalab/Kat/DSS-hLZ/screened/ -F

echo "Done converting."
```

### 9. Remove chimeras

Script:
```
klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/chimera.detection$ cat chimera.detection.sh
#!/bin/bash
##

#SBATCH -o /share/magalab/Kat/DSS-hLZ/chimera.detection/chimeras.out
#SBATCH -e /share/magalab/Kat/DSS-hLZ/chimera.detection/chimeras.err
#SBATCH -n 4
#SBATCH -p gc
#SBATCH -t 8-00:00:00
#SBATCH --mem-per-cpu=10000

echo "Setting up output and error directories, running 4 ntasks, partitioned to animal_sciences computer."

echo "Loading QIIME"

module load qiime/1.9.1

echo "Checking for chimeras in feces and contents samples using usearch61 method with default parameters"

time identify_chimeric_seqs.py -i /share/magalab/Kat/DSS-hLZ/screened/dsshlz.contents.screened.fna -r /share/magalab/Kat/DSS/qiime/gg_13_8_otus/rep_set/97_otus.fasta -m usearch61 --usearch61_minh 0.28 --usearch61_xn 8.0 -o /share/magalab/Kat/DSS-hLZ/chimera.detection/contents_chimera_output

time identify_chimeric_seqs.py -i /share/magalab/Kat/DSS-hLZ/screened/dsshlz.feces.screened.fna -r /share/magalab/Kat/DSS/qiime/gg_13_8_otus/rep_set/97_otus.fasta -m usearch61 --usearch61_minh 0.28 --usearch61_xn 8.0 -o /share/magalab/Kat/DSS-hLZ/chimera.detection/feces_chimera_output

echo "Filtering chimeric sequences."

time filter_fasta.py  -f /share/magalab/Kat/DSS-hLZ/screened/dsshlz.contents.screened.fna -o /share/magalab/Kat/DSS-hLZ/chimera.detection/contents.screened.nochimeras.fna -s /share/magalab/Kat/DSS-hLZ/chimera.detection/contents_chimera_output/chimeras.txt -n

time filter_fasta.py  -f /share/magalab/Kat/DSS-hLZ/screened/dsshlz.feces.screened.fna -o /share/magalab/Kat/DSS-hLZ/chimera.detection/feces.screened.nochimeras.fna -s /share/magalab/Kat/DSS-hLZ/chimera.detection/feces_chimera_output/chimeras.txt -n

echo "Done Filtering."
```

Output:

DSS-hLZ Feces
```
klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/chimera.detection/feces_chimera_output$ cat identify_chimeric_seqs.log
input_seqs_fp	/share/magalab/Kat/DSS-hLZ/screened/dsshlz.feces.screened.fna
output_dir	/share/magalab/Kat/DSS-hLZ/chimera.detection/feces_chimera_output
reference_seqs_fp	/share/magalab/Kat/DSS/qiime/gg_13_8_otus/rep_set/97_otus.fasta
suppress_usearch61_intermediates	False
suppress_usearch61_ref	False
suppress_usearch61_denovo	False
split_by_sampleid	False
non_chimeras_retention	union
usearch61_minh	0.28
usearch61_xn	8.0
usearch61_dn	1.4
usearch61_mindiffs	3
usearch61_mindiv	0.8
usearch61_abundance_skew	2.0
percent_id_usearch61	0.97
minlen	64
word_length	8
max_accepts	1
max_rejects	8
HALT_EXEC	False

ref_non_chimeras	591083
ref_chimeras	42058
denovo_chimeras	29870
denovo_non_chimeras	603271
```

DSS-hLZ Contents
```
klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/chimera.detection/contents_chimera_output$ ls
chimeras.txt					     dsshlz.contents.screened.fna_consensus_with_abundance.fasta
dsshlz.contents.screened.fna_chimeras_denovo.log     dsshlz.contents.screened.fna_consensus_with_abundance.uc
dsshlz.contents.screened.fna_chimeras_denovo.uchime  dsshlz.contents.screened.fna_smallmem_clustered.log
dsshlz.contents.screened.fna_chimeras_ref.log	     identify_chimeric_seqs.log
dsshlz.contents.screened.fna_chimeras_ref.uchime     non_chimeras.txt
dsshlz.contents.screened.fna_consensus_fixed.fasta
klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/chimera.detection/contents_chimera_output$ cat identify_chimeric_seqs.log
input_seqs_fp	/share/magalab/Kat/DSS-hLZ/screened/dsshlz.contents.screened.fna
output_dir	/share/magalab/Kat/DSS-hLZ/chimera.detection/contents_chimera_output
reference_seqs_fp	/share/magalab/Kat/DSS/qiime/gg_13_8_otus/rep_set/97_otus.fasta
suppress_usearch61_intermediates	False
suppress_usearch61_ref	False
suppress_usearch61_denovo	False
split_by_sampleid	False
non_chimeras_retention	union
usearch61_minh	0.28
usearch61_xn	8.0
usearch61_dn	1.4
usearch61_mindiffs	3
usearch61_mindiv	0.8
usearch61_abundance_skew	2.0
percent_id_usearch61	0.97
minlen	64
word_length	8
max_accepts	1
max_rejects	8
HALT_EXEC	False

ref_non_chimeras	618531
ref_chimeras	40603
denovo_chimeras	15423
denovo_non_chimeras	643711
```

### 10. Trim off reads which are too long, using MOTHUR

Run summary.seqs to determine distribution of read lengths, and use screen.seqs to filter out long reads.

Batch files to run summary.seqs
```
klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/mothur$ cat contents.summary.bat
summary.seqs(fasta=/share/magalab/Kat/DSS-hLZ/chimera.detection/contents.screened.nochimeras.fna)
```
```
klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/mothur$ cat feces.summary.bat
summary.seqs(fasta=/share/magalab/Kat/DSS-hLZ/chimera.detection/feces.screened.nochimeras.fna)
```

Outputs:
```
Batch Mode


mothur > summary.seqs(fasta=/share/magalab/Kat/DSS-hLZ/chimera.detection/feces.screened.nochimeras.fna)

Using 1 processors.

		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	50	50	0	2	1
2.5%-tile:	1	84	84	0	3	15722
25%-tile:	1	198	198	0	4	157219
Median: 	1	301	301	0	4	314438
75%-tile:	1	302	302	0	5	471657
97.5%-tile:	1	302	302	0	5	613154
Maximum:	1	439	439	0	241	628875
Mean:	1	251.533	251.533	0	4.10946
# of Seqs:	628875

Output File Names:
/share/magalab/Kat/DSS-hLZ/chimera.detection/feces.screened.nochimeras.summary

It took 4729 secs to summarize 628875 sequences.

mothur > quit()
```
```
Batch Mode


mothur > summary.seqs(fasta=/share/magalab/Kat/DSS-hLZ/chimera.detection/contents.screened.nochimeras.fna)

Using 1 processors.

		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	50	50	0	2	1
2.5%-tile:	1	84	84	0	3	16425
25%-tile:	1	198	198	0	4	164245
Median: 	1	301	301	0	4	328489
75%-tile:	1	302	302	0	5	492733
97.5%-tile:	1	302	302	0	5	640553
Maximum:	1	476	476	0	231	656977
Mean:	1	251.277	251.277	0	4.1701
# of Seqs:	656977

Output File Names:
/share/magalab/Kat/DSS-hLZ/chimera.detection/contents.screened.nochimeras.summary

It took 6728 secs to summarize 656977 sequences.

mothur > quit()
```

Batch files to run screen.seqs
```
klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/mothur$ cat feces.screen.bat
screen.seqs(fasta=/share/magalab/Kat/DSS-hLZ/chimera.detection/feces.screened.nochimeras.fna, minlength=50, maxlength=350, maxambig=0)
```
```
klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/mothur$ cat contents.screen.bat
screen.seqs(fasta=/share/magalab/Kat/DSS-hLZ/chimera.detection/contents.screened.nochimeras.fna, minlength=50, maxlength=350, maxambig=0)
```

Outputs from screening long reads:
```
Batch Mode


mothur > screen.seqs(fasta=/share/magalab/Kat/DSS-hLZ/chimera.detection/contents.screened.nochimeras.fna, minlength=50, maxlength=350, maxambig=0)

Using 1 processors.

Output File Names:
/share/magalab/Kat/DSS-hLZ/chimera.detection/contents.screened.nochimeras.good.fna
/share/magalab/Kat/DSS-hLZ/chimera.detection/contents.screened.nochimeras.bad.accnos


It took 3074 secs to screen 656977 sequences.

mothur > quit()
```
```
Batch Mode


mothur > screen.seqs(fasta=/share/magalab/Kat/DSS-hLZ/chimera.detection/feces.screened.nochimeras.fna, minlength=50, maxlength=350, maxambig=0)

Using 1 processors.

Output File Names:
/share/magalab/Kat/DSS-hLZ/chimera.detection/feces.screened.nochimeras.good.fna
/share/magalab/Kat/DSS-hLZ/chimera.detection/feces.screened.nochimeras.bad.accnos


It took 802 secs to screen 628875 sequences.

mothur > quit()
```

### 11. Pick OTUs using Open Reference Method
```
klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/pickotus$ cat pickotus.sh
#!/bin/bash
##

#SBATCH -o /share/magalab/Kat/DSS-hLZ/pickotus/pickotus.out
#SBATCH -e /share/magalab/Kat/DSS-hLZ/pickotus/pickotus.err
#SBATCH -n 4
#SBATCH -p gc
#SBATCH -t 8-00:00:00

echo "Setting up output and error directories, running 4 ntasks."

echo "Loading QIIME"

module load qiime/1.9.1

echo "Picking OTUs for feces and content samples using open reference method."

time pick_open_reference_otus.py -i /share/magalab/Kat/DSS-hLZ/chimera.detection/feces.screened.nochimeras.good.fna -o /share/magalab/Kat/DSS-hLZ/pickotus/feces_otus/ -f -a -O 8 --min_otu_size 3 -r /share/magalab/Kat/DSS/qiime/gg_13_8_otus/rep_set/97_otus.fasta

time pick_open_reference_otus.py -i /share/magalab/Kat/DSS-hLZ/chimera.detection/contents.screened.nochimeras.good.fna -o /share/magalab/Kat/DSS-hLZ/pickotus/contents_otus/ -a -O 8 --min_otu_size 3 -r /share/magalab/Kat/DSS/qiime/gg_13_8_otus/rep_set/97_otus.fasta

echo "Done picking OTUs."

echo "Summarizing the number of OTUs per sample"

time biom summarize-table -i /share/magalab/Kat/DSS-hLZ/pickotus/feces_otus/otu_table_mc3_w_tax_no_pynast_failures.biom -o /share/magalab/Kat/DSS-hLZ/pickotus/feces_otus/results_summary.txt

time biom summarize-table -i /share/magalab/Kat/DSS-hLZ/pickotus/contents_otus/otu_table_mc3_w_tax_no_pynast_failures.biom -o /share/magalab/Kat/DSS-hLZ/pickotus/contents_otus/results_summary.txt

echo "Done summarizing."
```
Summary: Feces OTUs
```
klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/pickotus/feces_otus$ cat results_summary.txt
Num samples: 78
Num observations: 7120
Total count: 510303
Table density (fraction of non-zero values): 0.072

Counts/sample summary:
 Min: 0.0
 Max: 21868.0
 Median: 6251.000
 Mean: 6542.346
 Std. dev.: 3913.832
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

Counts/sample detail:
 449: 0.0
 606: 53.0
 488: 64.0
 601: 72.0
 603: 450.0
 602: 689.0
 604: 788.0
 445: 1882.0
 444: 2354.0
 500: 2421.0
 600: 2535.0
 492: 2634.0
 584: 2878.0
 605: 2919.0
 484: 3777.0
 505: 3911.0
 441: 4515.0
 597: 4526.0
 453: 4533.0
 442: 4601.0
 590: 4605.0
 493: 4621.0
 592: 4781.0
 581: 4857.0
 486: 4934.0
 451: 4971.0
 479: 5094.0
 497: 5473.0
 594: 5633.0
 443: 5783.0
 470: 5789.0
 496: 5899.0
 499: 5924.0
 587: 5960.0
 595: 5980.0
 447: 6071.0
 456: 6096.0
 465: 6193.0
 452: 6224.0
 491: 6278.0
 477: 6288.0
 494: 6368.0
 501: 6431.0
 489: 6458.0
 495: 6466.0
 485: 6576.0
 454: 6708.0
 502: 6800.0
 455: 6929.0
 483: 6962.0
 582: 7392.0
 588: 7503.0
 503: 7546.0
 490: 7555.0
 450: 7646.0
 446: 7796.0
 598: 7853.0
 473: 8127.0
 589: 8253.0
 498: 8299.0
 585: 8656.0
 474: 8686.0
 463: 8699.0
 593: 8722.0
 580: 8729.0
 599: 9095.0
 607: 9109.0
 586: 9217.0
 481: 9627.0
 583: 9699.0
 591: 10372.0
 471: 11741.0
 596: 12154.0
 467: 12618.0
 466: 15303.0
 464: 16480.0
 458: 18804.0
 468: 21868.0
 ```

 Summary: Contents OTUs
 ```
 klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/pickotus/contents_otus$ cat results_summary.txt
Num samples: 80
Num observations: 5504
Total count: 523195
Table density (fraction of non-zero values): 0.067

Counts/sample summary:
 Min: 221.0
 Max: 12347.0
 Median: 6902.000
 Mean: 6539.938
 Std. dev.: 2557.887
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

Counts/sample detail:
 510: 221.0
 516: 1221.0
 517: 1645.0
 507: 1761.0
 509: 1817.0
 584: 2544.0
 514: 2903.0
 508: 2909.0
 518: 3098.0
 515: 3431.0
 512: 3433.0
 506: 3550.0
 521: 3705.0
 555: 3847.0
 520: 3901.0
 539: 4085.0
 592: 4102.0
 577: 4114.0
 563: 4329.0
 590: 4339.0
 513: 4445.0
 581: 4710.0
 511: 4716.0
 547: 4864.0
 556: 5569.0
 564: 5575.0
 587: 5656.0
 524: 5703.0
 595: 5779.0
 594: 5779.0
 561: 6104.0
 529: 6140.0
 531: 6277.0
 565: 6436.0
 566: 6624.0
 553: 6691.0
 582: 6701.0
 588: 6772.0
 532: 6785.0
 523: 6877.0
 540: 6927.0
 541: 6932.0
 548: 7009.0
 545: 7085.0
 528: 7090.0
 525: 7212.0
 557: 7226.0
 579: 7540.0
 593: 7601.0
 534: 7640.0
 549: 7689.0
 533: 7892.0
 519: 7895.0
 567: 7984.0
 538: 8023.0
 535: 8095.0
 589: 8099.0
 530: 8238.0
 580: 8267.0
 551: 8268.0
 585: 8275.0
 558: 8280.0
 560: 8413.0
 583: 8678.0
 526: 8906.0
 559: 8911.0
 522: 8942.0
 542: 8965.0
 537: 9009.0
 586: 9087.0
 546: 9106.0
 527: 9277.0
 550: 9799.0
 578: 10002.0
 562: 10068.0
 543: 10252.0
 591: 10737.0
 552: 10773.0
 544: 11498.0
 554: 12347.0
```

### 12. Filter samples with no OTUs (They cannot be entered into Phyloseq later)

```
klfurtado@cabernet:/share/magalab/Kat/DSS-hLZ/pickotus/feces_otus$ cat filter.feces.sh
#!/bin/bash

#SBATCH -o /share/magalab/Kat/DSS-hLZ/pickotus/feces_otus/filter.out
#SBATCH -e /share/magalab/Kat/DSS-hLZ/pickotus/feces_otus/filter.err
#SBATCH -n 4
#SBATCH -p gc
#SBATCH -t 8-00:00:00

echo "Loading QIIME v 1.9.1"

module load qiime/1.9.1

echo "Filtering feces .biom file to remove samples with no OTUs."

time filter_samples_from_otu_table.py -i /share/magalab/Kat/DSS-hLZ/pickotus/feces_otus/otu_table_mc3_w_tax_no_pynast_failures.biom -o /share/magalab/Kat/DSS-hLZ/pickotus/feces_otus/otu_table_mc3_w_tax_no_pynast_failures_filtered.biom -n 1 -m /share/magalab/Kat/DSS-hLZ/mapping.files/MappingFiles_DSS-hLZ_Feces.txt --output_mapping_fp /share/magalab/Kat/DSS-hLZ/mapping.files/MappingFiles_DSS-hLZ_Feces_Filtered.txt

echo "Done Filtering, try importing new .biom file into R."
```


#### Take the final .biom files and transfer them to the local hard drive. Then import into R using Phyloseq, along with metadata and tree files. See R script for details. 
