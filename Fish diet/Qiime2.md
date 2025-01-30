## An amplicon Data Analysis Pipeline

1. Import sequence
``` bash 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.txt \
  --output-path sequences.qza \
  --input-format PairedEndFastqManifestPhred33V2
```

2. View the quality of sequences
``` bash 
qiime demux summarize --i-data sequences.qza \
  --o-visualization seq-qc.qzv
```
3. Merging Paired-End Reads
``` bash 
qiime vsearch merge-pairs \
  --i-demultiplexed-seqs sequences.qza \
  --o-merged-sequences  \
  --o-unmerged-sequences sequences-unjoined.qza
```
5. Denoising reads
``` bash 
qiime deblur denoise-16S \
  --i-demultiplexed-seqs sequences-joined.qza \
  --p-trim-length 470 \
  --p-sample-stats \
  --o-representative-sequences Results/denoise/rep-sequences.qza \
  --o-table Results/denoise/table.qza \
  --o-stats Results/denoise/denoising-stats.qza
```

``` bash 
qiime feature-table summarize \
  --i-data Results/denoise/table.qza \
  --o-visualization Results/denoise/table.qzv
qiime feature-table tabulate-seqs \
  --i-data Results/denoise/rep-sequences.qza \
  --o-visualization Results/denoise/rep-sequences.qzv
```

6. De-replicating Reads
 >   Dereplication is the process of collapsing identical sequences into a single representative sequence and *recording their abundances*. This is a common step in amplicon sequencing analysis to reduce redundancy and simplify downstream analyses. 
 >   
 >   *Denoising tools like **Deblur** and **DADA2** inherently handle dereplication as part of their error-correction process. If you're using these tools, **you do not need to perform dereplication separately**. These tools will:*
 >   
``` bash 
qiime vsearch dereplicate-sequences \
  --i-sequences Results/denoise/rep-sequences.qza \
  --o-dereplicated-table Results/dereplicated/table.qza \
  --o-dereplicated-sequences Results/dereplicated/rep-sequences.qza
```

``` 
```
7. Generate a tree for phylogenetic diversity analyses
``` bash 
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences Results/denoise/rep-sequences.qza \
  --o-alignment Results/tree/aligned-rep-seqs.qza \
  --o-masked-alignment Results/tree/masked-aligned-rep-seqs.qza \
  --o-tree Results/tree/unrooted-tree.qza \
  --o-rooted-tree Results/tree/rooted-tree.qza
```
