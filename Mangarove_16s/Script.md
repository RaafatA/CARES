Importing Sequencing Data 
``` shell 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.txt \
  --output-path sequences.qza \
  --input-format PairedEndFastqManifestPhred33V2
```

```bash 
qiime vsearch merge-pairs \
  --i-demultiplexed-seqs sequences.qza \
  --o-merged-sequences sequences-joined.qza \
  --o-unmerged-sequences sequences-unjoined.qza
```

``` bash
qiime deblur denoise-16S --i-demultiplexed-seqs sequences.qza --p-trim-length -1 --p-sample-stats --o-representative-sequences denoise/rep-sequences.qza --o-table denoise/table.qza --o-stats denoise/denoising-stats.qza
```

``` bash 
qiime vsearch dereplicate-sequences --i-sequences denoise/rep-sequences.qza --o-dereplicated-table dereplicated/table.qza --o-dereplicated-sequences dereplicated/rep-sequences.qza
```


## 2. GreenGenes preprocessing 
- Extract the V3-V4 Region Using Primers
``` bash
qiime rescript extract-reads \ --i-sequences greengenes-sequences.qza \ --p-f-primer CCTACGGGNGGCWGCAG \ --p-r-primer GACTACHVGGGTATCTAATCC \ --p-trunc-len 250 \ --o-reads greengenes-v3v4-sequences.qza
```

``` bash 
qiime taxa filter-table \
  --i-table greengenes-taxonomy.qza \
  --m-metadata-file greengenes-taxonomy.qza \
  --p-include "Bacteria,Archaea" \
  --o-filtered-table filtered-taxonomy.qza
```

4. “Culling” low-quality sequences with cull-seqs
```bash 
qiime rescript cull-seqs \
    --i-sequences 2022.10.seqs.fna.qza \
    --o-clean-sequences greengenes2-2022.10-seqs-cleaned.qza
```

``` bash 
qiime rescript trim-alignment \ --i-sequences greengenes2-2022.10-seqs-cleaned.qza \ --p-primer-f CCTACGGGNGGCWGCAG \ --p-primer-r GACTACHVGGGTATCTAATCC \ --o-trimmed-sequences trimmed-v3v4-sequences.qza
```
5. Filtering sequences by length and taxonomy
``` bash 
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences greengenes2-2022.10-seqs-cleaned.qza \
    --i-taxonomy 2022.10.taxonomy.id.tsv.qza \
    --p-labels Archaea Bacteria \
    --p-min-lens 900 1200 \
    --o-filtered-seqs greengenes2-2022.10-seqs-filt.qza \
    --o-discarded-seqs greengenes2-2022.10-seqs-discard.qza
```

dereplication 
``` bash
qiime rescript dereplicate \
    --i-sequences greengenes2-2022.10-seqs-cleaned.qza \
    --i-taxa 2022.10.taxonomy.id.tsv.qza \
    --p-perc-identity 0.97 \
    --p-mode 'lca' \
    --o-dereplicated-sequences gg2-c97-seqs-derep-lca.qza \
    --o-dereplicated-taxa gg2-c97-tax-derep-lca.qza
```



## I have downloaded the Full-backbone file

#### 1. “Culling” low-quality sequences with cull-seqs
``` bash 
qiime rescript cull-seqs \
    --i-sequences 2022.10.backbone.full-length.fna.qza \
    --o-clean-sequences References/gg2-backbone.full-length-2022.10-seqs-cleaned.qza --p-n-jobs 4
```

#### 2. 
``` bash 
qiime rescript dereplicate \
    --i-sequences References/gg2-backbone.full-length-2022.10-seqs-cleaned.qza \
    --i-taxa 2022.10.backbone.tax.qza \
    --p-perc-identity 0.97 \
    --p-mode 'lca' \
    --o-dereplicated-sequences References/gg2-c97-seqs-derep-lca.qza \
    --o-dereplicated-taxa References/gg2-c97-tax-derep-lca.qza
```

### 3. Clustering using Vsearch- closed reference :
``` bash
qiime vsearch cluster-features-closed-reference \
  --i-table denoise/table.qza \
  --i-sequences denoise/rep-sequences.qza \
  --i-reference-sequences References/gg2-c97-seqs-derep-lca.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table Cluster/table-cr-97.qza \
  --o-clustered-sequences Cluster/rep-seqs-cr-97.qza \
  --o-unmatched-sequences Cluster/unmatched-cr-97.qza
```

``` bash
qiime vsearch cluster-features-open-reference \
  --i-table denoise/table.qza \
  --i-sequences denoise/rep-sequences.qza \
  --i-reference-sequences References/gg2-c97-seqs-derep-lca.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table open-clust/table-cr-97.qza \
  --o-clustered-sequences open-clust/rep-seqs-cr-97.qza \
  --o-new-reference-sequences open-clust/new-ref-seqs-or-97.qza
```

### 4. Classifying and assigning taxonomy using the pre-trained model
``` bash 
qiime feature-classifier classify-sklearn \
  --i-classifier 2022.10.backbone.full-length.nb.sklearn-1.4.2.qza \
  --i-reads Cluster/rep-seqs-cr-97.qza \
  --o-classification taxonomy.qza
```

For the open reference

``` bash 
qiime feature-classifier classify-sklearn \
  --i-classifier 2022.10.backbone.full-length.nb.sklearn-1.4.2.qza \
  --i-reads open-clust/rep-seqs-cr-97.qza \
  --o-classification open-clust/open-taxonomy.qza
```

``` bash 
qiime metadata tabulate \
  --m-input-file open-clust/open-taxonomy.qza \
  --o-visualization open-clust/open-taxonomy.qzv
```


### 5. Export the qza frequency table 
using qiime tools export  command
``` bash 
# for the closed-reference-clustering 
qiime tools export --input-path Cluster/table-cr-97.qza --output-path Cluster/table-cr-97.biom

# For the Open-reference-clustering 
qiime tools export --input-path open-clust/table-cr-97.qza --output-path open-clust/table-cr-97.biom


qiime tools export --input-path open-clust/open-taxonomy.qza --output-path open-clust/taxonomy.biom
```


#### - Convert BIOME to a Text file
``` bash 

# Closed 
biom convert \
  --input-fp Cluster/table-cr-97.biom/feature-table.biom \
  --output-fp Cluster/table-cr-closed-otu-abundance.txt \
  --to-tsv
# Open 
biom convert \
  --input-fp open-clust/table-cr-97.biome/feature-table.biom \
  --output-fp open-clust/table-cr-open-otu-abundance.txt \
  --to-tsv

# Taxonomy 
biom convert \
  --input-fp taxonomy.biom/feature-table.biom \
  --output-fp taxonomy.txt \
  --to-tsv

```
