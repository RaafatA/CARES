
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest-2.txt \
  --output-path sequences-2.qza \
  --input-format PairedEndFastqManifestPhred33V2


qiime dada2 denoise-paired   --i-demultiplexed-seqs sequences-2.qza   --p-trunc-len-f 250   --p-trunc-len-r 250   --p-n-threads 0  --o-table Results/denoise-dada2/table-2.qza   --o-representative-sequences Results/denoise-dada2/rep-seqs-2.qza   --o-denoising-stats Results/denoise-dada2/denoising-stats-2.qza




qiime feature-table summarize \
  --i-table Results/denoise-dada2/table-2.qza \
  --o-visualization Results/denoise-dada2/table-2.qzv \


qiime feature-table tabulate-seqs \
  --i-data Results/denoise-dada2/rep-seqs-2.qza \
  --o-visualization Results/denoise-dada2/rep-seqs-2.qzv
qiime metadata tabulate \
  --m-input-file Results/denoise-dada2/denoising-stats-2.qza \
  --o-visualization Results/denoise-dada2/denoising-stats-2.qzv


export TMPDIR=/media/raafat/927d53ff-9823-4e2d-9c98-cee61adf83f7/home/rafat
qiime feature-classifier classify-sklearn \
  --i-classifier GG2024/2024.09.backbone.full-length.nb.qza \
  --i-reads Results/denoise-dada2/rep-seqs-2.qza \
  --o-classification Results/denoise-dada2/taxonomy--2.qza

qiime tools export --input-path Results/denoise-dada2/rep-seqs-2.qza --output-path Results/denoise-dada2/sequences-ASVs-2.biom
# Exported Results/denoise-dada2/rep-seqs-2.qza as DNASequencesDirectoryFormat to directory Results/denoise-dada2/sequences-ASVs-2.biom
