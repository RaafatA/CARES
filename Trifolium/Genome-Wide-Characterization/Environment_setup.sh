# Clover Virtual Environment Creation
conda create -n Clover python=3.10
conda activate Clover

# Sequence search & domain analysis
conda install -c bioconda blast hmmer interproscan diamond
conda install -c bioconda seqkit emboss

# Extract and manipulate genomic regions
conda install -c bioconda bedtools samtools seqtk

