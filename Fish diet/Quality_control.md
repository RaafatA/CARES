# Step 1: Download and prepare the reference genome
wget https://ftp.ensembl.org/pub/release-113/fasta/oreochromis_niloticus/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa.gz

# Step 2: Index the reference genome using BWA
bwa index Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.toplevel.fa.gz

# Step 3: Define the list of samples
samples=("NT5" "NT11" "NT14" "NT16" "NT17" "NT20")

# Step 4: Loop through each sample and align the reads
for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"
    
    # Define forward and reverse read files
    forward_reads="${sample}.raw_1.fastq.gz"
    reverse_reads="${sample}.raw_2.fastq.gz"
    
    # Align paired-end reads to the reference genome using BWA
    bwa mem -t 4 nile_tilapia_genome.fasta "$forward_reads" "$reverse_reads" > "${sample}_aligned_reads.sam"
    
    # Convert SAM to BAM and sort
    samtools view -bS "${sample}_aligned_reads.sam" > "${sample}_aligned_reads.bam"
    samtools sort "${sample}_aligned_reads.bam" -o "${sample}_aligned_reads.sorted.bam"
    
    # Remove host DNA contamination
    samtools view -b -f 4 "${sample}_aligned_reads.sorted.bam" > "${sample}_unaligned_reads.bam"
    
    # Extract unaligned reads
    samtools fastq -f 4 "${sample}_unaligned_reads.bam" > "${sample}_unaligned_reads.fastq"
done
