# Create a database for blast
makeblastdb -in egyptian_clover_genome.fasta -dbtype nucl -out trifolium_db

blastn \
  -query sod-Arabidopsis_thaliana.fasta \
  -db trifolium_db \
  -out trifolium_SOD_hits_with_seq.tsv \
  -outfmt "7 qseqid sseqid pident length mismatch gapopen evalue bitscore qseq sseq" \
  -evalue 1e-5 \
  -word_size 7 \
  -reward 2 \
  -penalty -3 \
  -gapopen 5 \
  -gapextend 2 \
  -num_threads 8
