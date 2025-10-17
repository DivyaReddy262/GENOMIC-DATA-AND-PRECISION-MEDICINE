Name: Divya Reddy Konda
Course: Genomic data analysis & precision medicine
Assignment: Genome Assembly (Velvet vs Oases)
Date: 10/17/2025
Environment: IU HPC (BigRed200)

Dataset:
SRA Accession: SRR21904868
Data type: Illumina short-read DNA sequencing
Reference: E. coli genome (~4.6 Mb)

Environment Setup:
# Load Anaconda module on IU HPC
module load anaconda

# Create and activate environment
conda create -n ecasm -y fastqc trimmomatic sra-tools velvet oases python
conda activate ecasm

# Create working directories
mkdir -p ~/Genomics/{data/raw,data/trimmed,results/velvet}

STEP 1: Download SRA Data
# Use SRA-Toolkit
cd ~/Genomics/data/raw
prefetch SRR21904868
fasterq-dump SRR21904868

# Resulting files:
  SRR21904868_1.fastq
  SRR21904868_2.fastq

STEP 2: Quality Control
# Check quality of raw reads
fastqc SRR21904868_1.fastq SRR21904868_2.fastq -o .

STEP 3: Trimming
# Remove adapters and low-quality bases
cd ~/Genomics/data
mkdir trimmed
trimmomatic PE \
  raw/SRR21904868_1.fastq raw/SRR21904868_2.fastq \
  trimmed/1_paired.fq trimmed/1_unpaired.fq \
  trimmed/2_paired.fq trimmed/2_unpaired.fq \
  SLIDINGWINDOW:4:20 MINLEN:50


STEP 4: Genome Assembly - Velvet
cd ~/Genomics/results/velvet
# Run Velvet for multiple k-mer sizes
for k in 31 41 51 61 71 81 91; do
  velveth vel_k$k $k -shortPaired -fastq ../../data/trimmed/1_paired.fq ../../data/trimmed/2_paired.fq
  velvetg vel_k$k
done

STEP 5: Transcriptome Assembly - Oases
# Oases requires Velvet output directories
for k in 31 41 51 61 71 81 91; do
  oases vel_k$k
done

STEP 6: Evaluate Assemblies
# Summarize N50, total bases, number of contigs
python3 compare_assemblies.py   # (script provided)

# Best assembly results:
# Velvet: k = 91, N50 = 112,588 bp, contigs = 230, total length = 4,785,447 bp
# Oases:  k = 61, N50 = 751,529 bp, transcripts = 1,434, total bases = 13,975,317 bp

STEP 7: Contig Length Reporting
# Extract length of each contig / transcript
# Velvet (best K=91)
awk '/^>/{if(L){print L}; L=0; next} {L+=length($0)} END{if(L)print L}' \
  vel_k91/contigs.fa > vel_k91/contig_lengths.txt

# Oases (best K=61)
awk '/^>/{if(L){print L}; L=0; next} {L+=length($0)} END{if(L)print L}' \
  vel_k61/transcripts.fa > vel_k61/oases_transcript_lengths.txt

STEP 9: Results Summary
Velvet:
- Best K-mer: 91
- N50: 112,588 bp
- Contigs: 230
- Total assembly length: 4,785,447 bp
- Largest contig: 366,094 bp

Oases:
- Best K-mer: 61
- N50: 751,529 bp
- Transcripts: 1,434
- Total bases: 13,975,317 bp
- Largest transcript: 1,204,768 bp

STEP 10: Final Files
- velvet_table.tsv
- oases_table.tsv
- vel_k91/contig_lengths.txt
- vel_k61/oases_transcript_lengths.txt

Comparision:

• Genome assembly was performed using both Velvet and Oases at multiple k-mer values 
(31, 41, 51, 61, 71, 81, and 91). 
• For Velvet, the optimal assembly was achieved at k = 91, with an N50 of 112,588 bp, 230 
contigs, a total assembly length of 4.79 Mb, and the largest contig measuring 366,094 bp. 
For Oases, the best performance was obtained at k = 61, yielding a much higher N50-like 
of 751,529 bp, but with 1,434 transcripts and a total base count of 13.98 Mb, with the 
largest transcript at 1,204,768 bp. 
• Although Oases produced longer N50 values, it is designed primarily for transcriptome 
assembly and generated many more fragments, inflating total assembly length beyond the 
expected E. coli genome size. In contrast, Velvet generated a more accurate, coherent 
genome assembly closer to the known size (4.6 Mb). 
• Therefore, Velvet is the more suitable assembler for this genome dataset.

Notes:
* Velvet is better suited for genome assembly.
* Oases is transcriptome-focused and inflated genome size in this dataset.
* Final assembly selected: Velvet at k=91.