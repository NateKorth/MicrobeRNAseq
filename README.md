# Maize Microbe RNAseq Analysis Pipeline
This pipeline is a series of unix commands built to be applied in a IBM bsub submission system as employed by the North Carolina State University High Performance Computing center. 
The Conda environment(s) can be found in the environments folder and can be built independently with the most current versions of the software by installing the Required Packages to your own Conda environment. Some scripts for visualization/stats are written in R and can be easily applied 

## Repo information
This pipeline is designed to analyze RNAseq data. It seperates reads from plants, bacteria, and fungi and attempts to assign functionality / pathway information to each read
To use it, first clone the directory by entering the following command in your terminal: 
```
git clone git@github.com:NateKorth/MicrobeRNAseq.git
```

## Required Packages
The following can be installed to a conda environment with the command conda install
* ncbi-genome-download
* fastqc
* hisat2
* pytorch
* ribodetector
* DEseq2
* Kraken2
* eggNog-mapper

## Step0 Download Required Genomes
The human genome can be found here: https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/

The Maize reference genome can be found here (Get both the genome and gff annotation file): https://www.maizegdb.org/assembly

Many of the following steps should be batched to a high performance cpu or they will take a very long time to run, here is an example header for the NCSU HPC cluster (IBM bsub): - just make sure the number of threads and amount of memory you request here are in line to what each command you run requires:
```
#!/bin/bash
#BSUB -n 16
#BSUB -R "rusage[mem=65GB]"
#BSUB -W 75:00
#BSUB -q gage
#BSUB -J ProjectName
#BSUB -o out.%J
#BSUB -e err.%J
#BSUB -R "span[hosts=1]"
```

## Step 1 Trim/QC with fastqc and trimmomatic
```
conda activate myenvironment
#Make directory for fastqc ouput:
mkdir ./RawData/fastqc
#Run fastQC
for fastq_F in ./RawData/*R1_001.fastq.gz; do
    #Derive Reverse Read
    fastq_R=${fastq_F/R1_001.fastq/R2_001.fastq}
    Name=$(basename "${fastq_F}" _R1_001.fastq)
    fastqc -t 2 "${fastq_F}" "${fastq_R}" --outdir ./RawData/fastqc
done
```
After visually inspecting fastqc output decide on trimming, minimal trimming may be required
```
# Trimming with Trimmomatic
for fastq_F in ./RawData/*R1_001.fastq; do
    #Derive Reverse Read
    fastq_R=${fastq_F/R1_001.fastq/R2_001.fastq}
    Name=$(basename "${fastq_F}" _R1_001.fastq)

    #See the trimmomatic manual for specifics in the following command:

trimmomatic PE -phred33 -trimlog "${Name}_Trimlog.txt" "${fastq_F}" "${fastq_R}" \
    "./RawData/${Name}_paired_F.fastq" "./RawData/${Name}_unpaired_F.fastq" \
    "./RawData/${Name}_paired_R.fastq" "./RawData/${Name}_unpaired_R.fastq" \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:50

done
```

## Step 2 Align reads to human genome and remove from downstream analysis
```
# Build human Index:
gzip -cd GRCh38_latest_genomic.fna.gz > HumanRef.fa
hisat2-build -p 16 HumanRef.fa HumanRefIndex
# now that the index is built, remove the unzipped fasta file to save space:
rm HumanRef.fa

# Loop over all fastq files in the input directory
human_index="./RawData/Index/HumanRefIndex"

for fastq_F in "${input_dir}"/*R1_001.fastq; do
    #Derive Reverse Read
    fastq_R=${fastq_F/R1_001.fastq/R2_001.fastq}
    Output=$(basename "${fastq_F}" R1_001.fastq)

    # Run HISAT2
    hisat2 -p 32 -x "${human_index}" -1 "${fastq_F}" -2 "${fastq_R}" --no-mixed --no-discordant -S "./Output/${Output}_mapped2human.sam"

    # Sort the SAM file and convert to BAM
    samtools sort -o "${output_dir}/${Output}_mapped2human.sorted.bam" "${output_dir}/${Output}_mapped2human.sam" -@ 32

    # Extract unmapped reads
    samtools view -b -f 4 "${output_dir}/${Output}_mapped2human.sorted.bam" > "./Output/${Output}_humanremoved.bam" -@ 32

    # Convert BAM to FASTQ
    samtools bam2fq "./Output/${Output}_humanremoved.bam" > "./Output/${Output}_humanremoved.fastq" -@ 32
    bedtools bamtofastq -i "./Output/${Output}_humanremoved.bam" -fq "./Output/${Output}_humanremoved_F.fastq" -fq2 "./Output/${Output}_humanremoved_R.fastq"

done

# A few file deletions/conversions to save space:
rm ./Output/*.sam
```

## Step 3 Remove rRNA reads
If possible do this step on a gpu, at NCstate the header looks like this:
```
#!/bin/bash
#BSUB -n 2
#BSUB -R "rusage[mem=125GB]"
#BSUB -W 990
#BSUB -R "span[hosts=1]"
#BSUB -q gpu
#BSUB -gpu "num=1:mode=shared:mps=no"
#BSUB -J ribodetect2
#BSUB -o out.%J
#BSUB -e err.%J
```
Ribodetetctor is a deep learning method to identify rRNA, read more about it and how to install it as a conda environment here: https://github.com/hzi-bifo/RiboDetector
```
conda activate /usr/local/usrapps/gage/njkorth/ribodetect

#Loops to run ribodetector on all fastq files in Output Folder:
for fastq_F in ./Output/*humanremoved_F.fastq; do
    #Derive Reverse Read
    fastq_R=${fastq_F/F.fastq/R.fastq}
    Name=$(basename "${fastq_F}" _humanremoved_F.fastq)

    # Run Ribodetect
    ribodetector -l 151 -t 2  -i "${fastq_F}" "${fastq_R}" -m 20  -e rrna -o "./Output/${Name}_riboremoved_F.fastq" "./Output/${Name}_riboremoved_R.fastq" --log "${Name}_log.txt"
done

```
## Step 4 Align reads to maize genome
```
# Build Maize Index:
gzip -cd ./Zm-B73-REFERENCE-NAM-5.0.fa.gz > B73Refv5.fa
hisat2-build -p 16 B73Refv5.fa B73Index
# now that the index is built, remove the unzipped fasta file to save space:
rm BB73Refv5.fa

# Loop over all fastq files in the input directory
for fastq_F in "${input_dir}"/*humanremoved_F.fastq; do
    #Derive Reverse Read
    fastq_R=${fastq_F/humanremoved_F.fastq/humanremoved_R.fastq}
    Name=$(basename "${fastq_F}" _humanremoved_F.fastq)

    # Run HISAT2
    hisat2 -p 32 -x "${maize_index}" -1 "${fastq_F}" -2 "${fastq_R}" --no-mixed --no-discordant -S "./Output/${Name}_mapped2maize.sam"

    # Sort the SAM file and convert to BAM
    samtools sort -n -o "./Output/${Name}_mapped2maize.sorted.bam" "./Output/${Name}_mapped2maize.sam" -@ 32

    # Extract mapped reads
    samtools view -b -F 4 "./Output/${Name}_mapped2maize.sorted.bam" > "./Output/${Name}_mapped2maize.sorted1.bam" -@ 32

    # Convert BAM to FASTQ
    samtools bam2fq "./Output/${Name}_mapped2maize.sorted1.bam" > "./Output/${Name}_mapped2maize.fastq" -@ 32
    #samtools bam2fq -n "./Output/${Name}_mapped2maize.sorted.bam" -1 "./Output/${Name}_mapped2maize_F.fastq" -2 "./Output/${Name}_mapped2maize_R.fastq"
    bedtools bamtofastq -i "./Output/${Name}_mapped2maize.sorted1.bam" -fq "./Output/${Name}_mapped2maize_F.fastq" -fq2 "./Output/${Name}_mapped2maize_R.fastq"

    # Extract unmapped reads
    samtools view -b -f 4 "./Output/${Name}_mapped2maize.sorted.bam" > "./Output/${Name}_maizeremoved.bam" -@ 32

    # Convert BAM to FASTQ
    #samtools bam2fq "./Output/${Name}_maizeremoved.bam" > "./Output/${Name}_maizeremoved.fastq" --threads 32
    bedtools bamtofastq -i "./Output/${Name}_maizeremoved.bam" -fq "./Output/${Name}_maizeremoved_F.fastq" -fq2 "./Output/${Name}_maizeremoved_R.fastq"
done
```
## Step 4a Process maize reads
Calculate tpm
```

```
Conduct DEseq

## Step 5 Assign Bacterial / Fungal taxonomy to remaining reads

## Contact
For clarification on code missing annotation contact:
* Nate Korth: njkorth@ncsu.edu / nate.korth@gmail.com
* Joe Gage: jlgage@ncsu.edu
