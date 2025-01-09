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
Multiple Genomes and datasets will be required to complete the pipeline:
The human genome can be found here: https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/

The Maize reference genome can be found here (Get both the genome and gff annotation file): https://www.maizegdb.org/assembly

Download all required diamond databases via Eggnog (requires installation) and currate custom databases with:
```
download_eggnog_data.py --data_dir /path/to/data/storage/Eggnog2

# create any custom, taxa specific databases with:
create_dbs.py -m diamond --dbname Microbes --taxa Fungi,Bacteria,Archaea --data_dir /path/to/data/storage/Eggnog2
```


## Example Header
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

# Trimming with Trimmomatic
for fastq_F in ./RawData/*R1_001.fastq.gz; do
    #Derive Reverse Read
    fastq_R=${fastq_F/R1_001.fastq/R2_001.fastq.gz}
    Name=$(basename "${fastq_F}" _R1_001.fastq.gz)

trimmomatic PE -phred33 -trimlog "${Name}_Trimlog.txt" "${fastq_F}" "${fastq_R}" \
    "./RawData/${Name}_paired_F.fastq.gz" "./RawData/${Name}_unpaired_F.fastq.gz" \
    "./RawData/${Name}_paired_R.fastq.gz" "./RawData/${Name}_unpaired_R.fastq.gz" \
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

for fastq_F in .RawData/*_paired_F.fastq.gz; do
    #Derive Reverse Read
    fastq_R=${fastq_F/_paired_F.fastq/_paired_R.fastq}
    Name=$(basename "${fastq_F}" _paired_F.fastq)

    # Run HISAT2
    hisat2 -p 32 -x "${human_index}" -1 "${fastq_F}" -2 "${fastq_R}" --no-mixed --no-discordant -S "./Output/${Name}_mapped2human.sam"

    # Sort the SAM file and convert to BAM
    samtools sort -o "./Output/${Name}_mapped2human.sorted.bam" "./Output/${Name}_mapped2human.sam" -@ 32

    # Extract unmapped reads
    samtools view -b -f 4 "./Output/${Name}_mapped2human.sorted.bam" > "./Output/${Name}_humanremoved.bam" -@ 32

    # Convert BAM to FASTQ
#    samtools bam2fq "./Output/${Name}_humanremoved.bam" > "./Output/${Name}_humanremoved.fastq.gz" -@ 32
    bedtools bamtofastq -i "./Output/${Name}_humanremoved.bam" -fq "./Output/${Name}_humanremoved_F.fastq.gz" -fq2 "./Output/${Name}_humanremoved_R.fastq.gz"
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

maize_index="./RawData/Index/B73Index"

for fastq_F in "${input_dir}"/*_riboremoved_F.fastq; do
    #Derive Reverse Read
    fastq_R=${fastq_F/F.fastq/R.fastq}
    Name=$(basename "${fastq_F}" _riboremoved_F.fastq)

    # Run HISAT2
    hisat2 -p 32 -x "${maize_index}" -1 "${fastq_F}" -2 "${fastq_R}" --no-mixed --no-discordant -S "./Output/${Name}_mapped2maize.sam"

    # Sort the SAM file and convert to BAM
    samtools sort -n -o "./Output/${Name}_mapped2maize.sorted.bam" "./Output/${Name}_mapped2maize.sam" -@ 32

    # Extract mapped reads
    samtools view -b -F 4 "./Output/${Name}_mapped2maize.sorted.bam" > "./Output/${Name}_mapped2maize.sorted1.bam" -@ 32

    # Convert BAM to FASTQ - use samtools to export a single paired output file and bedtools to create 2 paired-end files - which you'll need depends on your downstream purposes
    samtools bam2fq "./Output/${Name}_mapped2maize.sorted1.bam" > "./Output/${Name}_mapped2maize.fastq" -@ 32
    bedtools bamtofastq -i "./Output/${Name}_mapped2maize.sorted1.bam" -fq "./Output/${Name}_mapped2maize_F.fastq" -fq2 "./Output/${Name}_mapped2maize_R.fastq"

    # Extract unmapped reads
    samtools view -b -f 4 "./Output/${Name}_mapped2maize.sorted.bam" > "./Output/${Name}_maizeremoved.bam" -@ 32

    # Convert BAM to FASTQ
    samtools bam2fq "./Output/${Name}_maizeremoved.bam" > "./Output/${Name}_maizeremoved.fastq" --threads 32
    bedtools bamtofastq -i "./Output/${Name}_maizeremoved.bam" -fq "./Output/${Name}_maizeremoved_F.fastq" -fq2 "./Output/${Name}_maizeremoved_R.fastq"
done
```
## Step 4a Process maize reads
Calculate tpm
```
# If you have a .gff3 file convert it to gtf with:
agat_convert_sp_gff2gtf.pl --gff ../RawData/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.59.chr.gff3 -o ../RawData/Zea_mays_B73_v5.1.gtf

# Calculate tpm using feature counts package:
for bam in ./Output/*mapped2maize.sorted.bam; do
    Name=$(basename "${bam}" mapped2maize.sorted.bam)

    #convert to feature counts table
    featureCounts -p -t gene -T 6 -a ./RawData/Zea_mays_B73_v5.1.gtf -o "./Output/${Name}_MaizeFeature_counts.txt" ${bam}
done
```
Conduct DEseq
```
#In R:

```
## Step 5 Assign Bacterial / Fungal taxonomy to remaining reads using Kraken2
```
# Download all relvant databases and build the Kraken Library:
kraken2-build --download-library bacteria --db /path/to/databases/Kraken2
kraken2-build --build --db ./Kraken2 --threads 28

# I've had problems with Kraken not using the number of threads I specify but fixed it by adding this line to my script:
export OMP_NUM_THREADS=28

# A loop to run Kraken:
for fastq_F in ./Output/*_maizeremoved_F.fastq; do
    # Derive Reverse Read
    fastq_R=${fastq_F/F.fastq/R.fastq}
    Name=$(basename "${fastq_F}" _maizeremoved_F.fastq)

    # Run Kraken2 on the bacterial reference
    kraken2 --db /share/gage/njkorth/GERMsD_TEST/RawData/Index/Kraken2 --threads 28 --paired "${fastq_F}" "${fastq_R}" \
    --output "./KrakenOut/${Name}_Kraken2_Out.txt" --report "./KrakenOut/${Name}_Kraken2_Report.txt" \
    --classified-out "./KrakenOut/${Name}_mapped2bacteria#.fastq"

    # Run Kraken2 on the fungal reference
    kraken2 --db /share/gage/njkorth/GERMsD_TEST/RawData/Index/K2Fungi --threads 12 --paired "${fastq_F}" "${fastq_R}" \
    --output "./KrakenOut/${Name}_Kraken2F_Out.txt" --report "./KrakenOut/${Name}_Kraken2F_Report.txt" --classified-out "./KrakenOut/${Name}_mapped2fungi#.fastq"

done
```

## Step 6, Assign functionality to reads using Eggnog:


## Contact
For clarification on code missing annotation contact:
* Nate Korth: njkorth@ncsu.edu / nate.korth@gmail.com
* Joe Gage: jlgage@ncsu.edu

Whatever parts of the pipeline you use, please remember to cite all relevant packages
