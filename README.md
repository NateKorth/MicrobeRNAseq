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

Sorghum reference genome:  https://ftp.sorghumbase.org/release-9/fasta/sorghum_bicolorv5/dna/Sorghum_bicolorv5.Sb-BTX623-REFERENCE-JGI-5.1.dna.toplevel.fa.gz
Sorghum gff annotation file:  https://ftp.sorghumbase.org/release-9/gff3/sorghum_bicolorv5/Sorghum_bicolorv5.Sb-BTX623-REFERENCE-JGI-5.1.gff3.gz

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
There are some pre-built Kraken2 indexes here:https://benlangmead.github.io/aws-indexes/k2

```
#example header for an array job:
#!/bin/bash
#BSUB -n 12
#BSUB -R "rusage[mem=60GB]"
#BSUB -W 150:00
#BSUB -J Kraken2[1-16]
#BSUB -o logs/K2.out.%J.%I
#BSUB -e logs/K2.err.%J.%I
#BSUB -R "span[hosts=1]"

#Activate a kraken2 environment or have Kraken2 installed on your workspace
conda activate Kraken2

# I've had problems with Kraken not using the number of threads I specify but fixed it by adding this line to the top of my scripts:
export OMP_NUM_THREADS=12

# Download all relvant databases and build the Kraken Library:

#gtdb full database (this is a very large set of bacterial genomes
aria2c -c -x 16 -s 16 -k 1M -o gtdb_genomes_reps.tar.gz "https://genome-idx.s3.amazonaws.com/kraken/k2_gtdb_genome_reps_20241109.tar.gz"

#Standard Database including fungi:
aria2c -c -x 16 -s 16 -k 1M -o k2_StdplusF.tar.gz "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20241228.tar.gz"

#Decompile it:
tar --use-compress-program=pigz -xvf k2_StdplusF.tar.gz -C .

# OR (if you have problems with the pre-built database, use k2 build:
k2 build --db /path/to/databases/Kraken2 --standard --threads 12

# Batch several Kraken2 jobs as an array:
# File processing
FASTQ_FILES=($(ls ./Output/*maizeremoved_F.fastq))
SAMPLE=${FASTQ_FILES[$((LSB_JOBINDEX-1))]}

# Derive Reverse Read
fastq_F=$SAMPLE
fastq_R=${fastq_F/F.fastq/R.fastq}
Name=$(basename "${fastq_F}" _maizeremoved_F.fastq)

#Loop through all fastqs with human seqs removed and run Kraken2:

#Use the --memory-mapping flag if the active memory you have avaliable is less than the size of the database you're using (Kraken will otherwise load the whole database into active memory)

k2 classify --db /path/to/databases/Kraken2 --threads 12 --paired "${fastq_F}" "${fastq_R}" \
    --output "./KrakenOut/${Name}_Kraken2_Out.txt" --report "./KrakenOut/${Name}_Kraken2_Report.txt" \
    --classified-out "./KrakenOut/${Name}_mapped2microbes#.fastq" --memory-mapping --quick

#rezip all the files:
gzip ./KrakenOut/${Name}_mapped2microbes_1.fastq
gzip ./KrakenOut/${Name}_mapped2microbes_2.fastq

conda deactivate

```

## Step 6, Assign functionality to reads using Eggnog:
This script is in the form of a loop but if you have a lot of samples, might change batch each job seperatly using an array as done in the previous section
```
# First load anaconda manager with eggnog installed and download eggnog database
# Make a database for just Microbes (Customize as you need):
create_dbs.py -m diamond --dbname Microbes --taxa Bacteria,Fungi,Archaea --data_dir /rs1/researchers/j/jlgage/users/njkorth/databases/Eggnog2

# Unzip any gzipped files before running:
gunzip ./Output/*maizeremoved*fastq.gz

#Combine fastqs into a single fasta
for fastq_F in ./Output/*maizeremoved_F.fastq; do
    #Derive Reverse Read
    fastq_R=${fastq_F/F.fastq/R.fastq}
    Name=$(basename "${fastq_F}" _maizeremoved_F.fastq)
    #Pair fastqs with Pear
    pear -f "$fastq_F" -r "$fastq_R" -o "./Output/${Name}" -j 30 -y 50G -p 0.05 -v 4 -m 450 -g 2
    #Combine the assemble and unassembled reads
    cat "./Output/${Name}.assembled.fastq" "./Output/${Name}.unassembled.forward.fastq" "./Output/${Name}.unassembled.reverse.fastq" > "./Output/${Name}.fastq"

#Annotate function with diamond:
#Using a lower e value cutoff to control for any plant reads being miss-annotated as fungi:

    emapper.py -m diamond --no_annot --no_file_comments --cpu 32 --data_dir /path/to/databases/Eggnog2 --seed_ortholog_evalue 1e-3 \
    -i "./Output/${Name}.fastq" -o "${Name}_e3" --itype CDS --output_dir ./Output/Annotation --dmnd_db /path/to/databases/Eggnog2/Microbes.dmnd

#Add gene names

    emapper.py --annotate_hits_table "./Output/Annotation/${Name}_e3.emapper.seed_orthologs" --data_dir path/to/databases/Eggnog2 --seed_ortholog_evalue 1e-3 \
    --no_file_comments -o "${Name}_e3" --dbmem --output_dir ./Output/Annotation --cpu 32 --dmnd_db path/to/databases/Eggnog2/Microbes.dmnd

done

#cleanup:
rm ./Output/*discarded*fastq
gzip ./Output/*fastq
```
## Step 6.1 Kraken taxonomy assignment
```
#In Bash:
#Install Kraken2 - v2.14 (I use conda)
conda activate Kraken2 
# Download kraken database: This is a already built kraken database, can also usce kraken2-build --download library to build your database
## Can find some pre built indexes here:https://benlangmead.github.io/aws-indexes/k2****
aria2c -c -x 16 -s 16 -k 1M -o gtdb_genomes_reps.tar.gz "https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20241228.tar.gz"
# Deompile it
tar --use-compress-program=pigz -xvf gtdb_genomes_reps.tar.gz -C .

# File processing (This script uses an array type submission on an ibm submission system at the NC State hpc)
FASTQ_FILES=($(ls ./Output/*maizeremoved_F.fastq))
SAMPLE=${FASTQ_FILES[$((LSB_JOBINDEX-1))]}

# Derive Reverse Read
fastq_F=$SAMPLE
fastq_R=${fastq_F/F.fastq/R.fastq}
Name=$(basename "${fastq_F}" _maizeremoved_F.fastq)

#Loop through all fastqs with human seqs removed and run Kraken2:

    # Run Kraken2 on the bacterial reference
    kraken2 --db path/to/your/databases/Kraken2 --threads 10 --paired "${fastq_F}" "${fastq_R}" \
    --output "./KrakenOut/${Name}_Kraken2_Out.txt" --report "./KrakenOut/${Name}_Kraken2_Report.txt" \
    --classified-out "./KrakenOut/${Name}_mapped2bacteria.fastq"

    # Run Kraken2 on the fungal reference
    kraken2 --db path/to/your/databases/Kraken2 --threads 12 --paired "${fastq_F}" "${fastq_R}" \
    --output "./KrakenOut/${Name}_Kraken2F_Out.txt" --report "./KrakenOut/${Name}_Kraken2F_Report.txt" --classified-out "./KrakenOut/${Name}_mapped2fungi#.fastq"

done

#rezip all the files:
gzip ./KrakenOut/${Name}_mapped2bacteria.fastq

#Read how many sequences are in each fastq:
for file in ./KrakenOut/*_F.fastq; do echo "$file: $(( $(grep -c '^@' "$file") / 4 )) reads";done

conda deactivate
```

## Step 6.2, assign pathway information to microbial reads and collapse by pathway and gene level
```
#In R:
# EggNOG Annotation Processing Script
# Workflow:
#   1. Read EggNOG annotation tables
#   2. Select necessary columns
#   3. Filter rows for Bacteria or Fungi
#   4. Collapse by taxa + name + desc, count occurrences
#   5. Save checkpoint CSV per sample
#   6. Combine all checkpoints into one big table

# --- Libraries ---
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(readxl)
  library(tidyverse)
  library(KEGGREST)
})

# --- Directories ---
indir  <- "D:/Data/Path/to/EggnogAnnotations"
outdir <- "./"

# --- Get input files ---
files <- list.files(indir, pattern="*annotations", full.names=TRUE)

# --- Get already processed checkpoint files ---
complete <- list.files(outdir, pattern="*_TaxaCounts.csv", full.names=TRUE)
complete <- gsub("_TaxaCounts.csv$", "", basename(complete))

# --- Filter out completed files (don’t reprocess) ---
files <- setdiff(basename(files), complete)

# --- Define processing function ---
process_annotations <- function(txt_file, indir, outdir) {
  infile <- file.path(indir, txt_file)
  message("Processing: ", infile)
  
  # Read only needed columns
  dt <- fread(
    infile,
    select = c("eggNOG_OGs", "max_annot_lvl", "Preferred_name", "Description", "seed_ortholog",
               "evalue","KEGG_ko","KEGG_Module","CAZy","EC","GOs","PFAMs")
  )
  
  # Rename for consistency
  setnames(dt, 
           c("eggNOG_OGs", "max_annot_lvl", "Preferred_name", "Description", "seed_ortholog",
             "evalue","KEGG_ko","KEGG_Module","CAZy","EC","GOs","PFAMs"),
           c("eggNOG_OGs", "taxa", "name", "desc","gene", "evalue","pathway_Kegg","pathway_KeggMod","pathway_CAZy","pathway_EC",
             "pathway_GO","pathway_PFAM"))
  
  # Keep only Bacteria or Fungi
  dt <- dt[grepl("Bacteria|Fungi", eggNOG_OGs, ignore.case = TRUE)]
  dt <- dt[evalue <= 1e-6]
  
  # Collapse
  gene_counts <- dt[, .(count = .N), by = gene]
  meta_lookup <- unique(dt[, .(eggNOG_OGs,gene, taxa, name, desc,pathway_Kegg,pathway_KeggMod,pathway_CAZy,pathway_EC,pathway_GO,pathway_PFAM)])
  dt_summary <- merge(gene_counts, meta_lookup, by = "gene", all.x = TRUE)
  setcolorder(dt_summary, c("eggNOG_OGs", "taxa", "name", "desc","gene","pathway_Kegg","pathway_KeggMod","pathway_CAZy","pathway_EC",
                            "pathway_GO","pathway_PFAM"))
  
  # Write checkpoint
  outfile <- file.path(outdir, paste0(txt_file, "_TaxaCounts.csv"))
  fwrite(dt_summary, outfile)
  
  return(outfile)
}

# --- Process all files (checkpoint stage) ---
checkpoint_files <- c()
for (f in files) {
  outfile <- process_annotations(f, indir, outdir)
  checkpoint_files <- c(checkpoint_files, outfile)
}

# --- Re-read all checkpoint files & combine ---
all_checkpoints <- list.files(outdir, pattern="*_TaxaCounts.csv", full.names=TRUE)

combined <- rbindlist(lapply(all_checkpoints, function(f) {
  dt <- fread(f)
  sample_id <- gsub("_TaxaCounts.csv$", "", basename(f))
  dt[, sample := sample_id]
  return(dt)
}), use.names = TRUE, fill = TRUE)

wide_df <- combined %>%
  dplyr::select(eggNOG_OGs,taxa, name, desc, gene,pathway_CAZy,pathway_EC,pathway_GO,pathway_Kegg,pathway_PFAM,pathway_KeggMod, sample, count) %>%
  pivot_wider(
    names_from = sample,
    values_from = count,
    values_fill = list(count = 0)
  )

sample_cols <- setdiff(colnames(wide_df), c("eggNOG_OGs","taxa","taxa2", "name", "desc", "gene","pathway_Kegg","pathway_KeggMod","pathway_CAZy","pathway_EC",
                                            "pathway_GO","pathway_PFAM","A_Test.emapper.annotations","B_Test.emapper.annotations"))

#Make a new column that is a sum of all reads per gene
wide_df <- wide_df %>%
  mutate(AllSamples = rowSums(across(all_of(sample_cols)), na.rm = TRUE))

wide_df <- wide_df %>%
  mutate(taxa2 = sub(".*\\|(.*)$", "\\1", eggNOG_OGs))

gc()

# Collapse the genes into one per taxa:
collapsed_df <- wide_df %>%
  # remove rows with "-" in name and desc
  filter(name != "-" & desc != "-") %>%
  # standardize to lowercase to ignore capitalization
  mutate(
    taxa2 = str_to_lower(taxa2),
    name = str_to_lower(name),
    pathway_EC=str_to_lower(pathway_EC)
  ) %>%
  # group by the unique combination of taxa, name, desc
  group_by(taxa2, name,pathway_EC) %>%
  # sum all sample columns
  summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop")

write.csv(collapsed_df,"Microbe_Genes_counts.csv")

# Collect pathway data for later use
gene2pathway <- collapsed_df %>%
  dplyr::select(name, pathway_EC, taxa2) %>%
  filter(!is.na(name) & name != "") %>%
  distinct()

# Split ECs into separate rows and sum counts per EC
df_ec <- gene2pathway %>%
  filter(pathway_EC != "-") %>%
  separate_rows(pathway_EC, sep = ",") %>%
  mutate(pathway_EC = trimws(pathway_EC)) %>%
  group_by(taxa2, name,pathway_EC) %>%
  summarise(across(starts_with("H"), sum), .groups = "drop")

# EC → Pathway mapping (character vector)
ec2path <- keggLink("pathway", "ec")

ec2path_tbl <- tibble(
  pathway_EC = sub("ec:", "", names(ec2path)),
  pathway    = basename(ec2path))%>%
  filter(!grepl("ec", pathway))

# Get unique pathway descriptions
path_ids <- unique(ec2path_tbl$pathway)

# safer to query in chunks
get_desc <- function(pid) {
  tryCatch({
    keggList(pid) %>% unname()
  }, error = function(e) NA_character_)
}

path_descs <- map_chr(path_ids, get_desc)

path_tbl <- tibble(
  pathway = path_ids,
  desc_p    = path_descs
)

# Join descriptions back
ec2path_tbl <- ec2path_tbl %>%
  left_join(path_tbl, by = "pathway")

# Join counts and collapse
df_final <- df_ec %>%
  left_join(ec2path_tbl, by = "pathway_EC") %>%
  filter(!is.na(desc_p)) %>%   # only keep map pathways
  group_by(taxa2, name,pathway_EC,pathway, desc_p) %>%
  summarise(across(starts_with("H"), sum), .groups = "drop")

# Save to CSV
write.csv(df_final, "MicrobeGene2PathwayA.csv", row.names = FALSE)


# collapse the genes into one per taxa:
collapsed_df_ByGene <- collapsed_df %>%
  # remove rows with "-" in name and desc
  filter(name != "-") %>%
  # standardize to lowercase to ignore capitalization
  mutate(
    name = str_to_lower(name),
    desc = str_to_lower(desc),
    pathway_EC=str_to_lower(pathway_EC),
    pathway_Kegg=str_to_lower(pathway_Kegg)
  ) %>%
  # group by the unique combination of taxa, name, desc
  group_by(name,desc,pathway_EC,pathway_Kegg) %>%
  # sum all sample columns
  summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop")

gene2pathway <- collapsed_df_ByGene %>%
  dplyr::select(name, desc, pathway_EC, pathway_Kegg) %>%
  filter(!is.na(name) & name != "") %>%
  distinct()

# Split ECs into separate rows and sum counts per EC
df_ec <- gene2pathway %>%
  filter(pathway_EC != "-") %>%
  separate_rows(pathway_EC, sep = ",") %>%
  mutate(pathway_EC = trimws(pathway_EC)) %>%
  group_by(name, desc, pathway_EC, pathway_Kegg) %>%
  summarise(across(starts_with("H"), sum), .groups = "drop")

# EC → Pathway mapping (character vector)
ec2path <- keggLink("pathway", "ec")

ec2path_tbl <- tibble(
  pathway_EC = sub("ec:", "", names(ec2path)),
  pathway    = basename(ec2path))%>%
  filter(!grepl("ec", pathway))

# Get unique pathway descriptions
path_ids <- unique(ec2path_tbl$pathway)

# safer to query in chunks
get_desc <- function(pid) {
  tryCatch({
    keggList(pid) %>% unname()
  }, error = function(e) NA_character_)
}

#path_descs <- map_chr(path_ids, get_desc)

path_tbl <- tibble(
  pathway = path_ids,
  desc_p    = path_descs
)

# Step : join descriptions back
ec2path_tbl <- ec2path_tbl %>%
  left_join(path_tbl, by = "pathway")

# Join counts and collapse
df_final <- df_ec %>%
  left_join(ec2path_tbl, by = "pathway_EC") %>%
  filter(!is.na(desc_p)) %>%   # only keep map pathways
  group_by(name, desc, pathway_EC, pathway_Kegg,pathway, desc_p) %>%
  summarise(across(starts_with("H"), sum), .groups = "drop")

# Save to CSV
write.csv(df_final, "MicrobeGene2PathwayB.csv", row.names = FALSE)

# Pick your favorite annotation method and
# collapse the genes into pathways:
collapsed_df_ByPathway_GO <- wide_df %>%
  # remove rows with "-" in name and desc
  filter(name != "-") %>%
  # group by the unique combination of taxa, name, desc
  group_by(pathway_GO) %>%
  # sum all sample columns
  summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop")

collapsed_df_ByPathway_PFAM <- wide_df %>%
  # remove rows with "-" in name and desc
  filter(name != "-") %>%
  # group by the unique combination of taxa, name, desc
  group_by(pathway_PFAM) %>%
  # sum all sample columns
  summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop")

collapsed_df_ByPathway_EC <- wide_df %>%
  # remove rows with "-" in name and desc
  filter(name != "-") %>%
  # standardize to lowercase to ignore capitalization
  # group by the unique combination of taxa, name, desc
  group_by(pathway_EC) %>%
  # sum all sample columns
  summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop")

collapsed_df_ByPathway_KeggMod <- wide_df %>%
  # remove rows with "-" in name and desc
  filter(name != "-") %>%
  # standardize to lowercase to ignore capitalization
  # group by the unique combination of taxa, name, desc
  group_by(pathway_KeggMod) %>%
  # sum all sample columns
  summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop")

collapsed_df_ByPathway_CAZy <- wide_df %>%
  # remove rows with "-" in name and desc
  filter(name != "-") %>%
  # standardize to lowercase to ignore capitalization
  # group by the unique combination of taxa, name, desc
  group_by(pathway_CAZy) %>%
  # sum all sample columns
  summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop")

collapsed_df_ByPathway_Kegg <- wide_df %>%
  # remove rows with "-" in name and desc
  filter(name != "-") %>%
  # standardize to lowercase to ignore capitalization
  # group by the unique combination of taxa, name, desc
  group_by(pathway_Kegg) %>%
  # sum all sample columns
  summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop")

#wide_df_filtered<-subset(collapsed_df,AllSamples>=1000)

# --- Write final combined output ---
write.csv(collapsed_df,"Microbe_Genes_counts.csv")
write.csv(collapsed_df_ByGene,"Microbe_Genes_collapsed_counts.csv")


message("Done! Final combined table saved to: ", final_outfile)

PART 2, add pathway info:

# Now read in pathway data:
# Read in bacterial pathway summaries
Pathway_files <- list.files("../../AnnotationOutputs/", pattern = "*_Bacteria.csv", full.names = TRUE)

Pathway_list <- lapply(Pathway_files, function(f) {
  df <- fread(f)
  # Add sample name (remove long suffix so it's just H1, H2, etc.)
  df[, Sample := sub("_Bacteria.csv", "", basename(f))]
  return(df)
})

# Combine into one long table
Pathway_combined <- rbindlist(Pathway_list, use.names = TRUE, fill = TRUE)

sample_cols_path <- setdiff(colnames(Pathway_wide), c("Pathways"))

# Now you can pivot wider
Pathway_wide <- dcast(Pathway_combined, Pathways ~ Sample, value.var = "TotalCounts")
Pathway_wide <- Pathway_wide %>%
  mutate(AllSamples = rowSums(across(all_of(sample_cols_path)), na.rm = TRUE))


# Histogram-style binning
hist_df_Path <- Pathway_wide %>%
  count(AllSamples, name = "N")


# Plot histograms to look at distributions of microbial genes/pathways at different levels
ggplot(collapsed_df, aes(x = AllSamples)) +
  geom_histogram(bins = 75, fill = "steelblue") +
  scale_x_log10() + 
  theme_classic(base_size = 14) +
  labs(
    x = "log10(Total Genes per taxa Occurrences Across Samples)",
    y = "Number of Genes",
    title = "Distribution of Total Gene Occurrences"
  )#+scale_y_break(c(4000, 10989000))

ggplot(collapsed_df_ByGene, aes(x = AllSamples)) +
  geom_histogram(bins = 75, fill = "steelblue") +
  scale_x_log10() + 
  theme_classic(base_size = 14) +
  labs(
    x = "log10(Total Gene Occurrences Across Samples)",
    y = "Number of Genes",
    title = "Distribution of Total Pathway Occurrences"
  )

ggplot(collapsed_df_ByPathway_GO, aes(x = AllSamples)) +
  geom_histogram(bins = 75, fill = "steelblue") +
  scale_x_log10() + 
  theme_classic(base_size = 14) +
  labs(
    x = "log10(Total Pathway Occurrences Across Samples)",
    y = "Number of Pathways",
    title = "Distribution of Total Pathway Occurrences"
  )

#################################
#Calculate Fungal and Bacterial Reads

reads_summary <- wide_df %>%
  # tag rows as "bacteria" or "fungi" if present
  mutate(group = case_when(
    str_detect(tolower(taxa), "bacteria") ~ "Bacteria",
    str_detect(tolower(taxa), "fungi") ~ "Fungi",
    TRUE ~ NA_character_
  )) %>%
  # keep only those two
  filter(!is.na(group)) %>%
  # collapse by group
  group_by(group) %>%
  summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop")

#Turn EC numbers into pathways:

# Split ECs into separate rows and sum counts per EC
df_ec <- collapsed_df_ByPathway_EC %>%
  filter(pathway_EC != "-") %>%
  separate_rows(pathway_EC, sep = ",") %>%
  mutate(pathway_EC = trimws(pathway_EC)) %>%
  group_by(pathway_EC) %>%
  summarise(across(starts_with("H"), sum), .groups = "drop")

# EC → Pathway mapping (character vector)
ec2path <- keggLink("pathway", "ec")

ec2path_tbl <- tibble(
  pathway_EC = sub("ec:", "", names(ec2path)),
  pathway    = basename(ec2path))%>%
  filter(!grepl("ec", pathway))

# get unique pathway descriptions
path_ids <- unique(ec2path_tbl$pathway)

# safer to query in chunks
get_desc <- function(pid) {
  tryCatch({
    keggList(pid) %>% unname()
  }, error = function(e) NA_character_)
}

path_descs <- map_chr(path_ids, get_desc)

path_tbl <- tibble(
  pathway = path_ids,
  desc    = path_descs
)

# Join descriptions back
ec2path_tbl <- ec2path_tbl %>%
  left_join(path_tbl, by = "pathway")

# Join counts and collapse
df_final <- df_ec %>%
  left_join(ec2path_tbl, by = "pathway_EC") %>%
  filter(!is.na(desc)) %>%   # only keep map pathways
  group_by(pathway, desc) %>%
  summarise(across(starts_with("H"), sum), .groups = "drop")

write.csv(df_final,"Microbe_Pathways_collapsed_EC_counts.csv")
```
# R markdowns for additional Steps:
## CompareSeqTypes 
### Rarefy and Filter datasets (Huge difference in read counts between transcriptome and amplicon)
### Normalize Data
### Quantify Shared Taxa between sequencing types
### Compare mean and coefficients of variation between soil types

## MaizeAndSorghumTranscriptAnalysis
### Ortholog Mapping - and merging orthologous genes
### Global Transcriptomic Variation (PERMANOVA)
### Differential Expression Analysis and visulization (DESeq2) 
### Functional Enrichment and GO Annotation

## MicrobeFunctionAnalysis
### Anlysis pipeline written to accomodate gene or pathway resolution
### Principal Component Analysis (PCA)
### ALDEx2 Differential Analysis

## MaizeAndMicrobeCoexpression
### Mantel Test and PERMANOVA
### Bipartite Network analysis (complete and subset by temperature)
### Differential Correlation Between Temperatures
### Elastic Net Regression (Biomass & Heat Stress Prediction)
### Random Forest analysis (Biomass & Heat Stress Prediction)
### Combination of Outputs to select important features

# Misc code
## Make table of orthologous genes and scan list of candidate genes for orthologs: 
```
#In bash:
#Download Maize protein fasta:
wget http://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.protein.fa.gz
#Download sorghum protein fasta:
wget https://ftp.sorghumbase.org/release-9/fasta/sorghum_bicolorv5/pep/Sorghum_bicolorv5.Sb-BTX623-REFERENCE-JGI-5.1.pep.all.fa.gz

#install orthofinder (can use conda install)
#run orthofinder: #Where -f is the folder containing the protein fasta files. See orthofinder documentation for more info: https://github.com/davidemms/OrthoFinder
orthofinder -f input/ -t 12 -a 2

#import ortholog file into R
```
# Amplicon Workflows (for amplicons from DNA and RNA)
```
# Work flows are essentially the same as in the dada2 tutorial found on Ben Callahan's github: https://benjjneb.github.io/dada2/tutorial_1_8.html
# Specifics for this dataset include truncation:
FORWARD_TRUNC <- 300 # determine from quality plots
REVERSE_TRUNC <- 285 # determine from quality plots
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(FORWARD_TRUNC,REVERSE_TRUNC),
                     trimLeft=c(17, 21), maxEE=c(4,6),
                     multithread=no_of_cores,
                     matchIDs=TRUE, compress=TRUE,
                     verbose=TRUE,n=1e6)
# This version of the SILVA database: Silva_nr99_v138.2_Genus.fa.gz
# Using Atchinson's distance for Beta diversity: Code by Isabella Borrero
#normalize to relative abundance
ps.prop <- transform_sample_counts(ps.clean.p1, function(OTU) OTU/sum(OTU))
#transform using CLR to account for compositionality of the data 
ps_tr_clr <- microbiome::transform(ps.prop, "clr")
sample_data(ps_tr_clr)$Genotype_Soil <- interaction(sample_data(ps_tr_clr)$Genotype, sample_data(ps_tr_clr)$Soil)

ord.pcoa.euc <- ordinate(ps_tr_clr, method="PCoA", distance="euclidean")

#PCoA (Aitchison's Distance) Active Fraction
plot_ordination(ps_tr_clr, ord.pcoa.euc, color = "Soil") +
  geom_point(aes(shape = Temp), size = 3) +
  #facet_wrap(~Soil) +  
  theme_classic() +
  theme(text = element_text(size = 20, family = "sans")) +
  #ggtitle("PCoA (Aitchison's Distance)") +
  ylab("PC2 (10.9%variance explained)")+
  xlab("PC1 (38% variance explained)")+
  scale_color_manual(values=c("brown3","purple1"))+
  labs(shape = "Temperature")

ps_dist_matrix_clr_euc <- phyloseq::distance(ps_tr_clr, method ="euclidean")
```

# Contact
For clarification on code missing annotation contact:
* Nate Korth: njkorthATncsu.edu or nate.korthATgmail.com
* Joe Gage: jlgageATncsu.edu

Whatever parts of the pipeline you use, please remember to cite all relevant packages
