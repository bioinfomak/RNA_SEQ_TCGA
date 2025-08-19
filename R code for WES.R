# Load Rsubread library
library(Rsubread)
setwd("D:/ALS/Sequencing/NGS_Data1/2_Raw_Data")
list.files()
ref <- "hg38.fa.gz"
#ref_als<- "chr9.fa.gz"  # for specific gene/ chr
buildindex(basename="./reference_index",reference=ref)
# Define the reference genome index (must be built beforehand)
reference_genome_index <- "./reference_index"

# Define input paired-end FASTQ files
fastq_files <- list.files( pattern = "fastq.gz$", full.names = TRUE)

# Loop through each pair of FASTQ files
for (i in seq(1, length(fastq_files), by = 2)) {
  # Check if there is a next file for paired-end
  if (i < length(fastq_files)) {
    readfile1 <- fastq_files[i]       # First read in the pair
    readfile2 <- fastq_files[i+1]   # Second read in the pair
    
    # Create output BAM file name
    output_bam_file <- sub("_.*", "_PE.bam", readfile1)  # Use first read file name for output BAM
  
    
    # Align the paired-end FASTQ files to the reference genome and generate BAM file
    align(index = reference_genome_index, 
          readfile1 = readfile1, 
          readfile2 = readfile2,
          output_file = output_bam_file,
          type="dna",
          minFragLength = 50,
          sortReadsByCoordinates=T, # generate .bai file
          PE_orientation = "fr",
          useAnnotation =T,
          annot.inbuilt="hg38",
          nthreads = 16,          # Specify the number of threads to use
          output_format = "BAM") # Output format
  } else {
    warning("Uneven number of FASTQ files. Last file will be ignored.")
  }
}

print("BAM files generation completed!")




