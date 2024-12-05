
##BiocManager::install("ShortRead")
library(ShortRead)
library(ggplot2)
library(here)

# Specify the path to the fastq.gz file
fastq_file <- here("data/LB_1_Q_L1_R1_001_AjEj3baKF3rB.fastq.gz")  # Replace with your file path
fastq_file2 <- here("data/LB_2_Q_L1_R1_001_IWYOPvW58d60.fastq.gz")  # Replace with your file path
fastq_file3 <- here("data/LB_3_Q_L1_R1_001_qojIwuVAjLqW.fastq.gz")  # Replace with your file path
fastq_file4 <- here("data/LB_4_Q_L1_R1_001_L8ZJjX60AHsX.fastq.gz")  # Replace with your file path


# Read the fastq.gz file
reads <- readFastq(fastq_file)
reads2 <- readFastq(fastq_file2)
reads3 <- readFastq(fastq_file3)
reads4 <- readFastq(fastq_file4)

# Explore the contents of the reads
print(reads)  # Displays a summary of the reads
head(sread(reads))  # Displays the sequences
head(quality(quality(reads)))  # Displays the quality scores
print(sread(reads))

# graph the quality scores
quals = quality(reads2)
numqscores = as(quals, 'matrix')
avgscore = rowMeans(numqscores, na.rm = T)
avgscore = as.data.frame(avgscore)

ggplot(avgscore) + 
  geom_histogram(aes(x=avgscore)) +
  theme_classic()
  
