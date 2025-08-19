library(VariantAnnotation)
library(VennDiagram)
library(UpSetR)
library(ggvenn)
set01 <- read.table("C:/Users/S Sarkar/Desktop/TGS ALS/VR1_S9_output_filterd.vcf", sep= "")
set2 <- read.table("C:/Users/S Sarkar/Desktop/TGS ALS/VR2_S10_output_filterd.vcf", sep="")
set3 <- read.table("C:/Users/S Sarkar/Desktop/TGS ALS/VR3_S11_output_filterd.vcf", sep="")
set4 <- read.table("C:/Users/S Sarkar/Desktop/TGS ALS/VR4_S12_output_filterd.vcf", sep="")
set5 <- read.table("C:/Users/S Sarkar/Desktop/TGS ALS/VR5_S13_output_filterd.vcf", sep="")
set6 <- read.table("C:/Users/S Sarkar/Desktop/TGS ALS/VR6_S14_output_filterd.vcf", sep="")
set7 <- read.table("C:/Users/S Sarkar/Desktop/TGS ALS/VR7_S15_output_filterd.vcf", sep="")
set8 <- read.table("C:/Users/S Sarkar/Desktop/TGS ALS/VR8_S16_output_filterd.vcf", sep="")


set1 <- as.vector(set1$V3)
set2 <- as.vector(set2$v3)
set3 <- as.vector(set3$V3)
set4 <- as.vector(set4$V3)
set5 <- as.vector(set5$V3)
set6 <- as.vector(set6$V3)
set7 <- as.vector(set7$V3)
set8 <- as.vector(set8$V3)

read_sets = list(set1_reads = set1,
                 set2_reads = set2,
                 set3_reads = set3,
                 set4_reads = set4,
                 set5_reads = set5,
                 set6_reads = set6,
                 set7_reads = set7,
                 set8_reads = set8)

ggvenn(read_sets)

x = list(read_sets)

x = list(S01=AS1_S1_output_filterd$V3, S02=AS2_S2_output_filterd$V3,
         S03=AS3_S3_output_filterd$V3,S04= AS4_S4_output_filterd$V3,
         S05=AS5_S5_output_filterd$V3, S06=AS6_S6_output_filterd$V3,
         S07=AS7_S7_output_filterd$V3,S08=AS8_S8_output_filterd$V3,
         S17=AC1_S17_output_filterd$V3, S18=AC2_S18_output_filterd$V3,
         S19=AC3_S19_output_filterd$V3,S20=AC4_S20_output_filterd$V3,
         S21=AC5_S21_output_filterd$V3, S22=AC6_S22_output_filterd$V3,
         S24=AC8_S24_output_filterd$V3
)

ggVennDiagram(x[c(1:8,15)])  # 4d venn




y<- process_region_data(Venn(x[c(1:8,15)])) 
y.common<- y[y$name=="S01/S02/S03/S04/S05/S06/S07/S08",]
y.common$item
y.common$item




upset(fromList(read_sets),
      sets = c("set1_reads", "set2_reads", "set3_reads", "set4_reads", "set5_reads", "set6_reads", "set7_reads", "set8_reads"),
      number.angles = 20, point.size = 2.5, line.size = 1.5,
      mainbar.y.label = "read intersection", sets.x.label = "read set size",
      text.scale = c(1.5, 1.5, 1.25, 1.25, 1.5, 1.5), mb.ratio = c(0.65, 0.35),
      order.by = "freq", keep.order = TRUE)

upset(
  fromList(read_sets),
  sets = c("set1_reads", "set2_reads", "set3_reads", "set4_reads", "set5_reads", "set6_reads", "set7_reads", "set8_reads"),
  number.angles = 20,
  point.size = 2.5,
  line.size = 1.5,
  mainbar.y.label = "Read Intersections",
  sets.x.label = "read set size",
  sets.bar.color = "skyblue",
  matrix.color = "darkred", order.by = "freq", keep.order = TRUE)

