

###########################################################
#                        PACKAGES
#                    ShortRead, ggplot2
###########################################################




## ---------------------------------------------
##          Packages & Libraries
## ---------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ShortRead")

library(ShortRead)
library(ggplot2)

## Downloaded FASTQ file using SRA Toolkit 
## FASTQ id of interest SRR13764788.fastq
## File stored in E:\Linux\SRA

## ---------------------------------------------
##          Read fastq file(s)
## ---------------------------------------------

setwd("E:/NGS")

# create a sub directory
dir.create('Fastq_read')
setwd("E:/NGS/Fastq_read")

## read fastq [The fastq file was retrieved from SRA using linux command line tools]
fq = readFastq('E:/Linux/SRA/SRR13764788.fastq')
head(fq)  
 
summary(fq) # Number of reads


## ---------------------------------------------
##  Look DNA sequences & graph quality scores
## ---------------------------------------------

## Quality scores: How was the quality of DNA obtained from NGS
## Did we get unamiguous base call? Did we get high quality seq data?

reads_fq = sread(fq)
head(reads_fq)


## This gives us a quick look at both ends of the read
## We can also quickly check that majority of the 
## characters are valid nt characters, which is a good sign 

## Pull out the read length

widths_fq = reads_fq@ranges@width
widths_fq #vector of width values ~ 2-5 hundred bp

## for ggplot we must have this as a dataframe
widths_fq = as.data.frame(reads_fq@ranges@width) 
write.csv(widths_fq, "widths_fq.csv", row.names = FALSE)

## Plot the widths
ggplot(widths_fq) + 
  geom_histogram(aes(x=reads_fq@ranges@width))  ## reads_fq1@ranges@width is the column name in the df

# The majority of the reads in reads_fq1 are in the range of ~500 nt

## Let's plot the quality score

quals_fq = quality(fq)
head(quals_fq)
numqscores_fq = as(quals_fq, 'matrix') #convert these quals into numeric scores
## Calculate average quality score of each read
## Remove NA by adding na.rm = TRUE
avgscore_fq = rowMeans(numqscores_fq, na.rm = TRUE)
avgscore_fq = as.data.frame(avgscore_fq)
## Calculate average quality score of each read
write.csv(avgscore_fq, "avgscore_fq.csv", row.names = FALSE)
                            
ggplot(avgscore_fq) +
  geom_histogram(aes(x=avgscore_fq)) ## avgscore_fq is the column name

## This is the distribution with ~39 as the quality score for majority of the reads
## This score indicates pretty good quality of the sequences 




