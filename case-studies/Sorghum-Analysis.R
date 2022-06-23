
#load dependencies
library("devtools")
library("data.table")
library("dplyr")
library("tidyr")
library("vcfR")
library("ggplot2")

#set file name of the VCF file to load
file <- "freebayes_D2.filtered.vcf"

#get the file from here: https://bss1innov1nafa1poc1.blob.core.windows.net/sample-container/Data-for-github/freebayes_D2.filtered.vcf


#set the samples
HighBulk <- "D2_F2_tt" #
LowBulk <- "D2_F2_TT" #


#Specify the chromosomes to be included in the analysis (i.e., exclude smaller contigs)
Chroms <- c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")


df <-
  importFromVCF(
    file = file,
    highBulk = HighBulk,
    lowBulk = LowBulk,
    chromList = Chroms
  )

ggplot(data = df) +
  geom_histogram(aes(x = DP.HIGH + DP.LOW)) +
  xlim(0,500) +
  geom_vline(xintercept=70, colour = "red") +
  geom_vline(xintercept=300, colour = "red")


#Plot Reference allele frequency
ggplot(data = df) +
  geom_histogram(aes(x = REF_FRQ))

#Plot Genotype quality
ggplot(data = df) +
  geom_histogram(aes(x = GQ.HIGH))
ggplot(data = df) +
  geom_histogram(aes(x = GQ.LOW))


df_filt <-
  filterSNPs(
    SNPset = df,
    #refAlleleFreq = 0.20,
    minTotalDepth = 100,
    maxTotalDepth = 400,
    minSampleDepth = 30,
    #depthDifference = 100,
    minGQ = 99,
    verbose = TRUE
  )

df_filt <-
  runGprimeAnalysis(
    SNPset = df_filt,
    windowSize = 5e6,
    outlierFilter = "deltaSNP",
    filterThreshold = 0.1
  )


#Plot only a subset of Chromosomes
plotQTLStats( SNPset = df_filt,
              var = "Gprime",
              plotThreshold = TRUE,
              q = 0.01,
 #             subset=c("NC_012873.2")
              subset=c("Chr04")
 )

#Run QTLseq analysis
#bulkSize: Number of individuals in the respective pools: c(HighBulk, LowBulk)!
df_filt <-
  runQTLseqAnalysis(
    SNPset = df_filt,
    windowSize = 5e6,
    popStruc = "F2",
    bulkSize = c(45, 38),
    replications = 10000,
    intervals = c(95, 99)
  )

plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE, subset=c("Chr04"))

