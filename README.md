
<!-- README.md is generated from README.Rmd. Please edit that file -->

# QTLseqr (PBGL version)

QTLseqr is an R package for Bulk Segregant Analysis from NGS data.

QTLseqr was developed and published by Ben N. Mansfeld and Rebecca Grumet 
and all credit should go to them. We forked [their github repository](https://github.com/bmansfeld/QTLseqr/) 
and made minor changes to their software to adopt it to our needs. We provide a brief example analysis below. For detailed instructions on usage and theory of the statistics please consult the [vignette of the original QTLseqr package](https://github.com/bmansfeld/QTLseqr/raw/master/vignettes/QTLseqr.pdf).

**If you use QTLseqr, please cite:**

> Mansfeld B.N. and Grumet R, QTLseqr: An R package for bulk segregant
> analysis with next-generation sequencing *The Plant Genome* (2018)
> [doi:10.3835/plantgenome2018.01.0006](https://dl.sciencesocieties.org/publications/tpg/abstracts/11/2/180006)


## Software and Methods

Bulk Segregant Analysis with Next Generation Sequencing (NGS-BSA) is an effective tool
for genetic mapping of causal loci, including QTLs. NGS-BSA is popular and
there exist several statistical approaches to the problem. Ben N. Mansfeld and Rebecca Grumet implemented two of them in R: 

* QTL-seq (deltaSNP index)
* G’ (G prime)

QTLseqr can import and filter SNP data, count relative allele frequencies, calculate SNP distributions, G’ values, log10(p-values), and the corresponding thresholds to establish  statistical significance. When using QTLseqr please make sure you also cite the publications that first described the respective method and statistic you work with: 

**QTL-seq (DeltaSNP-Index)**

> Takagi, H., Abe, A., Yoshida, K., Kosugi, S., Natsume, S., Mitsuoka,
> C., Uemura, A., Utsushi, H., Tamiru, M., Takuno, S., Innan, H., Cano,
> L. M., Kamoun, S. and Terauchi, R. (2013), QTL-seq: rapid mapping of
> quantitative trait loci in rice by whole genome resequencing of DNA
> from two bulked populations. *Plant J*, 74: 174–183.
> [doi:10.1111/tpj.12105](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.12105)

**G prime**

> Magwene PM, Willis JH, Kelly JK (2011).
> The Statistics of Bulk Segregant Analysis Using Next Generation Sequencing. *PLOS Computational Biology* 7(11): e1002255.
> [doi.org/10.1371/journal.pcbi.1002255](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255)


**Example Analysis**

The data for below example analysis and tutorial is drawn from a study on 
cold tolerance in rice by [Yang, Z et al. (2013)](https://doi.org/10.1371/journal.pone.0068433). It is the same data that Mansfeld and Grumet had used and we only slightly modified their instructions taken from their [github page](https://github.com/bmansfeld/QTLseqr/) and [vignette](https://github.com/bmansfeld/QTLseqr/raw/master/vignettes/QTLseqr.pdf).

>Yang Z, Huang D, Tang W, Zheng Y, Liang K, Cutler AJ, et al. (2013).
>Mapping of Quantitative Trait Loci Underlying Cold Tolerance in Rice Seedlings via High-Throughput Sequencing of Pooled Extremes. 
>PLoS ONE 8(7): e68433. https://doi.org/10.1371/journal.pone.0068433

The raw sequencing data for [Yang Z, et al. (2013)](https://doi.org/10.1371/journal.pone.0068433) is available from NCBI's short read archive in [BioProjet PRJNA198759](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA198759)

* Run SRR834931: Rice ET pool (385 extremely tolerant individuals) 
* Run SRR834927: Rice ES pool (430 extremely sensitive individuals)

We have downloaded the raw data (fasterq-dump SRR834931 SRR834927) and performed alignment (bwa) and variant calling (freebayes) against the Nipponbare reference genome (IRGSP-1.0) with our [PBGL snakemake workflow](https://github.com/pbgl/dna-proto-workflow). 
QTLseqr analysis on the [resulting VCF file produced by freebayes](https://bss1innov1nafa1poc1.blob.core.windows.net/sample-container/Data-for-github/wGQ-Filt-freebayes~bwa~IRGSP-1.0~both-segregant_bulks~filtered-default.vcf) can be performed as outlined below.


# Installation

Install the PBGL version of QTLseqr from github like so:

``` r
#install the dependencies
install.packages(c("data.table", 
                   "modeest",
                   "locfit",
                   "dplyr", 
                   "tidyr", 
                   "vcfR", 
                   "ggplot2"), dependencies=TRUE)

# install devtools
install.packages("devtools", dependencies=TRUE)

# use devtools to install QTLseqr from github
library("devtools")
#devtools::install_github("bmansfeld/QTLseqr")
#devtools::install_github("warthmann/QTLseqr")
devtools::install_github("pbgl/QTLseqr")

```
**Note:** When installing above R-packages, you may run into missing system dependencies. On clean Ubuntu systems these may include: libssl-dev liblapack-dev libopenblas-dev libcurl4-openssl-dev libxml2-dev libudunits2-dev libmariadbclient-dev gfortran libpng-dev libjpeg-dev.
Some functions of QTLseqr (e.g., counting SNPs) are implemented in C++. Hence, the install of QTLseqr from github requires the compiler. On Linux this should work out of the box, for Windows and Mac you will need Rtools or Xcode, respectively. You might need assistance from your system administrator.

Do the install from github if you can. Alternatively, we provide binaries of the pbgl version of the QTLseqr R-package for local install: [Linux](https://bss1innov1nafa1poc1.blob.core.windows.net/sample-container/Data-for-github/pbglQTLseqr-binaries/QTLseqr_0.7.5.2_R_x86_64-pc-linux-gnu.tar.gz), 
[Mac](https://bss1innov1nafa1poc1.blob.core.windows.net/sample-container/Data-for-github/pbglQTLseqr-binaries/QTLseqr_0.7.5.2.tgz), 
[Windows](https://bss1innov1nafa1poc1.blob.core.windows.net/sample-container/Data-for-github/pbglQTLseqr-binaries/QTLseqr_0.7.5.2.zip). They were compiled with R version 4.2.


# Example Analysis

**Expected Results**

This example figure is taken from Mansfeld et al. and serves for comparison to our/your own analysis.

![Example
figure](https://github.com/warthmann/QTLseqr/blob/master/all_plots.png
"Mansfeld Figure")

**Step-by-step instructions**

The below is functional R code. Copy/Paste into R-Studio.

``` r

#load dependencies
library("devtools")
library("data.table")
library("dplyr")
library("tidyr")
library("vcfR")
library("ggplot2")

#load the package
library("QTLseqr")

#Set samples
HighBulk <- "ET-pool-385" # SRA-run: SRR834931
LowBulk <- "ES-pool-430" # SRA-run: SRR834927 

#set file name of the VCF file to load
file <- "wGQ-Filt-freebayes~bwa~IRGSP-1.0~both-segregant_bulks~filtered-default.vcf"

#Specify which chromosomes should be included in the analysis (i.e., exclude smaller contigs)
Chroms <- c("NC_029256.1",
            "NC_029257.1",
            "NC_029258.1",
            "NC_029259.1",
            "NC_029260.1",
            "NC_029261.1",
            "NC_029262.1",
            "NC_029263.1",
            "NC_029264.1",
            "NC_029265.1",
            "NC_029266.1",
            "NC_029267.1")
```

Chromosome information can be found in the VCF file header.
Open the VCF file in a text editor and find the '##contig' lines. E.g.:

```console
$ less -S wGQ-Filt-freebayes~bwa~IRGSP-1.0~both-segregant_bulks~filtered-default.vcf
[...]
##reference=genomes_and_annotations/IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna 
##contig=<ID=NC_029256.1,length=43270923>  
##contig=<ID=NC_029257.1,length=35937250>  
##contig=<ID=NC_029258.1,length=36413819>  
##contig=<ID=NC_029259.1,length=35502694>  
##contig=<ID=NC_029260.1,length=29958434>  
##contig=<ID=NC_029261.1,length=31248787>  
##contig=<ID=NC_029262.1,length=29697621>  
##contig=<ID=NC_029263.1,length=28443022>  
##contig=<ID=NC_029264.1,length=23012720>  
##contig=<ID=NC_029265.1,length=23207287>  
##contig=<ID=NC_029266.1,length=29021106>  
##contig=<ID=NC_029267.1,length=27531856> 
[...]
```


``` r
#Import the data from the VCF file
#also calculates SNPindex and DeltaSNP
df <- 
    importFromVCF(
        file = file,
        highBulk = HighBulk,
        lowBulk = LowBulk,
        chromList = Chroms
     )

#Filter SNPs based on user-specified thresholds. 
#Remove the GQ filter if GQ is not present in the FORMAT field of your VCF file
df_filt <-
    filterSNPs(
        SNPset = df,
        refAlleleFreq = 0.20,
        minTotalDepth = 70,
        maxTotalDepth = 200,
        minSampleDepth = 30,
        #depthDifference = 100,
        minGQ = 150,
        verbose = TRUE
    )

#Run G' analysis on the filtered SNPs
df_filt <- 
    runGprimeAnalysis(
        SNPset = df_filt,
        windowSize = 1e6,
        outlierFilter = "deltaSNP",
        filterThreshold = 0.1
    )

#Run QTLseq analysis
#bulkSize: Number of individuals in the respective pools: c(HighBulk, LowBulk)!
df_filt <- 
    runQTLseqAnalysis(
        SNPset = df_filt,
        windowSize = 1e6,
        popStruc = "F2",
        bulkSize = c(385, 430), 
        replications = 10000,
        intervals = c(95, 99)
    )

#Plot SNP density
plotQTLStats(SNPset = df_filt, var = "nSNPs", plotIntervals = TRUE)

#Plot the statistics
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.02)
plotQTLStats(SNPset = df_filt, var = "negLog10Pval", plotThreshold = TRUE)

#Plot only a subset of Chromosomes
plotQTLStats( SNPset = df_filt, 
              var = "Gprime", 
              plotThreshold = TRUE, 
              q = 0.02, 
              subset=c("NC_029256.1","NC_029263.1")
)

?plotQTLStats

#get a list of significant Loci
getQTLTable(SNPset = df_filt, method = "Gprime", alpha=0.02)
getQTLTable(SNPset = df_filt, method = "QTLseq")

?getQTLTable

#export table as csv file
getQTLTable(SNPset = df_filt, 
            method = "Gprime", 
            alpha=0.01, 
            export=TRUE, 
            fileName= "my_first_BSA_result.csv"
)

```

**Determine Filtering Thresholds**

The below is functional R code. Copy/Paste into R-Studio. For details please consult the [vignette of the original QTLseqr package](https://github.com/bmansfeld/QTLseqr/raw/master/vignettes/QTLseqr.pdf).

``` r
library("ggplot2")

#Plot read depth(s)
ggplot(data = df) +
        geom_histogram(aes(x = DP.HIGH + DP.LOW)) +
        xlim(0,500)

ggplot(data = df) +
        geom_histogram(aes(x = DP.HIGH + DP.LOW)) +
        xlim(0,500) + 
        geom_vline(xintercept=70, colour = "red") +
        geom_vline(xintercept=200, colour = "red") 


#Plot Reference allele frequency 
ggplot(data = df) +
        geom_histogram(aes(x = REF_FRQ))


#Plot Genotype Qualities 
ggplot(data = df) +
        geom_histogram(aes(x = GQ.HIGH))
ggplot(data = df) +
        geom_histogram(aes(x = GQ.LOW))
        

#Plot SNP-Index
#When plotting the SNP-Index of segregating bulks of an F2 population we expect most of the SNPs approximiately normally distributed arround 0.5 and peaks in the tails.
ggplot(data = df) +
        geom_histogram(aes(x = SNPindex.HIGH))
ggplot(data = df) +
        geom_histogram(aes(x = SNPindex.LOW))

#Confirm the G' distribution to be log-normal 
plotGprimeDist(SNPset = df_filt, outlierFilter = "Hampel")
plotGprimeDist(SNPset = df_filt, outlierFilter = "deltaSNP", filterThreshold = 0.1)

``` 
        
