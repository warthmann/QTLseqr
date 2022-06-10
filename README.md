
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pbglQTLseqr

QTLseqr is an R package for Bulk Segregant Analysis.

QTLseqr was developed and published by Ben N. Mansfeld and Rebecca Grumet 
and all credit should go to them. We forked the software from 
[their github repository](https://github.com/bmansfeld/QTLseqr/) 
and made minor changes to adopt it to our needs. 

### **For more detailed instructions please read the vignette of the original  [here](https://github.com/bmansfeld/QTLseqr/raw/master/vignettes/QTLseqr.pdf)**

**If you use QTLseqr in published research, please cite:**

> Mansfeld B.N. and Grumet R, QTLseqr: An R package for bulk segregant
> analysis with next-generation sequencing *The Plant Genome*
> [doi:10.3835/plantgenome2018.01.0006](https://dl.sciencesocieties.org/publications/tpg/abstracts/11/2/180006)

QTLseqR in an R implementation of established methods and statistics. 
Make sure you also cite the respective method you work with: 

QTL-seq:

> Takagi, H., Abe, A., Yoshida, K., Kosugi, S., Natsume, S., Mitsuoka,
> C., Uemura, A., Utsushi, H., Tamiru, M., Takuno, S., Innan, H., Cano,
> L. M., Kamoun, S. and Terauchi, R. (2013), QTL-seq: rapid mapping of
> quantitative trait loci in rice by whole genome resequencing of DNA
> from two bulked populations. *Plant J*, 74: 174–183.
> [doi:10.1111/tpj.12105](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.12105)

G prime:

> Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk
> Segregant Analysis Using Next Generation Sequencing. *PLOS
> Computational Biology* 7(11): e1002255.
> [doi.org/10.1371/journal.pcbi.1002255](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255)

## Abstract

Next Generation Sequencing Bulk Segregant Analysis (NGS-BSA) is
efficient in detecting quantitative trait loci (QTL). Despite the
popularity of NGS-BSA and the R statistical platform, no R packages are
currently available for NGS-BSA. We present QTLseqr, an R package for
NGS-BSA that identifies QTL using two statistical approaches: QTL-seq
and G’. These approaches use a simulation method and a tricube smoothed
G statistic, respectively, to identify and assess statistical
significance of QTL. QTLseqr, can import and filter SNP data, calculate
SNP distributions, relative allele frequencies, G’ values, and
log10(p-values), enabling identification and plotting of QTL.


**Example Analysis**

The data for below example analysis and tutorial is drawn from a study on 
cold tolerance in rice by [Yang Z et al (2013)](https://doi.org/10.1371/journal.pone.0068433)

>Yang Z, Huang D, Tang W, Zheng Y, Liang K, Cutler AJ, et al. (2013) Mapping of Quantitative Trait Loci >Underlying Cold Tolerance in Rice Seedlings via High-Throughput Sequencing of Pooled Extremes. 
>PLoS ONE 8(7): e68433. https://doi.org/10.1371/journal.pone.0068433

The raw data is available from NCBI's short read archive in BioProjet [PRJNA198759](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA198759)
Run SRR834931: Rice ET pool (385 extremely tolerant individuals) 
Run: SRR834927 Rice ES pool (430 extremely sensitive individuals)

We have downloaded the raw data, performed alignment and variant calling against
the nipponbare reference genome IRGSP-1.0 with bwa and freebayes. QTLseqr analysis
on the 
[resulting VCF](https://bss1innov1nafa1poc1.blob.core.windows.net/sample-container/Data-for-github/wGQ-Filt-freebayes~bwa~IRGSP-1.0~both-segregant_bulks~filtered-default.vcf) 
file produced by freebayes can be performed as outlined below.


# Installation

Install the pbgl version of QTLseqr from github with:

``` r
# install devtools to download and install packages from github
install.packages("devtools")

#install the dependencies
install.packages(c("data.table", "dplyr", "tidyr", "vcfR", "ggplot2"), dependencies=TRUE)

# use devtools to install QTLseqr
#devtools::install_github("bmansfeld/QTLseqr")
devtools::install_github("warthmann/QTLseqr")
```

**Note:** Apart from regular package dependencies, there are some
Bioconductor tools that we use as well, as such you will be prompted to
install support for Bioconductor, if you haven’t already. QTLseqr makes
use of C++ to make some tasks significantly faster (like counting SNPs).
Because of this, in order to install QTLseqr from github you will be
required to install some compiling tools (Rtools and Xcode, for Windows
and Mac, respectively).


# Example

This example figure is from Mansfeld et al. and serves for comparison to our/your own analysis

![Example
figure](https://github.com/bmansfeld/QTLseqr/raw/master/all_plots.png
"Example figure")

![Example
figure](https://github.com/warthmann/QTLseqr/blob/master/all_plots.png
"Norman Figure")

This is a basic example which shows you how to import and analyze
NGS-BSA data.

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


#Set sample
HighBulk <- "ET-pool-385" # SRA-sample: SRR834931
LowBulk <- "ES-pool-430" # SRA-sample: SRR834927 

#set the VCF file name
file <- "wGQ-Filt-freebayes~bwa~IRGSP-1.0~both-segregant_bulks~filtered-default.vcf"

#Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
Chroms <- c("NC_029256.1","NC_029257.1","NC_029258.1","NC_029259.1","NC_029260.1","NC_029261.1","NC_029262.1","NC_029263.1","NC_029264.1","NC_029265.1","NC_029266.1","NC_029267.1")

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


df <- 
    importFromVCF(
        file = file,
        highBulk = HighBulk,
        lowBulk = LowBulk,
        chromList = Chroms
     )


#Filter SNPs based on some criteria
df_filt <-
    filterSNPs(
        SNPset = df,
        refAlleleFreq = 0.20,
        minTotalDepth = 100,
        maxTotalDepth = 400,
        minSampleDepth = 40,
        minGQ = 99,
        depthDifference = 100,
        verbose = TRUE
    )


#Run G' analysis
df_filt <- runGprimeAnalysis(
    SNPset = df_filt,
    windowSize = 1e6,
    outlierFilter = "deltaSNP",
    filterThreshold = 0.1
)

#Run QTLseq analysis
df_filt <- runQTLseqAnalysis(
    SNPset = df_filt,
    windowSize = 1e6,
    popStruc = "F2",
    bulkSize = c(385, 430),
    replications = 10000,
    intervals = c(95, 99)
)

#Plot
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
plotQTLStats(SNPset = df_filt, var = "nSNPs", plotIntervals = TRUE)
plotQTLStats(SNPset = df_filt, var = "negLog10Pval", plotIntervals = TRUE)

#Analytics
plotGprimeDist(SNPset = df_filt, outlierFilter = "Hampel")
plotGprimeDist(SNPset = df_filt, outlierFilter = "deltaSNP", filterThreshold = 0.1)

#export summary CSV
getQTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")
```
