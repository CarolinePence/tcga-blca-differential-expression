#### Section 1: Set options and environment ####

## Create an empty list object to store options inside of to keep options organized
opt <- list()

## Working directory for the project
# Relative path so the pipeline runs from wherever the repo is cloned
opt$wd <- "."
setwd(opt$wd)

## Create a character value defining the version of R and make it more path-friendly by removing special characters
opt$rString <- gsub(" ","_",gsub(" \\(.*","",R.version.string))

## Create a path to install the packages into
# Keeping packages in a project-local, R-version-specific folder avoids
# permission problems and stops this project's package versions from
# clashing with anything else installed on the machine
opt$yourPckDir <- file.path("./Rpackages", opt$rString)
dir.create(opt$yourPckDir, recursive = TRUE, showWarnings = FALSE)

## Tell R to look in that folder first when loading and installing packages
.libPaths(c(opt$yourPckDir, .libPaths()))

## Specify packages to be installed from the CRAN server (the typical repository of R packages)
opt$cranPackages <- c("BiocManager","ggplot2","ggrepel")

## Specify packages to be installed from the Bioconductor server (a biology specialized repository of packages)
opt$biocPackages <- c("TCGAbiolinks","DESeq2","SummarizedExperiment")

## Install any of the above that are not already available
opt$installed <- rownames(installed.packages())

opt$missingCran <- setdiff(opt$cranPackages, opt$installed)
if (length(opt$missingCran)) {
    install.packages(opt$missingCran, lib = opt$yourPckDir)
}

opt$missingBioc <- setdiff(opt$biocPackages, opt$installed)
if (length(opt$missingBioc)) {
    BiocManager::install(opt$missingBioc, lib = opt$yourPckDir,
                         ask = FALSE, update = FALSE)
}

# Make sure the package can be loaded
library(TCGAbiolinks)



#### Section 2: Load sample sheets and check quality ####
## Structure sample sheets files into data.frame objects
# The group data.frames result from searching the GDC and exporting "Sample sheets"
group1 <- read.csv("BLCA_non_smokers.csv", header = TRUE)
group2 <- read.csv("BLCA_reformed_smokers.csv", header = TRUE)


# The purpose of this section is to check for certain issues with the sample files.

## Check that group1 and group2 have the correct class
if( !"data.frame"%in%class(group1) ){ stop("Oh, no! Your group1 object is not a data.frame. Go fix this!") }
if( !"data.frame"%in%class(group2) ){ stop("Oh, no! Your group2 object is not a data.frame. Go fix this!") }


## Check that the needed column names are present in the sample sheet objects.
opt$checkedColumns <- c("Project.ID","Sample.ID","Data.Type")
if( !all(opt$checkedColumns%in%colnames(group1)) ){ stop("Oh, no! You are missing important column names. Check your input file / code for group 1.") }
if( !all(opt$checkedColumns%in%colnames(group2)) ){ stop("Oh, no! You are missing important column names. Check your input file / code for group 2.") }

## Check for repeated samples/cases that could distort sample comparisons
# This is needed when there are multiple files describing the same sample
if( length(unique(group1$Sample.ID)) < nrow(group1) ){stop("Oh, no! You have some duplicate values in the Sample.ID column of group 1; investigate and correct this. Usually this is due to ffpe samples which should be prioritized for removal.")}
if( length(unique(group2$Sample.ID)) < nrow(group2) ){stop("Oh, no! You have some duplicate values in the Sample.ID column of group 2; investigate and correct this. Usually this is due to ffpe samples which should be prioritized for removal.")}
if( length(unique(group1$Case.ID  )) < nrow(group1) ){stop("Oh, no! You have some duplicate values in the Case.ID column of group 1. This is nearly always in error; investigate and correct this if needed. The reason(s) for multiple files per case being present are typicallly described in the case ID pages on the GDC.")}
if( length(unique(group2$Case.ID  )) < nrow(group2) ){stop("Oh, no! You have some duplicate values in the Case.ID column of group 2. This is nearly always in error; investigate and correct this if needed. The reason(s) for multiple files per case being present are typicallly described in the case ID pages on the GDC.")}


## Check if there are so many samples that the user will possibly run into memory issues
if((nrow(group1)+nrow(group2))>250 ){
    warning("You have too many samples and might run into computational issues. ",
         "Strongly consider subsetting your sample sheet!")
}else if(max(c(nrow(group1),nrow(group2)))>100){
    warning("You have a large number of samples in one of your groups. ",
         "Consider randomly subsetting your sample sheets to to reduce their size ",
         "and avoid calculation issues.")
}

## Check if there are too few samples in either group.
# Considering how poorly controlled these experiments are (e.g. each case is a different patient),
#     there is little expectation of statistical reliability at <6 cases per group.
if(min(c(nrow(group1),nrow(group2)))<6){
    stop("Fewer than 6 samples in a group — statistical power is inadequate; revise sample selection.")
}

## Check that the data.type needed for differential expression analysis is described
if(!all(c(group1$Data.Type,group2$Data.Type)=="Gene Expression Quantification")){
    stop("This pipeline is only built to use Gene Expression Quantification data. ")
}




#### Section 3: Search for data to download ####
##### Search GDC for all files associated with target project #####
## Load the primary library used to search the GDC dataset
# This should have already been run in previous steps, but is included here for debugging.
library(TCGAbiolinks)

## Create a character vector that determines what projects to search through on the GDC based on the grouping files
opt$project <- unique(c(group1$Project.ID,group2$Project.ID))

## Do a broad search of the GDC that returns all relevant files
query1 <- GDCquery(
    project = opt$project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    access = "open"
)

##### Investigate query results to identify targeted samples ####
## Extract the data.frame from inside the initial
samDf <- query1$results[[1]]

## Returns TRUE if each row is in at least one of the groups
samDf$inGrp1      <- samDf$sample.submitter_id%in%group1$Sample.ID
samDf$inGrp2      <- samDf$sample.submitter_id%in%group2$Sample.ID
samDf$sampleLogic <- samDf$inGrp1 | samDf$inGrp2

## Visualize and explore your search results to double check it's working.
# Edit sample files to remove problematic samples if needed.
# Use the View() function to visualize directly:
if (interactive()) View(samDf[order(-samDf$sampleLogic),])

## Save the desired sample "barcodes" (a.k.a. cases) to use to filter another query of the database
opt$initBarcodes <- samDf$sample.submitter_id[samDf$sampleLogic]

##### Search GDC again but limit results to targeted samples ####
## Search the GDC and filter the results by the file barcodes
query2 <- GDCquery(
    project = opt$project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    access = "open",
    barcode = opt$initBarcodes #Note that this is the change from the first search.
)


##### Check group and query results #####
# The purpose of this section is to check for certain issues with the search results.

## Check to make sure grouping files have at least some data in them
length(group1$Sample.ID) #Should be >0
length(group2$Sample.ID) #Should be >0

## Check that the grouping files do not overlap:
sum(group1$Sample.ID%in%group2$Sample.ID) #Should be 0
sum(group2$Sample.ID%in%group1$Sample.ID) #Should be 0

## Check that the second search had some results
samDf2 <- query2$results[[1]]
nrow(samDf2) #Should be >0

## Determine which search results correspond to the samples in the group files
samDf2$inGroup1 <- samDf2$sample.submitter_id%in%group1$Sample.ID
samDf2$inGroup2 <- samDf2$sample.submitter_id%in%group2$Sample.ID
sum(samDf2$inGroup1) #Should be at least 6 per group
sum(samDf2$inGroup2) #Should be at least 6 per group

## Check that the second search had results that overlapped with BOTH group files
all(samDf2$inGroup1|samDf2$inGroup2)

## Check for any repeated patients
# Unless this project is modified to consider case, this should be FALSE.
any(duplicated(samDf2$cases.submitter_id))

## Visualize the different sample types in the search results
# This sometimes helps explain why there are duplicates
table(samDf2$sample_type)






#### Section 4: Obtain and structure the sequencing data into an R object ####
##### Download and prepare data into an object ####
## Download the data files from the GDC
GDCdownload(query = query2, method = "api", files.per.chunk = 10)

## Read the downloaded files and structure the sequencing data into a single object
dds1 <- GDCprepare(query = query2)

##### Check the quality of the prepared data ####
## Load the packages for modifying and analyzing the gene expression data structure
library(SummarizedExperiment) #Needed to access the function colData()
library(DESeq2) #Does the actual count-based analyses including Differential Expression (DE) analysis.

## Visualize the clinical data describing the samples to double check the preparation worked.
# This information is stored inside of a complex object that requires special functions to pull out.
# This function is named colData() because it described the "columns" of the count data (the samples).
# Once the data isolated, I have convert it to a data.frame so that it can be used in View()
if (interactive()) View(as.data.frame(colData(dds1)))

## See which sample IDs are duplicated (if any)
which(duplicated(dds1$sample_submitter_id)) #Should be none; remove with: dds1 <- dds1[,!duplicated(dds1$sample_submitter_id)]


## Additional notes about SummarizedExperiment based objects!
# These are multi-dimensional data set with multiple sections:
#   - colData   #This section houses a table describing each sample, accessed via $ or colData(dds)
#   - rowRanges #This section houses a table describing each locus (typically genes), accessed via rowRanges(dds)
#   - assays    #This is a list of data.frames describing the counts per locus (rows) and sample (columns) combination, the raw count data can be accessed via assays(dds)[[1]].
#
# The extraction tools change behavior for this data:
#   - $ and [[]] are used to access / modify / and create sample descriptors. I recommend $.
#   - [,] is used to subset the loci (using the rows argument) or the samples (using the column argument)
#
# Documentation: https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html

## Add new logical columns to the clinical data for use in downstream data cleaning
dds1$group1 <- dds1$sample_submitter_id%in%group1$Sample.ID
dds1$group2 <- dds1$sample_submitter_id%in%group2$Sample.ID

## Add a factor column to the clinical data to control how DESeq2 will group the samples
# This code ends up defining the log 2 fold change calculation as group 2 - group 1
# This definition ends up being really important when trying to interpret the results
dds1$comp   <- factor(x = dds1$group2,levels = c(FALSE,TRUE))

## Convert the sequencing data into DESeqDataSet class (a subclass of RangedSummarizedExperiment)
# This is also where we define what data DESeq2 should use to control the grouping.
dds1 <- DESeqDataSet(dds1,design = ~comp)

## Visualize the samples' clinical data and ensure they are what was expected
if (interactive()) View(as.data.frame(colData(dds1)))

## Some various checks of case descriptors
# These are are included to help understand / check the data
table(dds1$project_id)
table(dds1$sample_type)
table(dds1$ajcc_pathologic_stage)     #Currently broken in Spring 2026!
table(dds1$tissue_or_organ_of_origin) #Currently broken in Spring 2026!




#### Section 5: Analyze your count data #####
##### Detect bad samples based on PCA of normalized counts #####
## Create a function to visualize overall clustering of the samples via PCA
pcaFun <- function(ddsFun,plotToggle=T){
    ### Check arguments of pcaFun prior to analysis
    if(!"DESeq2"%in%installed.packages()){
        stop("DESeq2 needs to be available. Run the configuration file.")
    }
    if(class(ddsFun)=="DESeq2"){
      stop("ddsFun needs to be a DESeq2 class object")
    }
    if(!"comp"%in%colnames(colData(ddsFun))){
      stop("The column 'comp' needs to be the name of a factor column in ddsFun's sample descriptor section")
    }

    ### Calculate and extract normalized values
    ddsFun <- DESeq2::estimateSizeFactors(ddsFun)
    normCounts <- as.data.frame(counts(ddsFun, normalized = T))

    ### Calculate the variance of counts across each gene in order to subset data
    # This is needed because only rows with > 0 variance can be used in PCA
    # The apply loops function var() (variance) across rows of normCounts
    perLocusVar   <- apply(normCounts,MARGIN = 1,var)

    ### Determine which genes to use based on variance
    pcaLocusLogic <- perLocusVar>median(perLocusVar[perLocusVar>0])
    if(sum(pcaLocusLogic)==0){stop("Your data did not contain loci with variable counts")}

    ### Subsets and transposes the dataset for analysis with prcomp
    pcaInput <- t(normCounts[which(pcaLocusLogic),])

    ### Calculate principal components with scaling (to make each locus equal in weight)
    pca <- prcomp(pcaInput,scale. = T)

    ### Calculate the percent variance explained by each principal component
    pca_var <- pca$sdev^2
    pve <- pca_var / sum(pca_var)

    ### Plots a scatterplot of the samples visualized via principal component analysis
    if(plotToggle){
        plot(
            x = pca$x[,1],y = pca$x[,2],asp = 1,
            xlab = paste0("PC1: ",round(pve[1],3)*100,"%"),
            ylab = paste0("PC2: ",round(pve[2],3)*100,"%"),
            col  = colData(ddsFun)$comp,
            pch  = as.numeric(colData(ddsFun)$comp)+1
        )
    }

    ### Return a list of the data
    outLs <- list(x=pca$x,pve=pve)
    return(outLs)
}

## Calculate and visualize how the samples cluster via PCA
pca1Out <- pcaFun(ddsFun = dds1)


## Analyze the PCA results to determine undesirable samples / barcodes
# The goal is to remove samples that are so different from the rest, that they might increase error

## Subset the data to certain samples based on the SummarizedExperiment documentation
# Even if end up not needing to remove any samples,
#    should still create a dds2 object to keep from needing to modify
#    a bunch of downstream code.
# no extreme outliers so dds2 is just a copy of dds1
dds2 <- dds1
dim(dds2) 

## After initial sample filtering, repeat PCA to see if additional filtering is needed
# The following code looks good as the samples are generally near each other.
pca2Out <- pcaFun(dds2)


##### Detect bad loci based having high frequency of 0 in either grouping#####
## Calculate proportion with zero raw counts per grouping
locusCounts <- as.data.frame(counts(dds2, normalized = F))
locusCounts_grp1At0 <- rowMeans(locusCounts[,dds2$group1]==0)
locusCounts_grp2At0 <- rowMeans(locusCounts[,dds2$group2]==0)

## Visualize the proportion of the samples with 0 read depth per grouping
# The cutoff of 0.9 is somewhat arbitrary.
# The main problem is that loci with 0 in one group become Infinitely significant
# However, loci that are simply extremely low in one group also cause a problem!
# Such loci require modified methods to deal with.
hist(x = c(locusCounts_grp1At0,locusCounts_grp2At0),100)
opt$maxPct0Cutoff  <- 0.9
abline(v = opt$maxPct0Cutoff,col="red")

## Determine which loci meet the read depth based filter across both groups
locusCounts_logic <-
    locusCounts_grp1At0 < opt$maxPct0Cutoff &
    locusCounts_grp2At0 < opt$maxPct0Cutoff

## Subset the data to certain loci based on the SummarizedExperiment documentation
dds <- dds2[locusCounts_logic,]
dim(dds)



##### Check subset count data #####
## Calculate total gene and sample count
nrow(dds) #Expect well over 1,000 genes to remain
ncol(dds) #The number of samples should be at least 12, probably less than ~100, and definitely less than 200.

## Calculate the number of samples per grouping
class(dds$comp) #This needs to be a factor type object
table(dds$comp) #Should have two values, both more than 6.

## Calculate the minimum sample raw depth
min(colMeans(assays(dds)[[1]])) #Needs to be above 0

## Calculate the maximum proportion of samples with 0 counts per group
max(rowMeans(assays(dds)[[1]][,dds$group1]==0)) #Needs to be less than 100%
max(rowMeans(assays(dds)[[1]][,dds$group2]==0)) #Needs to be less than 100%





#### Section 6: Analyze the data with DESeq2 #####
##### Do the differential expression analysis ####

## Repeat the conversion to DESeqDataSet class object to make sure it's still well formatted
# This will produce messages if the object is no longer well structured
dds <- DESeqDataSet(dds,design = ~comp)

## Recalculate the size factors used in depth normalization after all the previous filtering
# This is needed because we have significantly changed the loci used
dds <- estimateSizeFactors(dds)

## Run the actual differential expression calculations.
# In summary:
# - It starts off by calculating size factors to do a TPM-like calculation and normalize the count data.
# - Then it calculates various measures of "dispersion" (the distributions of read depth) intended to help it interpret variability in expression across different distributions.
# - It then fits a negative binomial model to the count data (a model that better describes how read depths are distributed than something like the normal distribution)
# - It tests this model and checks for outliers.
# - Finally, it adjusts for the outliers, recalculates dispersion, refits the model, and performs final calculations.
# These measurements are incorporated into the input data and a modified version of the input is returned.
dds <- DESeq(dds)

## Organize the results of the differential expression calculations into a table.
res <- results(dds)

## Combine the results of differential expression with descriptions of the loci
resOutput <- cbind(as.data.frame(res),as.data.frame(rowRanges(dds)))

## Remove a duplicated column that causes problems downstream otherwise
resOutput <- resOutput[,!duplicated(colnames(resOutput))]

##### Analyze the results ####

## Plot a simple scatterplot / volcano plot of the analyzed genes
plot(resOutput$log2FoldChange,-log10(resOutput$padj))

## Save an image at the end of the pipeline for easy testing.
# This produces a file that can be loaded in order to recreate
#   the objects that existed when the  save.image() was run.
save.image("EndOfPipeline.rimage") 
#load("EndOfPipeline.rimage")


### Supplementary files
## Samples.csv
# create csv file of samples after filtering
samplesOut <- as.data.frame(colData(dds))

samplesOut <- samplesOut[ ,c("barcode", "group2", "project_id", "sample_type",
              "tissue_or_organ_of_origin", "primary_diagnosis",
              "age_at_diagnosis", "ajcc_pathologic_stage", "gender",
              "pack_years_smoked", "tobacco_smoking_quit_year",
              "tobacco_smoking_status", "paper_Tobacco.smoking.history",
              "paper_Number.pack.years.smoked")]

write.csv(samplesOut, "Samples.csv")


## Loci.csv
# create csv file of genes analyzed after filtering
write.csv(resOutput, "Loci.csv")



###### End of Pipeline ########

## Follow-up code analyzing the results.
source("functions.R")
library(ggplot2)
library(ggrepel)
nicePCA(pca1Out,dds1$group2)   # My first PCA plot before outlier removal
nicePCA(pca2Out,dds2$group2)    # My second PCA plot after outlier removal
niceVolcano(resOutput,0.05,1.5)
resSummary(resOutput) #This produces several message() and print() outputs describing the values that should go into my results section.


