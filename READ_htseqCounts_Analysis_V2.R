#------points to the directory in which the htseq-count output files are located-----------#
directory = "C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/CRC_Subtypes_FS21/READ"

# specify which files to read in using list.files, and select those files which contain the string "htseq" using grep
sampleFiles = grep("htseq",list.files(directory),value=TRUE)

samplePhenotype = read.csv("C:/Users/cmk_8/OneDrive/Dokumente/Master ACLS/Master Thesis/CRC_Subtypes_FS21/READ/READ_pheno.csv", header = TRUE, sep = ",")

sampleFileName = sub(".htseq.counts","",sampleFiles)
#sampleFileName2 <- sub(".htseq.counts","",samplePhenotype)
sampleCmsGroup = samplePhenotype$CMS_final
sampleID = samplePhenotype$sample_id
sampleSubmitterID = samplePhenotype$submitter_id
sampleCaseID = samplePhenotype$case_id
sampleTissue = samplePhenotype$sample_type

sampleTable = data.frame(sampleName = sampleFileName,
                          fileName = sampleFiles,
                          condition = sampleTissue,
                          cmsGroup = sampleCmsGroup,
                          submitterID = sampleSubmitterID,
                          sampleID = sampleID,
                          caseID = sampleCaseID,
                          fileID = sampleFileName)
#fileName2 = sampleFileName2)





# Then we build the DESeqDataSet using the following function
library("DESeq2")

ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition + cmsGroup)
ddsHTSeq

# Pre-filter low count genes 
keep = rowSums(counts(ddsHTSeq)) >= 10
ddsHTSeq = ddsHTSeq[keep,]

# Store the count data
countdata = assay( ddsHTSeq )
head( countdata )
# Store the meta data
coldata = colData( ddsHTSeq )

# ddsFullcountTable
#ddsFullCountTable = DESeqDataSetFromMatrix(countData= countdata, colData = coldata, design= ~ condition + cmsGroup)
#ddsFullCountTable

########################################################################################################
## Running the DESeq2 Pipeline - Chapter 3 Beginner Guide
########################################################################################################
#----- Analyze the Normal Tissue VS Primary Tumor of same submitter ID

# 1 ) Prepare the Relevant columns of interest for analysis
dds_normal = ddsHTSeq[ , ddsHTSeq$condition == "Solid Tissue Normal"]
as.data.frame( colData(dds_normal) )

# 2 ) Running the Pipeline
dds = DESeq(ddsHTSeq)
# --> Maybe run the Pipline with exaclty the same colnames --> Change this and try again

# 3) Inspecting the results table
res = results(dds)
res


























