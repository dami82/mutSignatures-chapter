# Load libs
#
library(mutSignatures)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(kableExtra)

# Set working dir
#
setwd('~/mutSignatures_chapter/03_method/')

# import mutation data from a list of VCF files
#
vcf_ext <- 'vcf.gz'
input_fld <- '../02_build_input_files/VCF/'
my_vcfs <- grep(pattern = vcf_ext, dir(input_fld), value = TRUE)
my_vcfs <- paste0(input_fld, my_vcfs)
head(my_vcfs)

# Read a list of VCF files via importVCFfiles()
#
x <- mutSignatures::importVCFfiles(vcfFiles = my_vcfs)
nrow(x)
head(x)

# Pre-processing: filter, attach nucleotide context, count mutation types
#

# Filter non SNV
x <- filterSNV(dataSet = x,  seq_colNames = c("REF", "ALT"))

# Optional: additional filtering (eg. minimum QUALITY of 10)
# Optional: adjust SAMPLE identifiers
x <- x %>%
  dplyr::filter(QUAL >= 10) %>%
  mutate(SAMPLEID = sub('__simulated$', '', .$SAMPLEID)) %>%
  mutate(SAMPLEID = gsub('[[:punct:]]', '-', toupper(.$SAMPLEID)))

nrow(x)
head(x)

# Visualize head
head(x) %>% kable() %>% kable_styling(bootstrap_options = "striped")

# Assure that the sequence Names match
# If seqNames do NOT match, an error will be returned
head(BSgenome::seqnames(BSgenome.Hsapiens.UCSC.hg19))

# Attach 3-nt context
x <- mutSignatures::attachContext(mutData = x, 
                                  BSGenomeDb = BSgenome.Hsapiens.UCSC.hg19, 
                                  chr_colName = 'CHROM', start_colName = 'POS', 
                                  end_colName = 'POS', nucl_contextN = 3, 
                                  context_colName = 'context')

# Remove mismatches (positions whose VCF reference does not match the genome sequence)
x <- removeMismatchMut(mutData = x,                  # input data.frame
                       refMut_colName = "REF",       # column name for ref base
                       context_colName = "context",  # column name for context
                       refMut_format = "N")    

# Compute mutType (according to the standards proposed by the Sanger Institute, 
# with pyrimidines as reference base)
x <- attachMutType(mutData = x,                      # as above
                   ref_colName = "REF",              # column name for ref base
                   var_colName = "ALT",              # column name for mut base
                   context_colName = "context") 
head(x)

# remove blacklist
#x <- x[!x$SAMPLEID %in% c("TCGA-A2-A0T6-01", "TCGA-E2-A15I-01"), ]
blkl <- c('TCGA-A8-A08H-01', 'TCGA-BH-A5J0-01', 'TCGA-A2-A0T6-01', 'TCGA-EW-A1OV-01', 'TCGA-E2-A2P6-01', 'TCGA-AO-A0J9-01', 'TCGA-A2-A0SY-01', 
          'TCGA-BH-A18J-01', 'TCGA-A8-A08H-01', 'TCGA-E2-A2P6-01', 'TCGA-A8-A0A7-01', 'TCGA-BH-A5J0-01', 'TCGA-EW-A1OV-01', 'TCGA-AO-A0J9-01', 
          'TCGA-A2-A0T5-01', 'TCGA-A8-A0A6-01', 'TCGA-D8-A1XK-01', 'TCGA-E9-A22G-01', 'TCGA-A8-A07L-01', 'TCGA-AN-A0AT-01', 'TCGA-AR-A0U1-01', 
          'TCGA-D8-A1JA-01', 'TCGA-A1-A0SI-01', 'TCGA-AN-A0XW-01', 'TCGA-AO-A03M-01', 'TCGA-C8-A1HM-01', 'TCGA-A2-A0SY-01', 'TCGA-BH-A18J-01',
          'TCGA-A2-A0T6-01', 'TCGA-E2-A15I-01')


head(blkl)
head(unique(x$SAMPLEID))

x <- x[x$SAMPLEID %in% blkl, ]

length(unique(x$SAMPLEID))


# Now the data.frame includes a column informing about the mutType and one about the sample identifier.
# We shall count mutation types per sample.

# Count mutation types per sample
#


my_counts <- mutSignatures::countMutTypes(mutTable = x, 
                                          mutType_colName = 'mutType', 
                                          sample_colName = 'SAMPLEID')

# my_counts is a mutationCounts object (class in mutSignatures)
# and is one of the inputs of the analysis
my_counts


# Define parameters for the non-negative matrix factorization procedure.
# you should parallelize if possible
my_params <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = 5,       # num signatures to extract
    num_totIterations = 50,           # bootstrapping: usually 500-1000
    num_parallelCores = 1)            # number of cores to use 

# De Novo Extraction of Mutational Signatures (may take a while)
#
blca_analysis <- 
  decipherMutationalProcesses(input = my_counts,
                              params = my_params)

# Inspect Results
#

# Visualize
msigPlot(blca_analysis$Results$signatures, signature = 1)
msigPlot(blca_analysis$Results$signatures, signature = 2)
msigPlot(blca_analysis$Results$signatures, signature = 3)
msigPlot(blca_analysis$Results$signatures, signature = 4)
msigPlot(blca_analysis$Results$signatures, signature = 5)

matchSignatures(blca_analysis$Results$signatures)

# Cosmic
cosmix <- getCosmicSignatures()
blca_cosmix <- cosmix[c(1,2,5,13)]
msigPlot(blca_cosmix, signature = 1)


# Check
zz <- mutSignatures::resolveMutSignatures(mutCountData = my_counts, signFreqData = blca_cosmix, byFreq = TRUE)
ww <- t(mutSignatures::as.data.frame(x = zz$results$freq.result))

head(sort(ww[, 1], decreasing = T), 7)

head(sort(ww[, 2], decreasing = T), 7)

head(sort(ww[, 3], decreasing = T), 7)

head(sort(ww[, 4], decreasing = T), 7)








matchSignatures(mutSign = blca_analysis$Results$signatures)

msigPlot(blca_analysis$Results$exposures)

zz <- t(mutSignatures::as.data.frame(x = blca_analysis$Results$exposures))
ww <- zz / apply(zz, 1, sum) 
apply(ww, 1, sum)

head(ww[order(ww[, 2], decreasing = TRUE), ])

zz <- cbind(A=1:10, B=c(1:10)/4, C=c(1:10)/10)

zz / 1:10



blkl <- c('TCGA-A2-A0T6-01', 'TCGA-E2-A15I-01')


cosmic_sig <- getCosmicSignatures()

matchSignatures(blca_analysis$Results$signatures)#, reference = cosmic_sig[c(1, 3, 5, 13)])

# Retrieve signatures (results)
blca.sig <- blca.analysis$Results$signatures

# Retrieve exposures (results)
blca.exp <- blca.analysis$Results$exposures

# Plot signature 1 (standard barplot, you can pass extra args such as ylim)
msigPlot(blca.sig, signature = 1, ylim = c(0, 0.10))

# Plot exposures (ggplot2 object, you can customize as any other ggplot2 object)
msigPlot(blca.exp, main = "BLCA samples") + 
  scale_fill_manual(values = c("#1f78b4", "#cab2d6", "#ff7f00", "#a6cee3"))


Visualize match
# signature 1 is similar to COSMIC ; 
# signatures 2 and 3 are similar to COSMIC
# Here, we should probably extract only 2 mutational signatures
hm1 <- msig1$plot + ggtitle("Match to COSMIC signs.")
hm2 <- msig2$plot + ggtitle("Match to known BLCA signs.")

# Show
grid.arrange(hm1, hm2, ncol = 2)
