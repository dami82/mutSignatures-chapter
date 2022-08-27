# Load libs
#
library(mutSignatures)
library(dplyr)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(kableExtra)

# Set working dir
#
setwd('~/mutSignatures_chapter/03_method/')

# import mutation data from a list of VCF files
#
vcf_ext <- 'vcf.gz'
de_novo_fld <- 'VCF_de_novo/'
my_vcfs <- grep(pattern = vcf_ext, dir(de_novo_fld), value = TRUE)
my_vcfs <- paste0(de_novo_fld, my_vcfs)
head(my_vcfs)

# Read a list of VCF files via importVCFfiles()
#
x <- mutSignatures::importVCFfiles(vcfFiles = my_vcfs)
dim(x)

# Pre-processing: filter, attach nucleotide context, count mutation types
#

# Optional: custom filtering (eg. minimum QUALITY of 10)
# Optional: adjust SAMPLE identifiers
x <- x %>%
  dplyr::filter(QUAL >= 10) %>%
  mutate(SAMPLEID = sub('__simulated$', '', .$SAMPLEID)) %>%
  mutate(SAMPLEID = gsub('[[:punct:]]', '-', toupper(.$SAMPLEID)))
dim(x)

# Filter non SNV
x <- filterSNV(dataSet = x,  seq_colNames = c("REF", "ALT"))
dim(x)

# Visualize head
head(x) %>% kable() %>% kable_styling(bootstrap_options = "striped")

# Assure that the sequence Names match
# If seqNames do NOT match, an error will be returned
head(BSgenome::seqnames(BSgenome.Hsapiens.UCSC.hg19))

# Attach 3-nt context
x <- attachContext(mutData = x, 
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

# Now the data.frame includes a column informing about the mutType and one about the sample identifier.
# We shall count mutation types per sample.

# Count mutation types per sample
#
blca_counts <- countMutTypes(mutTable = x, 
                             mutType_colName = 'mutType', 
                             sample_colName = 'SAMPLEID')

# my_counts is a mutationCounts object (class in mutSignatures)
# and is one of the inputs of the analysis
blca_counts

# Define parameters for the non-negative matrix factorization procedure.
# you should parallelize if possible
blca_params <- setMutClusterParams( 
    num_processesToExtract = 5,       # num signatures to extract
    num_totIterations = 15,           # bootstrapping: usually 500-1000
    num_parallelCores = 1, 
    seed = 12345)            # number of cores to use 

# De Novo Extraction of Mutational Signatures (may take a while)
#
#set.seed(12345)
blca_analysis <- decipherMutationalProcesses(input = blca_counts,
                                             params = blca_params)

# Results
de_novo_sigs <- blca_analysis$Results$signatures[c(3, 4, 5, 2, 1)]
de_novo_sigs <- setSignatureNames(de_novo_sigs, paste0('blca_', LETTERS[1:5]))


# Visualize
msigPlot(de_novo_sigs, signature = 1)
msigPlot(de_novo_sigs, signature = 2)
msigPlot(de_novo_sigs, signature = 3)
msigPlot(de_novo_sigs, signature = 4)
msigPlot(de_novo_sigs, signature = 5)

# Cosmic
cosmix <- getCosmicSignatures()
blca_cosmix <- cosmix[c(1,2,5,13)]
blca_new_vs_cosmix <- matchSignatures(de_novo_sigs, blca_cosmix)
blca_new_vs_cosmix$plot

# Exposures
my_pal <- c('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854')
msigPlot(blca_analysis$Results$exposures) +
  scale_fill_manual(values = my_pal)



##
##
##


# Deconvolute using known mutational signatures
#

vcf_ext <- 'vcf.gz'
deconv_fld <- 'VCF_deconv/'
my_vcfs <- grep(pattern = vcf_ext, dir(deconv_fld), value = TRUE)
my_vcfs <- paste0(deconv_fld, my_vcfs)

# Read a list of VCF files via importVCFfiles()
#
blca_deco_counts <- mutSignatures::importVCFfiles(vcfFiles = my_vcfs) %>%
  dplyr::filter(QUAL >= 10) %>%
  mutate(SAMPLEID = sub('__simulated$', '', .$SAMPLEID)) %>%
  mutate(SAMPLEID = gsub('[[:punct:]]', '-', toupper(.$SAMPLEID))) %>%
  
  filterSNV(seq_colNames = c("REF", "ALT")) %>%
  attachContext(BSGenomeDb = BSgenome.Hsapiens.UCSC.hg19, 
                chr_colName = 'CHROM', start_colName = 'POS', 
                end_colName = 'POS', nucl_contextN = 3, 
                context_colName = 'context') %>%
  
  removeMismatchMut(refMut_colName = "REF",       # column name for ref base
                    context_colName = "context",  # column name for context
                    refMut_format = "N")    %>%
  
  attachMutType(ref_colName = "REF",              # column name for ref base
                var_colName = "ALT",              # column name for mut base
                context_colName = "context") %>%
  
  countMutTypes(mutType_colName = 'mutType', sample_colName = 'SAMPLEID')

# Deconvolute using resolveMutSignatures()
blca_deco_expos <- resolveMutSignatures(mutCountData = blca_deco_counts,
                                        signFreqData = de_novo_sigs,
                                        byFreq = TRUE)

#  Build plot and Visualize
my_pal <- c('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854')
msigPlot(blca_deco_expos$results$count.result) +
  scale_fill_manual(values = my_pal)

# Export as data.frame
deco_expo <- mutSignatures::as.data.frame(blca_deco_expos$results$count.result, 
                                          transpose = TRUE)
head(deco_expo)

# Export freq.version as data.frame
freq_expo <- mutSignatures::as.data.frame(blca_deco_expos$results$freq.result, 
                                          transpose = TRUE)

# Sort by abs num of mutations
sample_ord <- apply(deco_expo, 1, sum) %>% sort(decreasing = TRUE) %>% names()

tmp <- rbind(
  reshape2:::melt.matrix(as.matrix(deco_expo), varnames = c('ID', 'Sign')) %>%
    mutate(panel = 'Abs Mut Counts)'),
  reshape2:::melt.matrix(as.matrix(freq_expo), varnames = c('ID', 'Sign')) %>%
    mutate(panel = 'Frequency')) %>%
  
  mutate(ID = factor(.$ID, levels = sample_ord))
  
ggplot(tmp, aes(x = ID, y = value, fill = Sign)) +
  geom_bar(stat = 'identity', position = position_stack()) +
  theme_bw() +
  facet_grid(rows = vars(panel), scales = 'free', space = 'fixed') +
  ylab('Mutations') + xlab('') +
  scale_fill_manual('Signatures', values = my_pal) +
  theme(axis.text.x = element_text(angle = 60, size = 7, hjust = 1, vjust = 1), 
        legend.position = 'top')
  
  
  
  

  

