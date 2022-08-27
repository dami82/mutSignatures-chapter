
# working firectory
#
setwd('~/mutSignatures_chapter/02_build_input_files/')

# load libs
#
library(dplyr)

# custom f(x)
#

# write_tcga_vcf
write_tcga_vcf <- function(x, meta_li, sample_colname = 'sample', 
                           min_n = 20, dest_dir = '.') {
  
  # work with flds
  cur_dir <- getwd()
  if (!dir.exists(dest_dir))
    try(dir.create(dest_dir))
  
  tryCatch({setwd(dest_dir)}, error = function(e) { 
    stop('Could not change to destination directory')})
  
  x <- x[!is.na(x[, sample_colname]), ]
  z <- split(x = x, f = x[, sample_colname])
  
  smpls <- names(z)
  y <- list()
  
  for (si in smpls) {
    
    flnm <- paste0(gsub('[[:punct:]]', '_', si), '__simulated.vcf')
    zi <- z[[si]]
    zi <- zi[, -which(colnames(zi) == sample_colname)]
    zi <- zi[as.numeric(order(zi$POS)), ]
    zi <- zi[order(zi$POS), ]
    zi <- zi[order(zi$CHROM), ]
    rownames(zi) <- NULL
    
    try(unlink(x = flnm), silent = TRUE)
    
    if (nrow(zi) >=  min_n) {
      
      y[[length(y) + 1]] <- flnm
      
      # write meta lines 
      for (li in meta_li) {
        if (grepl('^##', x = li))
          write(x = li, file = flnm, append = TRUE)
      }
      
      # write header
      hea <- paste0('#', paste(colnames(zi), collapse = '\t'))
      write(x = hea, file = flnm, append = TRUE)
      
      # write data
      write.table(x = zi, file = flnm, append = TRUE, quote = FALSE, 
                  sep = '\t', row.names = FALSE, col.names = FALSE)
    }
  }
  
  y <- do.call(c, y)
  
  setwd(cur_dir)
  return(y)
} 

# fix_indels
fix_indels <- function(x){
  
  # in (ref): REF
  # out (alt): ALT
  y <- list()
  keep <- nchar(x$REF) == 1
  my_cols <- colnames(x)
  
  if (sum(keep) > 0)
    y[[length(y) + 1]] <- x[keep, ]
  
  if (sum(!keep) > 0) {
    xx <- x[!keep, ]
    rownames(xx) <- NULL
    
    for (i in 1:nrow(xx)) {
      tmp <- NULL
      
      if(xx$ALT[i] == '-'){
        
        # explode the deletion
        tmp <- data.frame(
          CHROM = xx$CHROM[i], 
          POS = seq(xx$POS[i], (xx$POS[i] + nchar(xx$REF[i]) - 1), by = 1), 
          ID = xx$ID[i], 
          REF = strsplit(xx$REF[i], split = '')[[1]], 
          ALT = xx$ALT[i], 
          QUAL = xx$QUAL[i], 
          FILTER = xx$FILTER[i], 
          INFO = xx$INFO[i], 
          stringsAsFactors = FALSE)
        
        for (ci in my_cols[!my_cols %in% colnames(tmp)]) {  
          tmp[, ci] <- xx[i, ci] 
        }
        
        y[[length(y) + 1]] <- tmp
        
      } else if (nchar(xx$ALT[i]) == nchar(xx$REF[i])) {
        
        # explode the mutation  
        tmp <- data.frame(
          CHROM = xx$CHROM[i], 
          POS = seq(xx$POS[i], (xx$POS[i] + nchar(xx$REF[i]) - 1), by = 1), 
          ID = xx$ID[i], 
          REF = strsplit(xx$REF[i], split = '')[[1]], 
          ALT = strsplit(xx$ALT[i], split = '')[[1]], 
          QUAL = xx$QUAL[i], 
          FILTER = xx$FILTER[i], 
          INFO = xx$INFO[i], 
          stringsAsFactors = FALSE)
        
        for (ci in my_cols[!my_cols %in% colnames(tmp)]) {  
          tmp[, ci] <- xx[i, ci] 
        }
        
        y[[length(y) + 1]] <- tmp
      }
    }
  }
  y  <- do.call(rbind, y)   
  rownames(y) <- NULL
  y <- y[order(as.numeric(y$POS)), ]
  y <- y[order(y$CHROM), ]
  
  return(y)  
  
}

##
## Perform Analysis
##

# Input MAF
#
url_1 <- 'https://gdac.broadinstitute.org/runs/analyses__2016_01_28/'
url_2 <- 'reports/cancer/BRCA-TP/MutSigNozzleReport2CV/'
my_maf <- 'BRCA-TP.final_analysis_set.maf'

# Download MAF file
#
download.file(url = paste0(url_1, url_2, my_maf), destfile = url_3)

# Read in and re-format data
#
x <- read.delim(my_maf, header = 1, as.is = TRUE)
x <- x[, c('Chromosome', 'Start_Position', 'Reference_Allele', 
           'Tumor_Seq_Allele2', 'patient', 'i_dbNSFP_CADD_phred')]

x <- data.frame(
  CHROM = paste0('chr', x$Chromosome),
  POS = as.numeric(x$Start_Position),
  ID = '.',
  REF = x$Reference_Allele, 
  ALT = x$Tumor_Seq_Allele2, 
  QUAL = as.numeric(
    sub('\\|[\\.0-9]+$', '', sub('^[\\.0-9]+\\|', '', x$i_dbNSFP_CADD_phred))), 
  FILTER = '.', 
  INFO= '.',
  patient = substr(x$patient, 1, 15),
  stringsAsFactors = FALSE)

x <- x %>%
  dplyr::filter(.$CHROM %in% paste0('chr', c(1:22, 'X', 'Y'))) %>%
  dplyr::filter(.$ALT != .$REF, !is.na(.$POS)) %>%
  dplyr::mutate(QUAL = ifelse(is.na(.$QUAL), '', .$QUAL)) %>%
  fix_indels()

# Declare meta-info lines (simulate a real VCF file)
#
meta_li <- list(
  '##fileformat=VCFv4.2',
  '##FILTER=<ID=PASS,Description="All filters passed">',
  '##samtoolsVersion=1.2+htslib-1.2.1',
  '##samtoolsCommand=samtools mpileup -f ../index/hg19.fa -v sample.bam', 
  '##reference=file://../index/hg19.fa', 
  '##contig=<ID=chr10,length=135534747>',
  '##contig=<ID=chr11,length=135006516>',
  '##contig=<ID=chr12,length=133851895>',
  '##contig=<ID=chr13,length=115169878>',
  '##contig=<ID=chr14,length=107349540>',
  '##contig=<ID=chr15,length=102531392>',
  '##contig=<ID=chr16,length=90354753>',
  '##contig=<ID=chr17,length=81195210>',
  '##contig=<ID=chr18,length=78077248>',
  '##contig=<ID=chr19,length=59128983>',
  '##contig=<ID=chr1,length=249250621>',
  '##contig=<ID=chr2,length=243199373>',
  '##contig=<ID=chr3,length=198022430>',
  '##contig=<ID=chr4,length=191154276>',
  '##contig=<ID=chr5,length=180915260>',
  '##contig=<ID=chr6,length=171115067>',
  '##contig=<ID=chr7,length=159138663>',
  '##contig=<ID=chr8,length=146364022>',
  '##contig=<ID=chr9,length=141213431>',
  '##contig=<ID=chrX,length=155270560>',
  '##contig=<ID=chrY,length=59373566>',
  '##ALT=<ID=X,Description="Represents allele(s) other than observed.">',
  '##bcftools_callVersion=1.2-22-gdf1d0a0+htslib-1.2.1-74-g15a3be5',
  '##bcftools_callCommand=call -cv -Ov sample.vcf.gz')

# Build simulated VCF files 
#
my_vcf_files <- write_tcga_vcf(x = x, meta_li = meta_li, 
                               sample_colname = 'patient', min_n = 30, 
                               dest_dir = 'VCF')

# Clean
#
unlink(x = my_maf)


