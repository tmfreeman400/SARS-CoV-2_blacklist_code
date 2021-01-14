#Main R script for covid manuscript, incorporating elements of get_aligner_comparison_stats_v2.R, load_solo_DB_AFs.R and generate_combined_blacklist.R
#In that order.

#Load all necessary packages

#install.packages('ggplot2')
#BiocManager::install("DEGreport")
#install.packages('devtools')
#install_github("ggobi/ggally")
#install.packages('seqinr')
#install.packages("cowplot")
#install_github("Financial-Times/ftplottools")
#install.packages('ggpubr')
#install.packages('DiagrammeR')
#install.packages('tidyverse')
#install.packages('plyr')
##Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1)
#install.packages("V8")
#install.packages('DiagrammeRsvg')
##system("sudo apt-get install -y librsvg2-dev")
#install.packages('rsvg')
#install.packages("VennDiagram")

library('devtools')
library(GGally)
library(ggplot2)
library(DEGreport)
library('seqinr')
library("cowplot")
library(ftplottools)
require(ggpubr)
library('DiagrammeR')
library('plyr')
library('tidyverse')
library("V8")
library('DiagrammeRsvg')
library('rsvg')
library(VennDiagram) 

#Set ggplot theme for plots
ggplot2::theme_set(ggplot2::theme_bw())

#Change this code to match the directories where your files are saved:
main_project_directory <- '/media/tim/DATA/MovedDownloads/coronavirus_experiments/' #This is the main directory for inputs and outputs
vcf_directory <- '/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/' #This is the directory containing all VCF files,
#which should be saved within nested folders in the form "variant_caller/basecaller/SHEF-yourpatientID.othertext.vcf"
blacklist_directory <- '/media/tim/DATA/MovedDownloads/coronavirus_experiments/blacklists/' #Where we will save blacklist files CONSIDER HAVING THIS SAME AS MAIN
soloDBs_directory <- '/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/soloDBs/' #Where soloDB files are saved for all patients

#Also provide paths to the following required files:
fasta_path = '/media/tim/DATA/MovedDownloads/coronavirus_experiments/nCoV-2019.reference.fasta' #Path to SARS-CoV-2 reference sequence fasta file
bull_et_al_blacklist_path <- '/media/tim/DATA/MovedDownloads/coronavirus_experiments/blacklists/ONT_blacklist_regions_deveson.covid.vcf'
demaio_et_al_blacklist_path <- '/media/tim/DATA/MovedDownloads/coronavirus_experiments/blacklists/problematic_sites_sarsCov2_updated_DEC_goldman.vcf'
artic_primers_path <- "/media/tim/DATA/MovedDownloads/coronavirus_experiments/artic_primers_V3.bed"
qPCR_primers_path <- "/media/tim/DATA/MovedDownloads/coronavirus_experiments/qPCR_primers.bed"
sheffield_sgRNA_path <- "/media/tim/DATA/MovedDownloads/coronavirus_experiments/all_novel.csv"

#All edits to this file should be above this line. No edits should be required after this point for reproducing code.

##SECTION 1: Load in all summarised VCF files and generate aggregate data frame showing which patient,
#basecaller and variant caller each variant is called with, and at what allelic fraction

#Set working directory to coronavirus folder subdirectory containing VCF files
setwd(vcf_directory)

#Predefined function used to glue all values together within a given row of a data frame
rowglue <- function(df){vectorresult <- rep(NA, nrow(df)); for(i in 1:nrow(df)){vectorresult[i] <- paste0(df[i,], collapse = "")}; return(vectorresult)}

#Predefined function to capitalise first letter of string
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

#For each variant-caller/basecaller combination, read in the variants from the corresponding ALL_summary3 file
vcallers <- c("medaka", "nanopolish")
basecallers <- c("hac", "hac3", "rle", "flipflop")

all_vdata <- data.frame(POS=NA, REF=NA, ALT=NA, PATIENT=NA, BASECALLER=NA, VCALLER=NA)
all_vdatab <- data.frame(POS=NA, REF=NA, ALT=NA, AF=NA, ZMETRIC=NA, PATIENT=NA, BASECALLER=NA, VCALLER=NA)
#Create example first row with NA values to be removed later
for(vcaller in vcallers){
  for(basecaller in basecallers){
    vdata <- read.table(paste0("./", vcaller, "/", basecaller, "/ALL_summary3"), sep = "\t")
    vdata <- vdata[which(vdata$V10>0),] #Filter out variants that lack supporting reads
    vdata_2 <- vdata[,c(2, 4, 5, 12)]
    if(basecaller == "hac"){basecaller <- "hac4"}
    vdata_3 <- cbind(vdata_2, basecaller, vcaller)
    vdata_3b <- cbind(vdata[,c(2, 4, 5, 10, 11, 12)], basecaller, vcaller)
    names(vdata_3) <- names(all_vdata)
    names(vdata_3b) <- names(all_vdatab)
    all_vdata <- rbind(all_vdata, vdata_3)
    all_vdatab <- rbind(all_vdatab, vdata_3b)
  }
}
all_vdata2 <- all_vdata[-1,]
all_vdata2b <- all_vdatab[-1,]

##END OF SECTION 1

##SECTION 2: Plot pairwise combinations of basecallers (correlograms) for each variant caller,
#showing the allelic fraction values of all called variants in the first set of plots
#and the population frequency values (within the 884 patient cohort) in the second set of plots

# Create data frame with only AF data for each patient/variant combination, 1st col patient, 2nd col variant, remaining columns basecallers
#Get list of all unique patient/variant combinations

for(temp_vcaller in c("nanopolish", "medaka")){
  
AF_set_all <- all_vdata2b[which(all_vdata2b$VCALLER == temp_vcaller),]
AF_set_glued <- rowglue(AF_set_all[,c(2,1,3,6)])
unique_pat_vars <- unique(AF_set_glued)
gg_AF_data <- data.frame(PATVAR=unique_pat_vars, hac3=NA, hac4=NA, rle=NA, flipflop=NA)
#For each unique patient/variant combination get AF for hac3, hac4, rle and flipflop
for(varnum in 1:length(unique_pat_vars)){
  temp_var <- unique_pat_vars[varnum]
  try(hac3index <- which(AF_set_glued==temp_var & AF_set_all$BASECALLER=="hac3"))
  if(length(hac3index)==1){gg_AF_data$hac3[varnum] <- AF_set_all$AF[hac3index]}else{gg_AF_data$hac3[varnum] <- 0}
  try(hac4index <- which(AF_set_glued==temp_var & AF_set_all$BASECALLER=="hac4"))
  if(length(hac4index)==1){gg_AF_data$hac4[varnum] <- AF_set_all$AF[hac4index]}else{gg_AF_data$hac4[varnum] <- 0}
  try(rleindex <- which(AF_set_glued==temp_var & AF_set_all$BASECALLER=="rle"))
  if(length(rleindex)==1){gg_AF_data$rle[varnum] <- AF_set_all$AF[rleindex]}else{gg_AF_data$rle[varnum] <- 0}
  try(flipflopindex <- which(AF_set_glued==temp_var & AF_set_all$BASECALLER=="flipflop"))
  if(length(flipflopindex)==1){gg_AF_data$flipflop[varnum] <- AF_set_all$AF[flipflopindex]}else{gg_AF_data$flipflop[varnum] <- 0}
}
gg_AF_data2 <- gg_AF_data[,2:5]

# Check correlations (as scatterplots), distribution and print Pearson correlation coefficients 
lowerfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping)+
    geom_point(alpha = 0.1, size=0.01)+
    geom_line(aes(x = ..y..), color="red") +
    #geom_cor(method = "pearson", hjust=0.5, vjust=0, size =2.5, col='red') +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1))
}
lowerfun_bigger_points <- function(data,mapping){
  ggplot(data = data, mapping = mapping)+
    geom_point(alpha = 0.2, size=2)+
    geom_abline(intercept = 0, color="red") +
    #geom_cor(method = "pearson", hjust=0.5, vjust=0, size =2.5, col='red') +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1))
}

#Create GGALLY plot using ggpairs()

p_ggpairs <- ggpairs(data = gg_AF_data2, title=paste0("Allele frequency correlogram between basecallers at variants called in ", temp_vcaller),
                     lower = list(continuous = wrap(lowerfun))) + theme(panel.spacing = unit(0.75, "lines"))

pdf(file=paste0(main_project_directory, "corr_bcaller_AF_", temp_vcaller,".pdf"), width = 7, height = 6)
print(p_ggpairs)
dev.off()

png(file=paste0(main_project_directory, "corr_bcaller_AF_", temp_vcaller,".png"), width = 525, height = 450)
print(p_ggpairs)
dev.off()

#Also need to plot correlogram showing the proportions of patients in which each variant was called for each basecaller
#For each patient in gg_AF_data replace AF with 1 if called and 0 if not called, then strip patient info and count number of each variant
gg_AF_data_ALT <- gg_AF_data
gg_AF_data_ALT[,2:5] <- as.numeric(gg_AF_data_ALT[,2:5]>0)
for(i in 1:nrow(gg_AF_data_ALT)){
  gg_AF_data_ALT$PATVAR[i] <- unlist(strsplit(gg_AF_data_ALT$PATVAR[i], split="SHEF"))[1]
}
#Sum together rows with identical values in first column and divide by total number of patients (884)
dup_ALT <- which(duplicated(gg_AF_data_ALT$PATVAR))
nondup_ALT <- which(!duplicated(gg_AF_data_ALT$PATVAR))
for(i in dup_ALT){
  values <- gg_AF_data_ALT[i,]
  first_rep_index <- which(gg_AF_data_ALT$PATVAR == as.character(values[1]))[1]
  gg_AF_data_ALT[first_rep_index, 2:5] <- gg_AF_data_ALT[first_rep_index, 2:5] + values[2:5]
  gg_AF_data_ALT[i, 2:5] <- c(0,0,0,0)
}
gg_PATPROP_data <- gg_AF_data_ALT[nondup_ALT,]
gg_PATPROP_data[2:5] <- gg_PATPROP_data[2:5]/884
#Can now plot corresponding correlogram
gg_PATPROP_data2 <- gg_PATPROP_data[,2:5]
p_ggpairs_ALT <- ggpairs(data = gg_PATPROP_data2, title=paste0(CapStr(temp_vcaller), " variant frequency correlogram between basecallers in 884 patients"),
                              lower = list(continuous = wrap(lowerfun_bigger_points)),
                              diag = list(continuous = "barDiag")) + theme(panel.spacing = unit(0.75, "lines"))

pdf(file=paste0(main_project_directory, "corr_bcaller_PF_", temp_vcaller,".pdf"), width = 7, height = 6)
print(p_ggpairs_ALT)
dev.off()

png(file=paste0(main_project_directory, "corr_bcaller_PF_", temp_vcaller,".png"), width = 525, height = 450)
print(p_ggpairs_ALT)
dev.off()

}

##END OF SECTION 2

##SECTION 3: Filter data from section 1 to generate data frame of only variants that are called across all basecallers and variant callers

#Generate data frame of all variants' AF values across all 8 methods, and filter for only variants that are conserved
all_vdata2b_glued <- rowglue(all_vdata2b[,c(2,1,3,6)])
unique_pat_vars <- unique(all_vdata2b_glued)
gg_AF_data_all <- data.frame(PATVAR=unique_pat_vars, nano_hac3=NA, nano_hac4=NA, nano_rle=NA, nano_flipflop=NA,
                             med_hac3=NA, med_hac4=NA, med_rle=NA, med_flipflop=NA)
#For each unique patient/variant combination get AF for hac3, hac4, rle and flipflop
for(varnum in 1:length(unique_pat_vars)){
  temp_var <- unique_pat_vars[varnum]
  try(nano_hac3index <- which(all_vdata2b_glued==temp_var & all_vdata2b$BASECALLER=="hac3" & all_vdata2b$VCALLER=="nanopolish"))
  if(length(nano_hac3index)==1){gg_AF_data_all$nano_hac3[varnum] <- all_vdata2b$AF[nano_hac3index]}else{gg_AF_data_all$nano_hac3[varnum] <- 0}
  try(nano_hac4index <- which(all_vdata2b_glued==temp_var & all_vdata2b$BASECALLER=="hac4" & all_vdata2b$VCALLER=="nanopolish"))
  if(length(nano_hac4index)==1){gg_AF_data_all$nano_hac4[varnum] <- all_vdata2b$AF[nano_hac4index]}else{gg_AF_data_all$nano_hac4[varnum] <- 0}
  try(nano_rleindex <- which(all_vdata2b_glued==temp_var & all_vdata2b$BASECALLER=="rle" & all_vdata2b$VCALLER=="nanopolish"))
  if(length(nano_rleindex)==1){gg_AF_data_all$nano_rle[varnum] <- all_vdata2b$AF[nano_rleindex]}else{gg_AF_data_all$nano_rle[varnum] <- 0}
  try(nano_flipflopindex <- which(all_vdata2b_glued==temp_var & all_vdata2b$BASECALLER=="flipflop" & all_vdata2b$VCALLER=="nanopolish"))
  if(length(nano_flipflopindex)==1){gg_AF_data_all$nano_flipflop[varnum] <- all_vdata2b$AF[nano_flipflopindex]}else{gg_AF_data_all$nano_flipflop[varnum] <- 0}
  
  try(med_hac3index <- which(all_vdata2b_glued==temp_var & all_vdata2b$BASECALLER=="hac3" & all_vdata2b$VCALLER=="medaka"))
  if(length(med_hac3index)==1){gg_AF_data_all$med_hac3[varnum] <- all_vdata2b$AF[med_hac3index]}else{gg_AF_data_all$med_hac3[varnum] <- 0}
  try(med_hac4index <- which(all_vdata2b_glued==temp_var & all_vdata2b$BASECALLER=="hac4" & all_vdata2b$VCALLER=="medaka"))
  if(length(med_hac4index)==1){gg_AF_data_all$med_hac4[varnum] <- all_vdata2b$AF[med_hac4index]}else{gg_AF_data_all$med_hac4[varnum] <- 0}
  try(med_rleindex <- which(all_vdata2b_glued==temp_var & all_vdata2b$BASECALLER=="rle" & all_vdata2b$VCALLER=="medaka"))
  if(length(med_rleindex)==1){gg_AF_data_all$med_rle[varnum] <- all_vdata2b$AF[med_rleindex]}else{gg_AF_data_all$med_rle[varnum] <- 0}
  try(med_flipflopindex <- which(all_vdata2b_glued==temp_var & all_vdata2b$BASECALLER=="flipflop" & all_vdata2b$VCALLER=="medaka"))
  if(length(med_flipflopindex)==1){gg_AF_data_all$med_flipflop[varnum] <- all_vdata2b$AF[med_flipflopindex]}else{gg_AF_data_all$med_flipflop[varnum] <- 0}
  
}
gg_AF_data_all2 <- gg_AF_data_all[,2:9]
#Conserved variants do not have 0 in any columns, only keep rows without 0 values
keeprow <- c()
for(i in 1:nrow(gg_AF_data_all)){
  if(!any(gg_AF_data_all2[i,]==0)){keeprow <- c(keeprow, i)}
}
gg_AF_data_all_recurrent <- gg_AF_data_all[keeprow,] # Slightly fewer than half of the rows are kept
#Rank rows by mean AF of columns 2-9 (i.e. lowest AF recurrent variants at the top)
gg_AF_data_recurrent_ordered <- gg_AF_data_all_recurrent[order(rowMeans(gg_AF_data_all_recurrent[,2:9])),]


##END OF SECTION 3

##SECTION 4: Identify variants that only exist in a single basecaller or variant caller for generating blacklists.

gg_AF_single_bc <- gg_AF_data_all
gg_AF_single_bc$singleBC <- FALSE
for(i in 1:nrow(gg_AF_single_bc)){
  if( all((gg_AF_single_bc[i,1:9]==0) == c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)) |
      all((gg_AF_single_bc[i,1:9]==0) == c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)) |
      all((gg_AF_single_bc[i,1:9]==0) == c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE)) |
      all((gg_AF_single_bc[i,1:9]==0) == c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE))
  ){
    gg_AF_single_bc$singleBC[i] <- TRUE
  }
}
gg_AF_single_bc_filt <- gg_AF_single_bc[which(gg_AF_single_bc$singleBC),]
#Get vector of all unique positions covered by these variants for use in fisher exact test (extract positions from PATVAR column)
single_bc_pos_vec <- unique(as.numeric(gsub("[^0-9.-]", "", (unlist(strsplit(gg_AF_single_bc_filt$PATVAR, split='SHEF'))[seq(from=1, by=2, length.out = length(gg_AF_single_bc_filt$PATVAR))]))))

#Also get variants that exist in all basecallers and both variant callers
gg_AF_all_bc <- gg_AF_data_all
gg_AF_all_bc$allBC <- FALSE
for(i in 1:nrow(gg_AF_all_bc)){
  if( all((gg_AF_all_bc[i,1:9]==0) == c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  ){
    gg_AF_all_bc$allBC[i] <- TRUE
  }
}
gg_AF_all_bc_filt <- gg_AF_all_bc[which(gg_AF_all_bc$allBC),]
#Get vector of all unique positions covered by these variants for use in fisher exact test
all_bc_pos_vec <- unique(as.numeric(gsub("[^0-9.-]", "", (unlist(strsplit(gg_AF_all_bc_filt$PATVAR, split='SHEF'))[seq(from=1, by=2, length.out = length(gg_AF_all_bc_filt$PATVAR))]))))

##END OF SECTION 4

##SECTION 5: Order called variants by mean AF, then get list of all variants called at low allelic fraction values.

#Identify recurrent called variants that are low AF
table(unlist(strsplit(head(gg_AF_data_recurrent_ordered$PATVAR, 100), "SHEF"))[seq(from=1, by=2, length.out = 100)])

#Extract mean non-zero gg_AF_data_all AF values for each PATVAR
gg_AF_data_extra <- gg_AF_data_all 
gg_AF_data_extra$meannonzeroAF <- NA
gg_AF_data_extra$PAT <- gg_AF_data_extra$VAR <-NA
for(i in 1:nrow(gg_AF_data_extra)){
  gg_AF_data_extra$meannonzeroAF[i] <- mean(as.numeric(gg_AF_data_extra[i,2:9][which(gg_AF_data_extra[i,2:9] > 0)]))
  #Add PAT column to allow identification of patients with unusually low AFs across multiple called variants
  gg_AF_data_extra$PAT[i] <- paste0("SHEF-", unlist(strsplit(gg_AF_data_extra$PATVAR[i], split="-"))[2])
  gg_AF_data_extra$VAR[i] <- unlist(strsplit(gg_AF_data_extra$PATVAR[i], split="SHEF"))[1]
}
#If many variants are present at low AFs, this could indicate co-infection. Check distribution of AFs for each patient
unique_pat_AFdistr <- data.frame(PAT=unique(gg_AF_data_extra$PAT))
for(i in 1:nrow(unique_pat_AFdistr)){
  unique_pat_AFdistr$min[i] <- min(gg_AF_data_extra$meannonzeroAF[which(gg_AF_data_extra$PAT == unique_pat_AFdistr$PAT[i])], na.rm = TRUE)
  unique_pat_AFdistr$mean[i] <- mean(gg_AF_data_extra$meannonzeroAF[which(gg_AF_data_extra$PAT == unique_pat_AFdistr$PAT[i])], na.rm = TRUE)
  unique_pat_AFdistr$q10[i] <- quantile(gg_AF_data_extra$meannonzeroAF[which(gg_AF_data_extra$PAT == unique_pat_AFdistr$PAT[i])], na.rm = TRUE, probs = 0.1)
  unique_pat_AFdistr$var_count[i] <- length(gg_AF_data_extra$meannonzeroAF[which(gg_AF_data_extra$PAT == unique_pat_AFdistr$PAT[i])])
}

#Order gg_AF_data_extra by meannonzeroAF column
gg_AF_data_extra <- gg_AF_data_extra[order(gg_AF_data_extra$meannonzeroAF),]

#Which patients have multiple recurrently-called low AF variants?
gg_AF_data_extra_recurrent <- gg_AF_data_all_recurrent 
gg_AF_data_extra_recurrent$meannonzeroAF <- NA
gg_AF_data_extra_recurrent$PAT <- gg_AF_data_extra_recurrent$VAR <-NA
for(i in 1:nrow(gg_AF_data_extra_recurrent)){
  gg_AF_data_extra_recurrent$meannonzeroAF[i] <- mean(as.numeric(gg_AF_data_extra_recurrent[i,2:9][which(gg_AF_data_extra_recurrent[i,2:9] > 0)]))
  #Add PAT column to allow identification of patients with unusually low AFs across multiple called variants
  gg_AF_data_extra_recurrent$PAT[i] <- paste0("SHEF-", unlist(strsplit(gg_AF_data_extra_recurrent$PATVAR[i], split="-"))[2])
  gg_AF_data_extra_recurrent$VAR[i] <- unlist(strsplit(gg_AF_data_extra_recurrent$PATVAR[i], split="SHEF"))[1]
}
#If many variants are present at low AFs, this could indicate co-infection. Check distribution of AFs for each patient
unique_pat_AFdistr_recurrent <- data.frame(PAT=unique(gg_AF_data_extra_recurrent$PAT))
for(i in 1:nrow(unique_pat_AFdistr_recurrent)){
  unique_pat_AFdistr_recurrent$min[i] <- min(gg_AF_data_extra_recurrent$meannonzeroAF[which(gg_AF_data_extra_recurrent$PAT == unique_pat_AFdistr_recurrent$PAT[i])], na.rm = TRUE)
  unique_pat_AFdistr_recurrent$mean[i] <- mean(gg_AF_data_extra_recurrent$meannonzeroAF[which(gg_AF_data_extra_recurrent$PAT == unique_pat_AFdistr_recurrent$PAT[i])], na.rm = TRUE)
  unique_pat_AFdistr_recurrent$q10[i] <- quantile(gg_AF_data_extra_recurrent$meannonzeroAF[which(gg_AF_data_extra_recurrent$PAT == unique_pat_AFdistr_recurrent$PAT[i])], na.rm = TRUE, probs = 0.1)
  unique_pat_AFdistr_recurrent$var_count[i] <- length(gg_AF_data_extra_recurrent$meannonzeroAF[which(gg_AF_data_extra_recurrent$PAT == unique_pat_AFdistr_recurrent$PAT[i])])
}
#View list of patients with lowest mean AF recurrent variants
unique_pat_AFdistr_recurrent[head(match(sort(unique_pat_AFdistr_recurrent$mean), (unique_pat_AFdistr_recurrent$mean)), 10),]
#Look at recurrent variants in the 6 patients with the lowest mean allelic fraction recurrent variants.
gg_AF_data_extra_recurrent[which(gg_AF_data_extra_recurrent$PAT %in% unique_pat_AFdistr_recurrent[
  head(match(sort(unique_pat_AFdistr_recurrent$mean), (unique_pat_AFdistr_recurrent$mean)), 10),1][1:6]),]
#SHEF-C7C6B has 4 recurrent variants, all of which are within the range of 0.49-0.59 mean AF - potentially suggesting a coinfection,
#if this is not a result of sequencing error

##END OF SECTION 5

##SECTION 6: Get summary data for each unique variant, including the total number of patients each variant was called in,
#the minimum, mean and max numbers of patients across different basecallers and variant callers, and the numbers of
#cases where basecallers and variant callers gave different results for whether that unique variant was called

#Get list of unique variants found in all combinations
unique_variants <- unique(all_vdata2[,1:3])

##For each of these variants, we want to annotate the following info:
#What is the total number of unique patients with this variant across all basecallers and variant callers?
unique_variants$TOTPATUNIQUE <- NA
#How many patients on average (per basecaller/variant-caller combination) have this variant?
unique_variants$NPATMEAN <- NA
#What is the minimum number of patients with this variant in any given basecaller/variant-caller combination?
unique_variants$NPATMIN <- NA
#What is the maximum number of patients with this variant for a single given basecaller/variant-caller combination?
unique_variants$NPATMAX <- NA
#In how many patient/variant caller combinations are there discrepancies between the basecallers' results?
unique_variants$BASECALLERDIS <- NA
#In how many patient/basecaller combinations are there discrepancies between the variant callers' results?
unique_variants$VCALLERDIS <- NA

#Precompute glued columns 1-3 of all_vdata2 for faster calculation in loop later
glueddata <- rowglue(all_vdata2[,1:3])

for(varnum in 1:nrow(unique_variants)){
  #Get the NPAT values from the all_vdata2 data frame
  given_variant <- unique_variants[varnum,]
  #For a given unique variant, identify which basecallers/variant-callers the variant callers the variant is found in
  given_vdata <- all_vdata2[which(glueddata==paste0(given_variant[1:3], collapse = "")),]
  unique_variants$TOTPATUNIQUE[varnum] <- length(unique(given_vdata[,4]))
  NPATvalues <- table(rowglue(given_vdata[-4]))
  if(length(NPATvalues)!=8){NPATvalues <- c(NPATvalues, 0)}
  unique_variants$NPATMEAN[varnum] <- sum(NPATvalues)/8
  unique_variants$NPATMIN[varnum] <- min(NPATvalues)
  unique_variants$NPATMAX[varnum] <- max(NPATvalues)
  
  given_table <- table(factor(given_vdata[,4], levels = unique(all_vdata[,4])),
                       factor(given_vdata[,5], levels = c(BASECALLER=c("hac4", "hac3", "rle", "flipflop"))),
                       factor(given_vdata[,6], levels = c(VCALLER=c("medaka", "nanopolish"))))
  
  #If there is more than 1 patient with a given variant:
  agreeing_patients <- which(rowSums(given_table[,,])%in%c(0,8))
  discrepancy_patients <- which((rowSums(given_table[,,])%%8)!=0)
  discrepancy_table <- given_table[discrepancy_patients,,]
  
  #Get number of patients where basecallers disagree on variant
  if(length(discrepancy_patients)==0){between_basecaller_discrepancy <- 0}
  if(length(discrepancy_patients)==1){between_basecaller_discrepancy <- sum((rowSums(rbind(discrepancy_table[,1],discrepancy_table[,2])) %% 4) != 0)}
  if(length(discrepancy_patients)>1){between_basecaller_discrepancy <- sum((rowSums(rbind(discrepancy_table[,,1],discrepancy_table[,,2])) %% 4) != 0)}
  
  #Get number of patients where vcallers disagree on variant
  if(length(discrepancy_patients)==0){between_vcaller_discrepancy <- 0}
  if(length(discrepancy_patients)==1){between_vcaller_discrepancy <- sum((rowSums(rbind(discrepancy_table[1,],discrepancy_table[2,],discrepancy_table[3,],discrepancy_table[4,])) %% 2) != 0)}
  if(length(discrepancy_patients)>1){between_vcaller_discrepancy <- sum((rowSums(rbind(discrepancy_table[,1,],discrepancy_table[,2,],discrepancy_table[,3,],discrepancy_table[,4,])) %% 2) != 0)}
  
  unique_variants$BASECALLERDIS[varnum] <- between_basecaller_discrepancy
  unique_variants$VCALLERDIS[varnum] <- between_vcaller_discrepancy
  
}

#Sort unique_variants by numbers of discrepancies
variants_ordered <- unique_variants[order(-unique_variants$BASECALLERDIS, -unique_variants$VCALLERDIS),]

#Label whether a variant is an SNV or INDEL:
variants_ordered$TYPE <- "INDEL"
variants_ordered$TYPE[which((variants_ordered$REF%in%c("A","C","G","T")) & (variants_ordered$ALT%in%c("A","C","G","T")))] <- "SNV"

#How, if at all, should BASECALLERDIS and VCALLERDIS be scaled?
#Scaling by the total number of patients with that variant reveals discrepancies per patient
variants_ordered$BASECALLERDISPERPAT <- variants_ordered$BASECALLERDIS/variants_ordered$TOTPATUNIQUE
variants_ordered$VCALLERDISPERPAT <- variants_ordered$VCALLERDIS/variants_ordered$TOTPATUNIQUE

#Create similar data frame of numbers of discrepancies, but ignoring the exact REF/ALT changes
variants_ordered_stripped <- variants_ordered[,c(1,4,8,9)]
#For rows with the same POS value, sum other columns together, since discrepancies exist at multiple different
#variants at the same position, in different patients
dup_index2 <- which(duplicated(variants_ordered_stripped$POS))
nondup_index2 <- which(!duplicated(variants_ordered_stripped$POS))
for(i in dup_index2){
  values <- variants_ordered_stripped[i,]
  first_rep_index <- which(variants_ordered_stripped$POS == as.numeric(values[1]))[1]
  variants_ordered_stripped[first_rep_index, 2:4] <- variants_ordered_stripped[first_rep_index, 2:4] + values[2:4]
  variants_ordered_stripped[i, 2:4] <- c(0,0,0)
}
variants_ordered_stripped <- variants_ordered_stripped[nondup_index2,]
#Calculate per patient statistics
variants_ordered_stripped$BASECALLERDISPERPAT <- variants_ordered_stripped$BASECALLERDIS/variants_ordered_stripped$TOTPATUNIQUE
variants_ordered_stripped$VCALLERDISPERPAT <- variants_ordered_stripped$VCALLERDIS/variants_ordered_stripped$TOTPATUNIQUE

##END OF SECTION 6

##SECTION 7: Create BED file of variants that are only called in one variant caller (medaka, nanopolish),
#also create BED file of variants that are not called in all basecallers (hac3, hac4, rle, flipflop)

all_rows <- variants_ordered[,1:4]
all_rows$only_called_with <- NA
names(all_rows)[4] <- 'FREQ_COUNT'
#Annotate VC in which variant is called and write to 2 separate BED files
for(i in 1:nrow(all_rows)){
  relevant_rows <- gg_AF_data_extra[which(gg_AF_data_extra$VAR == paste0(all_rows[i,c(2,1,3)], collapse = '')),]
  if(sum(relevant_rows[,2:5])==0){VCanno <- 'medaka'}else{
    if(sum(relevant_rows[,6:9])==0){VCanno <- 'nanopolish'}else{
      VCanno <- 'bothVC'}}
  all_rows$only_called_with[i] <- VCanno
}
#Majority of variants are only called in medaka: Plot Venn diagram showing this

venn.diagram(list(Medaka = which(all_rows$only_called_with %in% c("medaka", "bothVC")),
                  Nanopolish = which(all_rows$only_called_with %in% c("nanopolish", "bothVC"))), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, paste0(main_project_directory, "vcaller_overlap_venn_diagram.tiff"),
             main="Unique called variants", imagetype = 'tiff', margin=0.15, ext.text=FALSE,
             main.pos=c(0.55, 0.95), height = 2000, width = 2000, cat.default.pos="outer")

venn.diagram(list(Medaka = which(all_rows$only_called_with %in% c("medaka", "bothVC")),
                  Nanopolish = which(all_rows$only_called_with %in% c("nanopolish", "bothVC"))), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, paste0(main_project_directory, "vcaller_overlap_venn_diagram.png"),
             main="Unique called variants", imagetype = 'png', margin=0.15, ext.text=FALSE,
             main.pos=c(0.55, 0.95), height = 2000, width = 2000, cat.default.pos="outer")

all_rows <- all_rows[which(all_rows$only_called_with!="bothVC"),]
all_rows <- all_rows[order(all_rows[,1]),]
variants_only_called_in_medaka <- all_rows[which(all_rows$only_called_with == 'medaka'),]
variants_only_called_in_nanopolish <- all_rows[which(all_rows$only_called_with == 'nanopolish'),]
write.table(x=variants_only_called_in_medaka[,1:4], file = paste0(blacklist_directory,'only_called_in_medaka.bed'), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(x=variants_only_called_in_nanopolish[,1:4], file = paste0(blacklist_directory,'only_called_in_nanopolish.bed'), quote = FALSE, row.names = FALSE, sep = "\t")


#Create BED file of variants that are never called across all 4 basecallers, with both variant callers
variants_never_called_from_all_4_BC_bothVC <- variants_ordered[which(variants_ordered$BASECALLERDISPERPAT==2),1:4]
variants_never_called_from_all_4_BC_bothVC$ONLY_CALLED_WITH <- NA
variants_never_called_from_all_4_BC_bothVC$NOT_CALLED_WITH <- NA
names(variants_never_called_from_all_4_BC_bothVC)[4] <- 'FREQ_COUNT'
#Annotate which basecallers variant is/isn't called with
for(i in 1:nrow(variants_never_called_from_all_4_BC_bothVC)){
  relevant_rows <- gg_AF_data_extra[which(gg_AF_data_extra$VAR == paste0(variants_never_called_from_all_4_BC_bothVC[i,c(2,1,3)], collapse = '')),]
  BCanno1 <- paste0(colnames(relevant_rows[2:9])[which(colSums(relevant_rows[,2:9])>0)], collapse = ';')
  variants_never_called_from_all_4_BC_bothVC$ONLY_CALLED_WITH[i] <- BCanno1
  BCanno2 <- paste0(colnames(relevant_rows[2:9])[which(colSums(relevant_rows[,2:9])==0)], collapse = ';')
  variants_never_called_from_all_4_BC_bothVC$NOT_CALLED_WITH[i] <- BCanno2
}
variants_never_called_from_all_4_BC_bothVC <- variants_never_called_from_all_4_BC_bothVC[order(variants_never_called_from_all_4_BC_bothVC[,1]),]
write.table(x=variants_never_called_from_all_4_BC_bothVC, file = paste0(blacklist_directory,'variants_never_called_with_all_BC_bothVC.bed'), quote = FALSE, row.names = FALSE, sep = "\t")


#Create BED file of variants that are never called across all 4 basecallers, with one variant caller
variants_never_called_from_all_4_BC_oneVC <- variants_ordered[which(variants_ordered$BASECALLERDISPERPAT==1),1:4]
variants_never_called_from_all_4_BC_oneVC$inconsistently_called_with <- NA
variants_never_called_from_all_4_BC_oneVC$ONLY_CALLED_WITH <- NA
variants_never_called_from_all_4_BC_oneVC$NEVER_CALLED_WITH <- NA
names(variants_never_called_from_all_4_BC_oneVC)[4] <- 'FREQ_COUNT'
#Annotate VC in which variant is inconsistently called across BCs
for(i in 1:nrow(variants_never_called_from_all_4_BC_oneVC)){
  relevant_rows <- gg_AF_data_extra[which(gg_AF_data_extra$VAR == paste0(variants_never_called_from_all_4_BC_oneVC[i,c(2,1,3)], collapse = '')),]
  if(all(relevant_rows[,2:5]>0) | all(relevant_rows[,2:5]==0)){VCanno <- 'medaka'}else{
    if(all(relevant_rows[,6:9]>0) | all(relevant_rows[,6:9]==0)){VCanno <- 'nanopolish'}else{
      errorCondition(paste0('Check this value: ', i))}}
  variants_never_called_from_all_4_BC_oneVC$inconsistently_called_with[i] <- VCanno
  #Annotate which basecallers variant is called in/ not called in, in the VC which has inconsistent calling
  if(VCanno=='nanopolish'){
    BCanno1 <- BCanno2 <- c()
    if(sum(relevant_rows[,2])>0){BCanno1 <- c(BCanno1, 'HAC3')}else{BCanno2 <- c(BCanno2, 'HAC3')}
    if(sum(relevant_rows[,3])>0){BCanno1 <- c(BCanno1, 'HAC4')}else{BCanno2 <- c(BCanno2, 'HAC4')}
    if(sum(relevant_rows[,4])>0){BCanno1 <- c(BCanno1, 'RLE')}else{BCanno2 <- c(BCanno2, 'RLE')}
    if(sum(relevant_rows[,5])>0){BCanno1 <- c(BCanno1, 'FLIPFLOP')}else{BCanno2 <- c(BCanno2, 'FLIPFLOP')}
  }
  if(VCanno=='medaka'){
    BCanno1 <- BCanno2 <- c()
    if(sum(relevant_rows[,6])>0){BCanno1 <- c(BCanno1, 'HAC3')}else{BCanno2 <- c(BCanno2, 'HAC3')}
    if(sum(relevant_rows[,7])>0){BCanno1 <- c(BCanno1, 'HAC4')}else{BCanno2 <- c(BCanno2, 'HAC4')}
    if(sum(relevant_rows[,8])>0){BCanno1 <- c(BCanno1, 'RLE')}else{BCanno2 <- c(BCanno2, 'RLE')}
    if(sum(relevant_rows[,9])>0){BCanno1 <- c(BCanno1, 'FLIPFLOP')}else{BCanno2 <- c(BCanno2, 'FLIPFLOP')}
  }  
  BCanno1p <- paste0(BCanno1, collapse = ';')
  variants_never_called_from_all_4_BC_oneVC$ONLY_CALLED_WITH[i] <- BCanno1p
  BCanno2p <- paste0(BCanno2, collapse = ';')
  variants_never_called_from_all_4_BC_oneVC$NEVER_CALLED_WITH[i] <- BCanno2p
}
#Remove variants that are never called and write to 2 separate, ordered BED files
variants_never_called_from_all_4_BC_oneVC <- variants_never_called_from_all_4_BC_oneVC[which(variants_never_called_from_all_4_BC_oneVC$ONLY_CALLED_WITH!=''),]
variants_never_called_from_all_4_BC_oneVC <- variants_never_called_from_all_4_BC_oneVC[order(variants_never_called_from_all_4_BC_oneVC[,1]),]
variants_inconsistently_called_in_medaka <- variants_never_called_from_all_4_BC_oneVC[which(variants_never_called_from_all_4_BC_oneVC$inconsistently_called_with == 'medaka'),]
variants_inconsistently_called_in_nanopolish <- variants_never_called_from_all_4_BC_oneVC[which(variants_never_called_from_all_4_BC_oneVC$inconsistently_called_with == 'nanopolish'),]
write.table(x=variants_inconsistently_called_in_medaka[,c(1:4,6:7)], file = paste0(blacklist_directory,'inconsistently_called_in_medaka.bed'), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(x=variants_inconsistently_called_in_nanopolish[,c(1:4,6:7)], file = paste0(blacklist_directory,'inconsistently_called_in_nanopolish.bed'), quote = FALSE, row.names = FALSE, sep = "\t")

##END OF SECTION 7

##SECTION 8: Plot workflow diagrams for the manuscript and save these

#Colours are from the NIHR Brand Guidelines at 60% transparency
# https://cambridgebrc.nihr.ac.uk/wp-content/uploads/2019/04/NIHR-Brand-Guidelines_V1.1-March-2019.pdf

#Part 1: Workflow for generating blacklists
workflow1 <- DiagrammeR::grViz("digraph {
  graph [layout = dot, rankdir = TB, fontsize=10]
  node [shape = rectangle, style=filled]
  rec0 [label = 'Sequence SARS-CoV-2 samples\n from 884 patients', fillcolor = '#FFE5B2']
  rec1 [label = 'SARS-CoV-2 BAM files', fillcolor = '#FFE5B2', shape = folder]
  rec2a [label = 'Get allelic fraction values\n for all samples at all loci', fillcolor = '#F29E95']
  rec2b [label = 'Generate Incremental Database', fillcolor = '#F29E95']
  rec3 [label =  'Generate Monte Carlo \nmodel for 884 samples', fillcolor = '#F29E95']
  rec4 [label = 'Calculate systematic bias for all loci', fillcolor = '#F29E95']
  rec5 [label = 'SARS-CoV-2 VCF files', fillcolor = '#FFE5B2', shape = folder]
  rec6 [label = 'Get blacklist of variants \nwith systematic bias >5%', shape = folder, fillcolor = '#9BCBD0']
  rec7 [label = 'Get blacklist of non-DEL \nvariants with DEL AF >20%', shape = folder, fillcolor = '#9BCBD0']
  rec8 [label = 'Get blacklist of variants \nonly found in one variant caller', shape = folder, fillcolor = '#9BCBD0']
  # edge definitions with the node IDs
  rec0 -> rec1 -> rec2a -> rec2b -> rec4
  rec3 -> rec4 -> rec6
  rec0 -> rec5
  rec5 -> {rec6 rec7 rec8}
  rec2a -> {rec7}
  }")


print(workflow1) %>% 
  export_svg %>% charToRaw %>% rsvg_pdf(paste0(main_project_directory, "workflow1.pdf"))

print(workflow1) %>% 
  export_svg %>% charToRaw %>% rsvg_png(paste0(main_project_directory, "workflow1.png"))

#Part 2: Workflow for analysing co-infection possibility
workflow2 <- DiagrammeR::grViz("digraph {
  graph [layout = dot, rankdir = TB, fontsize=10]
  node [shape = rectangle, style=filled]
  rec0 [label = 'Sequence SARS-CoV-2 samples\n from 884 patients', fillcolor = '#FFE5B2']
  rec1 [label = 'SARS-CoV-2 BAM files', fillcolor = '#FFE5B2', shape = folder]
  rec2a [label = 'Get allelic fraction values\n for all samples at all loci', fillcolor = '#F29E95']
  rec2b [label = 'Generate Incremental Database', fillcolor = '#F29E95']
  rec3 [label =  'Generate Monte Carlo \nmodel for 884 samples', fillcolor = '#F29E95']
  rec4 [label = 'Identify high SD variants', fillcolor = '#A1CAA6']
  rec5 [label = 'SARS-CoV-2 VCF files', fillcolor = '#FFE5B2', shape = folder]
  rec6 [label = 'Identify clustered low AF\n variants with PyClone', fillcolor = '#A1CAA6']
  rec7 [label = 'Check variants are low AF across patient cohort', fillcolor = '#A1CAA6']
  rec8 [label = 'Compare against putative \nco-infection variants in the literature', fillcolor = '#A1CAA6']
  # edge definitions with the node IDs
  rec0 -> rec1 -> rec2a -> rec2b -> rec4
  rec3 -> rec4
  rec0 -> rec5 -> {rec4 rec6}
  {rec4 rec6} -> rec7 -> rec8
  }")

print(workflow2) %>% 
  export_svg %>% charToRaw %>% rsvg_pdf(paste0(main_project_directory, "workflow2.pdf"))

print(workflow2) %>% 
  export_svg %>% charToRaw %>% rsvg_png(paste0(main_project_directory, "workflow2.png"))

##SECTION 9: Read in suspect loci and annotate suspect variants

folder_containing_suspect_loci_folders = main_project_directory #This folder
#contains folders titled A, C, G, T, DEL, INS, corresponding to the suspect loci data sets contained therein.

#Load in suspect loci sets for all 4 basecallers
setwd(folder_containing_suspect_loci_folders)

alleles <- c('A', 'C', 'G', 'T')
bc_list <- c('flipflop', 'hac3', 'hac', 'rle')
i <- 0
j <- 0
recurrent_sus_list <- list()

for(allele in alleles){
  j <- j + 1
  temp_sus_list <- list()
  for(bc in bc_list){
    i <- i + 1
    sus_file <- read.table(paste0(allele, '/', allele, 'chr1_with_z_', bc, '.bed'))
    #Find positions which have a systematic bias greater than 5% in all 4 basecallers and less than 50%
    sus_file_5p <- sus_file[which(sus_file$V4>0.05 & sus_file$V4<0.5),]
    
    #If standard deviation is greater than 0.1 and also larger than the AAF, it is too close
    #to the minimum SD curve and should be filtered out to remove possible false positive
    sus_file_5p <- sus_file_5p[which(sus_file_5p$V5<0.1 | sus_file_5p$V5<sus_file_5p$V4),]
    
    #Check that standard deviation is low at these positions - SD should never be greater than 30%
    #print(sus_file_5p[which(sus_file_5p$V5==max(sus_file_5p$V5)),])
    
    temp_sus_list[[i]] <- sus_file_5p
  }
  #Get list of all suspect loci (3rd column) that are in all 4 basecallers with minimum 5% systematic bias
  recurrent_sus_list[[j]] <- intersect(intersect(temp_sus_list[[1]][,3], temp_sus_list[[2]][,3]),
                                       intersect(temp_sus_list[[3]][,3], temp_sus_list[[4]][,3]))
  i <- 0
}


####
#Get variant list and identify variants that occur at recurrent strong suspect loci
sus_var_temp <- c()
bcmeanthreshold_temp <- c()
cov_ref <- read.fasta(fasta_path, forceDNAtolower = FALSE, seqonly = TRUE)
cov_ref_vec <- c(unlist(strsplit(cov_ref[[1]], split = '')), "END")

for(i in 1:4){
  temp <- recurrent_sus_list[[i]]
  pretemp <- cov_ref_vec[temp]
  posttemp <- rep(c('A', 'C', 'G', 'T')[i], length(temp))
  sus_var_temp <- c(sus_var_temp, paste0(pretemp, temp, posttemp))
  #Get highest mean + 2SD threshold for systematic bias for each variant
  hacsustable <- read.table(paste0('./',c('A', 'C', 'G', 'T')[i],'/',c('A', 'C', 'G', 'T')[i],'chr1_with_z_hac.bed'))
  hac3sustable <- read.table(paste0('./',c('A', 'C', 'G', 'T')[i],'/',c('A', 'C', 'G', 'T')[i],'chr1_with_z_hac3.bed'))
  rlesustable <- read.table(paste0('./',c('A', 'C', 'G', 'T')[i],'/',c('A', 'C', 'G', 'T')[i],'chr1_with_z_rle.bed'))
  flipflopsustable <- read.table(paste0('./',c('A', 'C', 'G', 'T')[i],'/',c('A', 'C', 'G', 'T')[i],'chr1_with_z_flipflop.bed'))
  max_table <- rbind(hacsustable[match(temp, hacsustable$V3),9], hac3sustable[match(temp, hac3sustable$V3),9],
                     rlesustable[match(temp, rlesustable$V3),9], flipflopsustable[match(temp, flipflopsustable$V3),9])
  max_vec <- c(); for(mi in 1:ncol(max_table)){max_vec <- c(max_vec, max(max_table[,mi]))}
  bcmeanthreshold_temp <- c(bcmeanthreshold_temp, max_vec)
}
sus_var <- intersect(gg_AF_data_extra$VAR, sus_var_temp)
sus_var_thresh <- bcmeanthreshold_temp[which(sus_var_temp %in% sus_var)]
#With mean + 2SD value for systematic bias threshold, check that variant AF is greater (if not, variant is likely caused by systematic bias)
plus2z <- sus_var_thresh[match(gg_AF_data_extra[which(gg_AF_data_extra$VAR %in% sus_var),11], sus_var)]
gg_AF_data_SUSVAR <- cbind(gg_AF_data_extra[which(gg_AF_data_extra$VAR %in% sus_var),], plus2z)
gg_AF_data_SUSVAR[which(gg_AF_data_SUSVAR$plus2z > gg_AF_data_SUSVAR$meannonzeroAF),]
#Write BED files for all suspect variants, and separately for confirmed suspect variants
all_sus_var <- unique(gg_AF_data_SUSVAR$VAR)
confirmed_sus_var <- unique(gg_AF_data_SUSVAR$VAR[which(gg_AF_data_SUSVAR$plus2z > gg_AF_data_SUSVAR$meannonzeroAF)])
#Define function that writes variant vector to BED file
write_vtoBED <- function(input_var, BED_filename){
  input_var_BED <- data.frame(POS=0, REF='0', ALT='0')
  for(i in 1:length(input_var)){
    reftemp <- gsub("^(.+?)(?=\\d).*", "\\1", input_var[i], perl = TRUE)
    numtemp <- as.numeric(gsub("[^0-9.-]", "", input_var[i]))
    alttemp <- gsub("[0-9.-]", "", gsub(".*^(.+?)(?=\\d)", "", input_var[i], perl = TRUE))
    input_var_BED <- rbind(input_var_BED, c(numtemp, reftemp, alttemp))
  }
  input_var_BED <- input_var_BED[-1,]
  input_var_BED <- unique(input_var_BED[order(as.numeric(input_var_BED$POS)),])
  write.table(x=input_var_BED, file = BED_filename, quote = FALSE, row.names = FALSE, sep = "\t")
}
write_vtoBED(all_sus_var, paste0(blacklist_directory, 'suspect5_called_variants.bed'))
write_vtoBED(confirmed_sus_var, paste0(blacklist_directory, 'suspect5_called_variants_confirmed2SD.bed'))

setwd(vcf_directory)

##SECTION 10: Annotate variants that have no discrepancies in all patients where they occur, since these are more reliable
#And view the distribution of allelic fraction values for variants are consistently called in all patients where they occur,
#despite low allelic fraction values in some or many of those patients.

#670 unique variants had no discrepancies? 670/3123
variants_no_disc <- variants_ordered[which(variants_ordered$BASECALLERDIS==0 & variants_ordered$VCALLERDIS==0),]

#Could also annotate super-recurrent variants in variants_no_disc, which are recurrent in all patients where they occur
variants_no_disc$VAR <- NA
for(i in 1:nrow(variants_no_disc)){variants_no_disc$VAR[i] <- paste0(variants_no_disc$REF[i], variants_no_disc$POS[i], variants_no_disc$ALT[i])}
gg_AF_data_recurrent_ordered$super_recurrent <- NA
for(i in 1:nrow(gg_AF_data_recurrent_ordered)){
  VARtemp <- unlist(strsplit(gg_AF_data_recurrent_ordered$PATVAR[i], "SHEF"))[1]
  if(VARtemp %in% variants_no_disc$VAR){gg_AF_data_recurrent_ordered$super_recurrent[i] <- TRUE}else{
    gg_AF_data_recurrent_ordered$super_recurrent[i] <- FALSE}
}
#Most of the lowest AF recurrent variants are only low AF in low number of patients (is the AF higher in other patients with these variants?)
gg_AF_data_recurrent_ordered$PAT <- gg_AF_data_recurrent_ordered$VAR <- NA
for(i in 1:nrow(gg_AF_data_recurrent_ordered)){
  tempPATVAR <- unlist(strsplit(gg_AF_data_recurrent_ordered$PATVAR[i], split = 'SHEF'))
  gg_AF_data_recurrent_ordered$PAT[i] <- paste0("SHEF", tempPATVAR[2])
  gg_AF_data_recurrent_ordered$VAR[i] <- tempPATVAR[1]
}
#How much does AF vary for the same variant? Look at lowest AF recurrent variants. Can vary as much as ~0.2-0.8 AF
#Are variants called at lower AF in some patients than others? (possible evidence of recombination, etc.?) Yes, and main distribution looks left-skewed
#Low AF variants do not occur in more than 1 patient in most cases.
table(gg_AF_data_recurrent_ordered$VAR[1:100])

#How many rows indicate variants with AF less than 50% on average? 117
which(rowMeans((gg_AF_data_recurrent_ordered[,2:9])) < 0.5 )

#G11083T, however occurs in 73 of the lowest 100 AF cases, and generally shows a peak around AF of 0.4 across patients, but with wide SD
ggplot((gg_AF_data_recurrent_ordered[which(gg_AF_data_recurrent_ordered$VAR=="G11083T"),]), aes(x=med_flipflop)) + geom_histogram() +
  xlab('Called variant AF with medaka/flipflop') + ggtitle("Histogram of G11083T called variant AF values")

#C21575T (4), A9634T (4), G28857A (2) are the only other variants that are called at less than 0.5 AF in more than one patient
#ggplot((gg_AF_data_recurrent_ordered[which(gg_AF_data_recurrent_ordered$VAR=="C21575T"),]), aes(x=med_flipflop)) + geom_histogram() +
#  xlab('Called variant AF with medaka/flipflop') + ggtitle("Histogram of C21575T called variant AF values")

#ggplot((gg_AF_data_recurrent_ordered[which(gg_AF_data_recurrent_ordered$VAR=="A9634T"),]), aes(x=med_flipflop)) + geom_histogram() +
#  xlab('Called variant AF with medaka/flipflop') + ggtitle("Histogram of A9634T called variant AF values")

#ggplot((gg_AF_data_recurrent_ordered[which(gg_AF_data_recurrent_ordered$VAR=="G28857A"),]), aes(x=med_flipflop)) + geom_histogram() +
#  xlab('Called variant AF with medaka/flipflop') + ggtitle("Histogram of G28857A called variant AF values")

#Check if these are in blacklist regions:
#G11083T Yes, both
#C21575T Yes in De Maio, Not in Bull
#A9634T Not in De Maio, Yes in Bull
#G28857A Not in De Maio or Bull


##END OF SECTION 10

##SECTION 11: Investigation of potential co-infection loci examined

#Which variants appeared at AF values <50% in multiple patients?
variants_lt_50af <- table(gg_AF_data_recurrent_ordered$VAR[1:100])[which(table(gg_AF_data_recurrent_ordered$VAR[1:100])>1)]

#What is position of interest where we are checking if deletion is real/coinfection/systematic bias?
for(POI in as.numeric(gsub("[^0-9.-]", "", names(variants_lt_50af)))){
#Get list of called deletions:
AF_set_indels <- AF_set_all[which(nchar(AF_set_all$REF)>1),]
#Get maximum distance of DEL start position from specific SNV component
max_offset <- max(unique(nchar(AF_set_all$REF)))
#Check for region around POI
# print(AF_set_indels[which(AF_set_indels$POS %in% seq(from=POI-max_offset, to=POI+max_offset)),]) #Uncomment to view result
#No called deletions at this position, so check if DEL at POI is either coinfection or systematic bias:

#Load in DEL IncDBs - these 2 lines can be safely deleted
incdb_DEL_hac <- read.table(paste0(folder_containing_suspect_loci_folders, '/DEL/DELchr1_with_z_hac.bed'))
# print(incdb_DEL_hac[which(incdb_DEL_hac$V3 == POI),]) #Uncomment to view result
#Deletion at position of interest is not due to systematic presence of deletions
#Since V9 column is lower than all of the AF values for called indels in each case
}

#11083 DEL does not have low SD, so deletions are not a result of systematic error, suggesting co-infection could be the source of standard deviation.
#Likewise 21575, 9634, 28857 and 3037  could have co-infection

#(Checked in BASH) All 5 loci checked had very little correlation in deletion percentages between all pairs of loci across all patients
#Clustering of low AF variants does occur, but does not appear to affect these same variants across different patients
#So while co-infection may occur, it seems unlikely to be the source of this specific set of low AF variants occurring in more than 1 patient

#Conclusion 11083, 21575, 9634 and 28857 DELs are not the result of systematic bias towards DEL at that position in the genome,
#but they do not appear to be the result of co-infection either, since DEL fractions across patients would then cluster (be similar) between these loci,
#Since they would all correlate with the same 3rd factors such as time since infection, age, viral load etc.
#Therefore propose that these DELs are not due to co-infection. Could be a result of non-systematic error that is prevalent at those positions.
#Could this be explained by other factors such as recombination?

#In addition, PyClone was used to establish that variants sometimes clustered together at similar intermediate AF values in patients such as
#SHEF-C7C6B , which has 4 variants with AF 59-69%

##END OF SECTION 11

##SECTION 12: Show that INDELs have more discrepancies than SNVs in both basecallers and variant callers
#with plots and KS test

#Are there more discrepancies between basecallers or between variant callers? 
#(More discrepancies between variant callers, but twice as many opportunities for this in each VCALLER)
quantile(variants_ordered$BASECALLERDIS, probs = seq(from=0, to=1, by=0.1))
quantile(variants_ordered$VCALLERDIS, probs = seq(from=0, to=1, by=0.1))

#Plot violin plots for the numbers of discrepancies
ggplot(data = variants_ordered, mapping = aes(x="Aligner", y=BASECALLERDISPERPAT, fill='device')) + geom_violin() + geom_point(position = position_jitterdodge())
#BASECALLERDISPERPAT: 0 means basecallers all agree for all variants, 1 means basecallers have
#discrepancies in one variant caller on average, 2 means basecallers have discrepancies in both variant callers.
#Intermediate values mean basecallers agree in some patients but not others in the same variant caller

#VCALLERDISPERPAT: 0 means variant callers agree for all variants, 1 means variant callers give different results in 1 basecaller on average, etc.
ggplot(data = variants_ordered, mapping = aes(x="Variant_Caller", y=VCALLERDISPERPAT, fill='device')) + geom_violin() + geom_point(position = position_jitterdodge())

#Which loci have most discrepancies? Which basecallers/variant callers give different results most often?

#Which changes typically give discrepancies between basecallers' results? 
table(variants_ordered$TYPE[1:100])
#Only 14/100 variants with most basecaller discrepancies are SNVs, 86/100 are indels
#How does this compare with the positions with the fewest discrepancies?
table(variants_ordered$TYPE[(nrow(variants_ordered)-99):nrow(variants_ordered)])
#95% of variants with fewest basecaller discrepancies are SNVs by contrast! 
#But could this difference just be a result of differences in the total numbers of patients?
ggplot(data=variants_ordered, aes(x=TYPE, y=TOTPATUNIQUE)) + geom_violin()
ggplot(data=variants_ordered, aes(x=TYPE, y=NPATMEAN)) + geom_violin()

#Plot correlation using boxplots and statistical test
ggplot(data=variants_ordered, aes(x=TYPE, y=BASECALLERDISPERPAT)) + geom_violin()
ggplot(data=variants_ordered, aes(x=TYPE, y=VCALLERDISPERPAT)) + geom_violin()

#INDELs appear to genuinely have more discrepancies than SNVs, even when accounting for number of patients:
#Use Kolmogorov-Smirnov test (does not assume normality or equal variances)
ks.test(x=variants_ordered$BASECALLERDISPERPAT[which(variants_ordered$TYPE=="INDEL")],
        y=variants_ordered$BASECALLERDISPERPAT[which(variants_ordered$TYPE=="SNV")],
        alternative = "two.sided") #Post-hoc analysis of distributions shows INDELs have higher BASECALLERDISPERPAT

ks.test(x=variants_ordered$VCALLERDISPERPAT[which(variants_ordered$TYPE=="INDEL")],
        y=variants_ordered$VCALLERDISPERPAT[which(variants_ordered$TYPE=="SNV")],
        alternative = "two.sided") #Post-hoc analysis of distributions shows INDELs have higher VCALLERDISPERPAT


##END OF SECTION 12

##SECTION 13: Test if suspect loci enriched or depleted at positions with discrepancies and show that there is no correlation.

#Are suspect loci enriched or depleted at positions with discrepancies?
#E.g. Suspect loci are caused by systematic biases for an basecaller/vcaller combo. If these systematic biases are only present in some basecallers,
#but not others, then might see enrichment overall.
#However if systematic biases are similar across all combos, then would expect depletion of suspect loci at variants with more discrepancies:
#This would indicate that differences between basecallers/vcallers were not due to systematic effects only present in a single combination

#Read ranked variants into R (these are all <1 SD of systematic bias) and match to the above ordered variants where possible
suspect_table <- data.frame(POS=NA, REF=NA, ALT=NA, NUMPATSUS=NA, BASECALLER=NA, VCALLER=NA)
for(vcaller in vcallers){
  for(basecaller in basecallers){
    rankedchunk <- read.table(file = paste0(vcf_directory, vcaller, "/", basecaller, "/ALL_suspect_variants_rankedtable.vcf"))[,c(2,4,5,8)]
    rankedchunk$BASECALLER <- basecaller
    rankedchunk$VCALLER <- vcaller
    names(rankedchunk) <- names(suspect_table)
    suspect_table <- rbind(suspect_table, rankedchunk)
  }
}
suspect_table <- suspect_table[-1,]
suspect_table$BASECALLER[which(suspect_table$BASECALLER=="hac")] <- "hac4"

#Filter out variants that were suspect - do these have higher or lower levels of discrepancies?
suspect_variants_unique <- unique(rowglue(suspect_table[,c(2,1,3)]))
variants_ordered_in_suspect <- variants_ordered[which(rowglue(variants_ordered[,c(2,1,3)]) %in% suspect_variants_unique),]
variants_ordered_notin_suspect <- variants_ordered[which(! rowglue(variants_ordered[,c(2,1,3)]) %in% suspect_variants_unique),]

variants_ordered$SUSPECT <- FALSE
variants_ordered$SUSPECT[which(rowglue(variants_ordered[,c(2,1,3)]) %in% suspect_variants_unique)] <- TRUE

#Show discrepancy numbers using violin plots:
ggplot(data=variants_ordered, aes(x=SUSPECT, y=BASECALLERDISPERPAT)) + geom_violin()
ggplot(data=variants_ordered, aes(x=SUSPECT, y=VCALLERDISPERPAT)) + geom_violin()
#INDELs appear to genuinely have more discrepancies than SNVs, even when accounting for number of patients:
#Use Kolmogorov-Smirnov test (does not assume normality or equal variances)
ks.test(x=variants_ordered$BASECALLERDISPERPAT[which(variants_ordered$SUSPECT)],
        y=variants_ordered$BASECALLERDISPERPAT[which(variants_ordered$SUSPECT)],
        alternative = "two.sided") #No differences in basecaller discrepancies between suspect and non-suspect variants

ks.test(x=variants_ordered$VCALLERDISPERPAT[which(variants_ordered$SUSPECT)],
        y=variants_ordered$VCALLERDISPERPAT[which(variants_ordered$SUSPECT)],
        alternative = "two.sided") #No differences in variant caller discrepancies between suspect and non-suspect variants


#Can this suggest which basecaller/variant caller combination is likely to give the correct result when a discrepancy occurs?
#Do any basecallers or variant callers consistently perform better in this regard?

##END OF SECTION 13

##SECTION 14: Write code to load soloDB allelic fractions for all patients as a list with 884 elements,
#Use this data to annotate non-DEL variants where more than 20% of reads support a DEL variant at that position
#We will include these variants as one of the blacklists

setwd(soloDBs_directory)

#Get vector of all soloDB file names for given basecaller
for(bc in c('hac3', 'hac', 'rle', 'flipflop')){
  soloDB_filelist <- dir(path = bc)
  soloDB_patIDs <- paste0("SHEF-", unlist(strsplit(soloDB_filelist, split="_SHEF-"))[seq(from=2, by=3, length.out = length(soloDB_filelist))])
  temp_soloDB_list <- list()
  
  #Use for loop to read in each soloDB one at a time and add to list
  for(i in 1:length(soloDB_filelist)){
    temp_soloDB_list[[i]] <- read.table(paste0(bc, '/', soloDB_filelist[i]), col.names = c('A', 'C', 'G', 'T', 'DEL', 'INS'))
    #print(i)
  }
  
  ##Use above list of tables to get deletion AF for each patient/variant combination
  #Load in list of variants: gg_AF_data_extra (get_aligner_comparison_stats_v2.R)
  if(bc == 'flipflop'){gg_AF_data_extra$DEL_AF_FLIPFLOP <- NA}
  if(bc == 'rle'){gg_AF_data_extra$DEL_AF_RLE <- NA}
  if(bc == 'hac3'){gg_AF_data_extra$DEL_AF_HAC3 <- NA}
  if(bc == 'hac'){gg_AF_data_extra$DEL_AF_HAC4 <- NA}
  
  for(i in 1:nrow(gg_AF_data_extra)){
    var_coord <- as.numeric(gsub("[^0-9.-]", "", gg_AF_data_extra$VAR[i]))
    var_pat <- gg_AF_data_extra$PAT[i]
    var_pat_index <- which(soloDB_patIDs == var_pat)
    if(length(var_pat_index)==0){next()}
    
    if(bc == 'flipflop'){gg_AF_data_extra$DEL_AF_FLIPFLOP[i] <- temp_soloDB_list[[var_pat_index]][var_coord, 5]}
    if(bc == 'rle'){gg_AF_data_extra$DEL_AF_RLE[i] <- temp_soloDB_list[[var_pat_index]][var_coord, 5]}
    if(bc == 'hac3'){gg_AF_data_extra$DEL_AF_HAC3[i] <- temp_soloDB_list[[var_pat_index]][var_coord, 5]}
    if(bc == 'hac'){gg_AF_data_extra$DEL_AF_HAC4[i] <- temp_soloDB_list[[var_pat_index]][var_coord, 5]}
    
  }
}

#Get list of variants called with a high DEL AF in one or both variant callers when called (DEL AF > 20%) and write to BED file
all_called_PATVARS <- cbind(all_vdata2[,5:6], rowglue(all_vdata2[,c(2,1,3,4)]))
names(all_called_PATVARS)[3] <- "PATVAR"
#Also load in variant AFs and DEL AFs for all_called_PATVARS from gg_AF_data_extra
all_called_PATVARS$VAF <- NA
all_called_PATVARS$DELAF <- NA
for(i in 1:nrow(all_called_PATVARS)){
  use_col <- NA
  DEL_col <- NA
  if(all_called_PATVARS$VCALLER[i]=='nanopolish'){
    if(all_called_PATVARS$BASECALLER[i]=='hac3'){use_col <- 2; DEL_col <- 13}
    if(all_called_PATVARS$BASECALLER[i]=='hac4'){use_col <- 3; DEL_col <- 14}
    if(all_called_PATVARS$BASECALLER[i]=='rle'){use_col <- 4; DEL_col <- 15}
    if(all_called_PATVARS$BASECALLER[i]=='flipflop'){use_col <- 5; DEL_col <- 16}
  }
  if(all_called_PATVARS$VCALLER[i]=='medaka'){
    if(all_called_PATVARS$BASECALLER[i]=='hac3'){use_col <- 6; DEL_col <- 13}
    if(all_called_PATVARS$BASECALLER[i]=='hac4'){use_col <- 7; DEL_col <- 14}
    if(all_called_PATVARS$BASECALLER[i]=='rle'){use_col <- 8; DEL_col <- 15}
    if(all_called_PATVARS$BASECALLER[i]=='flipflop'){use_col <- 9; DEL_col <- 16}
  }
  tempindex <- which(gg_AF_data_extra$PATVAR == all_called_PATVARS$PATVAR[i])
  all_called_PATVARS$VAF[i] <- gg_AF_data_extra[tempindex, use_col]
  all_called_PATVARS$DELAF[i] <- gg_AF_data_extra[tempindex, DEL_col]
}

PATVAR_flipflop <- all_called_PATVARS$PATVAR[which(all_called_PATVARS$BASECALLER=="flipflop")]
PATVAR_rle <- all_called_PATVARS$PATVAR[which(all_called_PATVARS$BASECALLER=="rle")]
PATVAR_hac3 <- all_called_PATVARS$PATVAR[which(all_called_PATVARS$BASECALLER=="hac3")]
PATVAR_hac4 <- all_called_PATVARS$PATVAR[which(all_called_PATVARS$BASECALLER=="hac4")]

high_del_AF_vars_FLIPFLOP <- gg_AF_data_extra[which((gg_AF_data_extra$PATVAR %in% PATVAR_flipflop) & gg_AF_data_extra$DEL_AF_FLIPFLOP>0.2),"VAR"]
high_del_AF_vars_RLE <- gg_AF_data_extra[which((gg_AF_data_extra$PATVAR %in% PATVAR_rle) & gg_AF_data_extra$DEL_AF_RLE>0.2),"VAR"]
high_del_AF_vars_HAC3 <- gg_AF_data_extra[which((gg_AF_data_extra$PATVAR %in% PATVAR_hac3)>0 & gg_AF_data_extra$DEL_AF_HAC3>0.2),"VAR"]
high_del_AF_vars_HAC4 <- gg_AF_data_extra[which((gg_AF_data_extra$PATVAR %in% PATVAR_hac4) & gg_AF_data_extra$DEL_AF_HAC4>0.2),"VAR"]

#Create BED files with columns POS, REF, ALT, FREQ_COUNT
write_vtoBED_withcount <- function(input_var, BED_filename){
  input_var_BED <- data.frame(POS=0, REF='0', ALT='0')
  for(i in 1:length(input_var)){
    reftemp <- gsub("^(.+?)(?=\\d).*", "\\1", input_var[i], perl = TRUE)
    numtemp <- as.numeric(gsub("[^0-9.-]", "", input_var[i]))
    alttemp <- gsub("[0-9.-]", "", gsub(".*^(.+?)(?=\\d)", "", input_var[i], perl = TRUE))
    input_var_BED <- rbind(input_var_BED, c(numtemp, reftemp, alttemp))
  }
  input_var_BED <- input_var_BED[-1,]
  #Count replicates of each line and name new column FREQ_COUNT
  input_var_BED <- ddply(input_var_BED,.(input_var_BED$POS, input_var_BED$REF, input_var_BED$ALT),nrow)
  colnames(input_var_BED) <- c("POS", "REF", "ALT", "FREQ_COUNT")
  input_var_BED <- unique(input_var_BED[order(as.numeric(input_var_BED$POS)),])
  write.table(x=input_var_BED, file = BED_filename, quote = FALSE, row.names = FALSE, sep = "\t")
}
write_vtoBED_withcount(high_del_AF_vars_FLIPFLOP, paste0(blacklist_directory,'highDEL_called_variants_flipflop.bed'))
write_vtoBED_withcount(high_del_AF_vars_RLE, paste0(blacklist_directory,'highDEL_called_variants_rle.bed'))
write_vtoBED_withcount(high_del_AF_vars_HAC3, paste0(blacklist_directory,'highDEL_called_variants_hac3.bed'))
write_vtoBED_withcount(high_del_AF_vars_HAC4, paste0(blacklist_directory,'highDEL_called_variants_hac4.bed'))

#Also create BED file containing all of the variants that were called with one or more variant caller with a high DEL AF across all basecallers
high_del_AF_vars_ALLBC <- gg_AF_data_extra[which(
  (gg_AF_data_extra$PATVAR %in% PATVAR_flipflop) & gg_AF_data_extra$DEL_AF_FLIPFLOP>0.2 &
    (gg_AF_data_extra$PATVAR %in% PATVAR_rle) & gg_AF_data_extra$DEL_AF_RLE>0.2 &
    (gg_AF_data_extra$PATVAR %in% PATVAR_hac3) & gg_AF_data_extra$DEL_AF_HAC3>0.2 &
    (gg_AF_data_extra$PATVAR %in% PATVAR_hac4) & gg_AF_data_extra$DEL_AF_HAC4>0.2
),"VAR"]

#Note: this is not the same as the intersect of all of the other values, since every individual variant must fulfil all conditions in
#at least 1 individual, while the intersect gives variants that have high DEL with every basecaller, but those conditions can be fulfilled
#across multiple individuals

write_vtoBED_withcount(high_del_AF_vars_ALLBC, paste0(blacklist_directory,'highDEL_called_variants_allbc.bed'))

##END OF SECTION 14

##SECTION 15: Find allelic fraction value of major allele for each genomic position
#in order to compute the mean error rate across the genome for each basecaller

#Find AF of major allele in each row to compute the error rate across the genome
for(bc in c('hac3', 'hac', 'rle', 'flipflop')){
  soloDB_filelist <- dir(path = bc)
  soloDB_patIDs <- paste0("SHEF-", unlist(strsplit(soloDB_filelist, split="_SHEF-"))[seq(from=2, by=3, length.out = length(soloDB_filelist))])
  temp_soloDB_list <- list()
  
  #Use for loop to read in each soloDB one at a time and add to list
  for(i in 1:length(soloDB_filelist)){
    temp_soloDB_list[[i]] <- read.table(paste0(bc, '/', soloDB_filelist[i]), col.names = c('A', 'C', 'G', 'T', 'DEL', 'INS'))
  }
  
  mean_error_rates <- rep(NA, length(temp_soloDB_list))
  
  for(patnum in 1:length(mean_error_rates)){
    AF_mat <- temp_soloDB_list[[patnum]]
    y_indices <- max.col(AF_mat)
    major_AF <- rep(NA, nrow(AF_mat))
    for(x in 1:nrow(AF_mat)){
      major_AF[x] <- AF_mat[x, y_indices[x]] 
    }
    
    mean_error_rates[patnum] <- (1-mean(major_AF[which(major_AF!=0)]))
  }
  
  print(mean(mean_error_rates))
  boxplot(mean_error_rates)
}


##END OF SECTION 15

##SECTION 16: Generate single combined blacklist from all different blacklists calculated so far
#Also annotate which variants appear in other published blacklists (De Maio et al. 2020, Bull et al. 2020)

setwd(blacklist_directory)

#Read in all blacklist groups and merge each group

DEL_allbc <- read.table(paste0(blacklist_directory,'highDEL_called_variants_allbc.bed'), header = TRUE)
DEL_flipflop <- read.table(paste0(blacklist_directory,'highDEL_called_variants_flipflop.bed'), header = TRUE)
DEL_hac3 <- read.table(paste0(blacklist_directory,'highDEL_called_variants_hac3.bed'), header = TRUE)
DEL_hac4 <- read.table(paste0(blacklist_directory,'highDEL_called_variants_hac4.bed'), header = TRUE)
DEL_rle <- read.table(paste0(blacklist_directory,'highDEL_called_variants_rle.bed'), header = TRUE)
#Merge blacklists that have high DEL anywhere and remove duplicates
all_DEL_loci <- rbind(DEL_allbc[,1:3], DEL_flipflop[,1:3], DEL_hac3[,1:3], DEL_hac4[,1:3], DEL_rle[,1:3])
all_DEL_loci <- all_DEL_loci[order(all_DEL_loci[,1]),]
all_DEL_loci <- all_DEL_loci[which(!duplicated(all_DEL_loci)),]
#Annotate which BC each locus has a high DEL in
all_DEL_loci$HIGH_DEL_BASECALLERS <- NA
for(i in 1:nrow(all_DEL_loci)){
  temp_bc_list <- c()
  if(rowglue(all_DEL_loci[i,1:3]) %in% rowglue(DEL_flipflop[,1:3])){temp_bc_list <- c(temp_bc_list, 'FLIPFLOP')}
  if(rowglue(all_DEL_loci[i,1:3]) %in% rowglue(DEL_hac3[,1:3])){temp_bc_list <- c(temp_bc_list, 'HAC3')}
  if(rowglue(all_DEL_loci[i,1:3]) %in% rowglue(DEL_hac4[,1:3])){temp_bc_list <- c(temp_bc_list, 'HAC4')}
  if(rowglue(all_DEL_loci[i,1:3]) %in% rowglue(DEL_rle[,1:3])){temp_bc_list <- c(temp_bc_list, 'RLE')}
  all_DEL_loci$HIGH_DEL_BASECALLERS[i] <- paste0(temp_bc_list, collapse = ';')
}

#Next we want to read in variants with systematic bias and merge with the above
suspect_table <- read.table('suspect5_called_variants.bed', header = TRUE)
suspect_table2 <- read.table('suspect5_called_variants_confirmed2SD.bed', header = TRUE)
suspect_table_notin_DELtable <- suspect_table[which(!rowglue(suspect_table[,1:3]) %in% rowglue(all_DEL_loci[,1:3])), 1:3]
suspect_table_notin_DELtable$HIGH_DEL_BASECALLERS <- NA
all_DEL_loci <- rbind(all_DEL_loci, suspect_table_notin_DELtable)
all_DEL_loci <- all_DEL_loci[order(all_DEL_loci[,1]),]
all_DEL_loci$SYSTEMATIC_BIAS_OVER_5_PERCENT <- FALSE
all_DEL_loci$VARIANT_AF_HIGHER_THAN_SYSTEMATIC_BIAS <- TRUE
for(i in 1:nrow(all_DEL_loci)){
  #Check if variant is in suspect_table and suspect_table2
  if(rowglue(all_DEL_loci[i,1:3]) %in% rowglue(suspect_table[,1:3])){all_DEL_loci$SYSTEMATIC_BIAS_OVER_5_PERCENT[i] <- TRUE}
  if(rowglue(all_DEL_loci[i,1:3]) %in% rowglue(suspect_table2[,1:3])){all_DEL_loci$VARIANT_AF_HIGHER_THAN_SYSTEMATIC_BIAS[i] <- FALSE}
}

#Next we want to check if variants are only called in a single variant caller
only_medaka <- read.table('only_called_in_medaka.bed', header = TRUE)[,1:3]
only_medaka_glued <- rowglue(only_medaka)
only_nanopolish <- read.table('only_called_in_nanopolish.bed', header = TRUE)[,1:3]
only_nanopolish_glued <- rowglue(only_nanopolish)
both_VC <- intersect(only_medaka_glued, only_nanopolish_glued) #Sometimes variants are called in both VC, but separate individuals, so check this
#None of the above variants are called in both VC across multiple individuals

#Add rows to all_DEL_loci once more
only_medaka_notin_DELtable <- only_medaka[which(!only_medaka_glued %in% rowglue(all_DEL_loci[,1:3])), 1:3]
only_nanopolish_notin_DELtable <- only_nanopolish[which(!only_nanopolish_glued %in% rowglue(all_DEL_loci[,1:3])), 1:3]
only_medaka_notin_DELtable$HIGH_DEL_BASECALLERS <- only_nanopolish_notin_DELtable$HIGH_DEL_BASECALLERS <- NA
only_medaka_notin_DELtable$SYSTEMATIC_BIAS_OVER_5_PERCENT <- only_nanopolish_notin_DELtable$SYSTEMATIC_BIAS_OVER_5_PERCENT <- FALSE
only_medaka_notin_DELtable$VARIANT_AF_HIGHER_THAN_SYSTEMATIC_BIAS <- only_nanopolish_notin_DELtable$VARIANT_AF_HIGHER_THAN_SYSTEMATIC_BIAS <- TRUE
all_DEL_loci <- rbind(all_DEL_loci, only_medaka_notin_DELtable, only_nanopolish_notin_DELtable)
all_DEL_loci <- all_DEL_loci[order(all_DEL_loci[,1]),]

#For each row, annotate with variant caller if called with just one variant caller
all_DEL_loci$ONLY_CALLED_WITH_ONE_VARIANT_CALLER <- FALSE
for(i in 1:nrow(all_DEL_loci)){
  #Check if variant is in added tables
  if(rowglue(all_DEL_loci[i,1:3]) %in% only_medaka_glued){all_DEL_loci$ONLY_CALLED_WITH_ONE_VARIANT_CALLER[i] <- 'MEDAKA'}
  if(rowglue(all_DEL_loci[i,1:3]) %in% only_nanopolish_glued){all_DEL_loci$ONLY_CALLED_WITH_ONE_VARIANT_CALLER[i] <- 'NANOPOLISH'}
}

#Next we want to annotate whether basecallers agree on variant calls in medaka (NA for variants only called in nanopolish)
#Likewise for nanopolish
bc_disagree_medaka <- read.table('inconsistently_called_in_medaka.bed', header = TRUE, sep = '\t')
bc_disagree_medaka_glued <- rowglue(bc_disagree_medaka[,1:3])
bc_disagree_nanopolish <- read.table('inconsistently_called_in_nanopolish.bed', header = TRUE, sep = '\t')
bc_disagree_nanopolish_glued <- rowglue(bc_disagree_nanopolish[,1:3])
bc_disagree_bothvc <- read.table('variants_never_called_with_all_BC_bothVC.bed', header = TRUE)[,1:3]
bc_disagree_bothvc_glued <- rowglue(bc_disagree_bothvc)
intersect(bc_disagree_medaka_glued, bc_disagree_nanopolish_glued) #No variants in both bc_disagree sets, so can add both sets simultaneously
bc_disagree_medaka_notin_DELtable <- bc_disagree_medaka[which(!bc_disagree_medaka_glued %in% rowglue(all_DEL_loci[,1:3])), 1:3]
bc_disagree_nanopolish_notin_DELtable <- bc_disagree_nanopolish[which(!bc_disagree_nanopolish_glued %in% rowglue(all_DEL_loci[,1:3])), 1:3]
bc_disagree_bothvc_notin_DELtable <- bc_disagree_bothvc[which(!bc_disagree_bothvc_glued %in% rowglue(all_DEL_loci[,1:3])), 1:3]

bc_disagree_medaka_notin_DELtable$HIGH_DEL_BASECALLERS <- bc_disagree_nanopolish_notin_DELtable$HIGH_DEL_BASECALLERS <- bc_disagree_bothvc_notin_DELtable$HIGH_DEL_BASECALLERS <- NA
bc_disagree_medaka_notin_DELtable$SYSTEMATIC_BIAS_OVER_5_PERCENT <- bc_disagree_nanopolish_notin_DELtable$SYSTEMATIC_BIAS_OVER_5_PERCENT <- bc_disagree_bothvc_notin_DELtable$SYSTEMATIC_BIAS_OVER_5_PERCENT <- FALSE
bc_disagree_medaka_notin_DELtable$VARIANT_AF_HIGHER_THAN_SYSTEMATIC_BIAS <- bc_disagree_nanopolish_notin_DELtable$VARIANT_AF_HIGHER_THAN_SYSTEMATIC_BIAS <- bc_disagree_bothvc_notin_DELtable$VARIANT_AF_HIGHER_THAN_SYSTEMATIC_BIAS <- TRUE
bc_disagree_medaka_notin_DELtable$ONLY_CALLED_WITH_ONE_VARIANT_CALLER <- bc_disagree_nanopolish_notin_DELtable$ONLY_CALLED_WITH_ONE_VARIANT_CALLER <- bc_disagree_bothvc_notin_DELtable$ONLY_CALLED_WITH_ONE_VARIANT_CALLER <- FALSE
all_DEL_loci <- rbind(all_DEL_loci, bc_disagree_medaka_notin_DELtable, bc_disagree_nanopolish_notin_DELtable, bc_disagree_bothvc_notin_DELtable)
all_DEL_loci <- all_DEL_loci[order(all_DEL_loci[,1]),]

#For each row, annotate whether basecallers agree on variant calls in medaka and/or nanopolish (NA for variants only called in other vcaller)
all_DEL_loci$BASECALLERS_AGREE_ON_VARIANT_CALL_MEDAKA <- TRUE
all_DEL_loci$BASECALLERS_AGREE_ON_VARIANT_CALL_NANOPOLISH <- TRUE
all_DEL_loci_glued <- rowglue(all_DEL_loci[,1:3])
for(i in 1:nrow(all_DEL_loci)){
  #Check if variant is in added tables
  if(all_DEL_loci_glued[i] %in% bc_disagree_medaka_glued){all_DEL_loci$BASECALLERS_AGREE_ON_VARIANT_CALL_MEDAKA[i] <- FALSE}
  if(all_DEL_loci_glued[i] %in% bc_disagree_nanopolish_glued){all_DEL_loci$BASECALLERS_AGREE_ON_VARIANT_CALL_NANOPOLISH[i] <- FALSE}
  if(all_DEL_loci_glued[i] %in% bc_disagree_bothvc_glued){all_DEL_loci$BASECALLERS_AGREE_ON_VARIANT_CALL_NANOPOLISH[i] <- all_DEL_loci$BASECALLERS_AGREE_ON_VARIANT_CALL_MEDAKA[i] <- FALSE}
  if(all_DEL_loci$ONLY_CALLED_WITH_ONE_VARIANT_CALLER[i] == 'NANOPOLISH'){all_DEL_loci$BASECALLERS_AGREE_ON_VARIANT_CALL_MEDAKA[i] <- NA}
  if(all_DEL_loci$ONLY_CALLED_WITH_ONE_VARIANT_CALLER[i] == 'MEDAKA'){all_DEL_loci$BASECALLERS_AGREE_ON_VARIANT_CALL_NANOPOLISH[i] <- NA}
}

###Next we want to annotate variants that were called despite having zero supportive reads - 0% AF not supported by soloDBs
###0% AF originates from depth/ files in nanomed_884 directory
##zeroAF <- read.table('zeroAF_called_variants_allbc.bed', header = TRUE)[,1:3]
##zeroAF_glued <- rowglue(zeroAF)
##zeroAF_notin_DELtable <- zeroAF[which(!zeroAF_glued %in% rowglue(all_DEL_loci[,1:3])), 1:3]
##
##zeroAF_notin_DELtable$HIGH_DEL_BASECALLERS <- NA
##zeroAF_notin_DELtable$SYSTEMATIC_BIAS_OVER_5_PERCENT <- FALSE
##zeroAF_notin_DELtable$VARIANT_AF_HIGHER_THAN_SYSTEMATIC_BIAS <- TRUE
##zeroAF_notin_DELtable$ONLY_CALLED_WITH_ONE_VARIANT_CALLER <- FALSE
##zeroAF_notin_DELtable$BASECALLERS_AGREE_ON_VARIANT_CALL_MEDAKA <- TRUE
##zeroAF_notin_DELtable$BASECALLERS_AGREE_ON_VARIANT_CALL_NANOPOLISH <- TRUE
##
##all_DEL_loci <- rbind(all_DEL_loci, zeroAF_notin_DELtable)
##all_DEL_loci <- all_DEL_loci[order(all_DEL_loci[,1]),]
##
###For each row, annotate whether variant has zero supporting reads in primer trimmed BAM files (may have zero supporting reads in original BAM files too)
###Zero supporting reads seem to be caused by supporting reads being at very low coverage combined with mpileup not detecting SNV supporting reads within indels
###Since these are considered INS/DEL rather than ACGT
##all_DEL_loci$ZERO_SUPPORTING_READS <- FALSE
##all_DEL_loci_glued <- rowglue(all_DEL_loci[,1:3])
##for(i in 1:nrow(all_DEL_loci)){
##  #Check if variant is in added tables
##  if(all_DEL_loci_glued[i] %in% zeroAF_glued){all_DEL_loci$ZERO_SUPPORTING_READS[i] <- TRUE}
##}

#Next we want to annotate variants that appear in either the Deveson (Bull et al.) or Goldman (De Maio et al.) blacklists
bull <- read.table(bull_et_al_blacklist_path, header = FALSE, as.is = TRUE)[,c(2,4)]
bull$ALT <- '.' #No specific ALT being masked, so don't add to new list, just annotate columns of pre-existing blacklist variants
#for(i in 1:nrow(bull)){bull$V3[i] <- cov_ref_vec[bull$V2[i]]}
all_DEL_loci_glued2 <- rowglue(all_DEL_loci[,1:2])
all_DEL_loci$IN_DEVESON_BLACKLIST <- FALSE
all_DEL_loci$IN_DEVESON_BLACKLIST[which(all_DEL_loci_glued2 %in% rowglue(bull[,1:2]))] <- TRUE

demaio <- read.table(demaio_et_al_blacklist_path, as.is = TRUE, sep = '\t')[,c(2,4,7)]
demaio_mask <- demaio[which(demaio$V7 == 'mask'),1:2]
demaio_mask$ALT <- '.'
demaio_caution <- demaio[which(demaio$V7 == 'caution'),1:2]
demaio_caution$ALT <- '.'
all_DEL_loci$IN_GOLDMAN_MASK_LIST <- FALSE
all_DEL_loci$IN_GOLDMAN_MASK_LIST[which(all_DEL_loci_glued2 %in% rowglue(demaio_mask[,1:2]))] <- TRUE
all_DEL_loci$IN_GOLDMAN_CAUTION_LIST <- FALSE
all_DEL_loci$IN_GOLDMAN_CAUTION_LIST[which(all_DEL_loci_glued2 %in% rowglue(demaio_caution[,1:2]))] <- TRUE


#Create final column giving recommendation to mask or take caution with variant
#High DEL -> caution (variant may be real, but should be manually checked)
#Systematic bias over 5% -> caution
#Variant AF not greater than systematic bias -> mask
#Only called with one variant caller -> caution
#Basecallers disagree -> caution
#In either mask list -> mask
#In Goldman caution list -> caution
#Annotate caution first, then overwrite with mask if applicable
all_DEL_loci$MASKING_RECOMMENDATION <- NA
for(i in 1:nrow(all_DEL_loci)){
  if(!is.na(all_DEL_loci$HIGH_DEL_BASECALLERS[i])){all_DEL_loci$MASKING_RECOMMENDATION[i] <- 'caution'}
  if(all_DEL_loci$SYSTEMATIC_BIAS_OVER_5_PERCENT[i]){all_DEL_loci$MASKING_RECOMMENDATION[i] <- 'caution'}
  if(all_DEL_loci$ONLY_CALLED_WITH_ONE_VARIANT_CALLER[i] != FALSE){all_DEL_loci$MASKING_RECOMMENDATION[i] <- 'caution'}
  if(!is.na(all_DEL_loci$BASECALLERS_AGREE_ON_VARIANT_CALL_MEDAKA[i])){
    if(all_DEL_loci$BASECALLERS_AGREE_ON_VARIANT_CALL_MEDAKA[i] == FALSE){all_DEL_loci$MASKING_RECOMMENDATION[i] <- 'caution'}
  }
  if(!is.na(all_DEL_loci$BASECALLERS_AGREE_ON_VARIANT_CALL_NANOPOLISH[i])){
    if(all_DEL_loci$BASECALLERS_AGREE_ON_VARIANT_CALL_NANOPOLISH[i] == FALSE){all_DEL_loci$MASKING_RECOMMENDATION[i] <- 'caution'}
  }
  #if(all_DEL_loci$ZERO_SUPPORTING_READS[i] == TRUE){all_DEL_loci$MASKING_RECOMMENDATION[i] <- 'caution'}
  if(all_DEL_loci$IN_GOLDMAN_CAUTION_LIST[i] == TRUE){all_DEL_loci$MASKING_RECOMMENDATION[i] <- 'caution'}
  if(all_DEL_loci$VARIANT_AF_HIGHER_THAN_SYSTEMATIC_BIAS[i] == FALSE){all_DEL_loci$MASKING_RECOMMENDATION[i] <- 'mask'}
  if(all_DEL_loci$IN_GOLDMAN_MASK_LIST[i] == TRUE){all_DEL_loci$MASKING_RECOMMENDATION[i] <- 'mask'}
  if(all_DEL_loci$IN_DEVESON_BLACKLIST[i] == TRUE){all_DEL_loci$MASKING_RECOMMENDATION[i] <- 'mask'}
}
output_table <- all_DEL_loci

#Write header for blacklist table describing each column header in more depth. Write masking recommendation explanation in the main text.

#Write full blacklist table to Freeman_full_blacklist file
write.table(output_table, file = 'Freeman_full_blacklist.xlsx', sep = '\t', quote = FALSE, row.names = FALSE)

#Separate blacklist variants table into medaka and nanopolish specific tables for supplementary tables
medaka_output_table <- output_table[which(output_table$ONLY_CALLED_WITH_ONE_VARIANT_CALLER %in% c(FALSE, "MEDAKA")),c(-9)]
nanopolish_output_table <- output_table[which(output_table$ONLY_CALLED_WITH_ONE_VARIANT_CALLER %in% c(FALSE, "NANOPOLISH")),c(-8)]
write.table(medaka_output_table, file = 'Freeman_medaka_blacklist.xlsx', sep = '\t', quote = FALSE, row.names = FALSE)
write.table(nanopolish_output_table, file = 'Freeman_nanopolish_blacklist.xlsx', sep = '\t', quote = FALSE, row.names = FALSE)


##END OF SECTION 16

##SECTION 17: Visualise distribution of variants in blacklists against each other (remove plots not needed from here)
#Also look at summary stats of how frequently variants occur in different blacklists

#Set colours
navy=rgb(25,62,114,maxColorValue = 255)
coral=rgb(234,93,78,maxColorValue = 255)
orange=rgb(242,147,48,maxColorValue = 255)
yellow=rgb(254,212,122,maxColorValue = 255)
purple=rgb(102,103,173,maxColorValue = 255)
aqua=rgb(46,169,176,maxColorValue = 255)
green=rgb(70,168,108,maxColorValue = 255)
grey=rgb(172,188,195,maxColorValue = 255)

#Plot positions on covid genome
ggplot() +
  geom_histogram(data = output_table, mapping = aes(x = POS), bins = 100, color=purple, fill=purple) +
  ylab("Number of blacklist variants") + xlab("Genomic position")

#Generate similar plots for high DEL variants, variants with systematic bias, and variants that are not consistently called
ggplot() +
  geom_histogram(data = output_table[which(!is.na(output_table$HIGH_DEL_BASECALLERS)),], mapping = aes(x = POS), bins = 100, color=purple, fill=purple) +
  ylab("Number of variants with >20% DEL fraction") + xlab("Genomic position")
ggplot() +
  geom_histogram(data = output_table[which(output_table$SYSTEMATIC_BIAS_OVER_5_PERCENT),], mapping = aes(x = POS), bins = 100, color=purple, fill=purple) +
  ylab("Number of variants with >5% systematic bias") + xlab("Genomic position")
ggplot() +
  geom_histogram(data = output_table[which(output_table$ONLY_CALLED_WITH_ONE_VARIANT_CALLER!=FALSE),], mapping = aes(x = POS), bins = 100, color=purple, fill=purple) +
  ylab("Number of variants only called with one variant caller") + xlab("Genomic position")
ggplot() +
  geom_histogram(data = output_table[which(output_table$BASECALLERS_AGREE_ON_VARIANT_CALL_MEDAKA==FALSE | output_table$BASECALLERS_AGREE_ON_VARIANT_CALL_NANOPOLISH==FALSE),], mapping = aes(x = POS), bins = 100, color=purple, fill=purple) +
  ylab("Number of variants not called with all basecallers") + xlab("Genomic position")

#In conclusion, vast majority of blacklist points are from basecaller discrepancies, 
#which is why the overall histogram is very similar to the density plot for discrepancies.

#Variant caller discrepances, basecaller discrepancies both correlate strongly with each other and with high DEL fraction variants
#Systematic bias has less clear correlation with these. Test this with fisher test:
genome_length <- length(cov_ref_vec)
set1 <- rowglue(suspect_table[, 1:3])
set2 <- rowglue(all_DEL_loci[which(!is.na(all_DEL_loci$HIGH_DEL_BASECALLERS)), 1:3])

a <- length(intersect(set1, set2)) #In suspect loci & high DEL variants
b <- length(setdiff(set1, set2)) #In suspect loci, not high DEL variants
c <- length(setdiff(set2, set1)) #In high DEL variants, not suspect loci
d <- (genome_length - (a+b+c)) #Not in either

f.table <- matrix(c(a,b,c,d), nrow = 2, dimnames = list(c("In high DEL variants", "Not in high DEL variants"),c("In suspect loci", "Not in suspect loci")))
fisher.test(f.table) #Variants with systematic bias are 12X enriched in the set of high DEL variants, with a significant p value of 0.0137

#Therefore systematic bias variants, high DEL AF variants, and variants with basecaller and variant caller discrepancies are all correlated with each other
#Variants that are unreliable have an increased rate of all of these features, and they are not independent of each other.

#Get values for table of numbers of variants with different blacklist features
sum(output_table$SYSTEMATIC_BIAS_OVER_5_PERCENT)
sum(!output_table$VARIANT_AF_HIGHER_THAN_SYSTEMATIC_BIAS)
sum(!is.na(output_table$HIGH_DEL_BASECALLERS))
sum(output_table$ONLY_CALLED_WITH_ONE_VARIANT_CALLER!=FALSE)
length(which((!output_table$BASECALLERS_AGREE_ON_VARIANT_CALL_MEDAKA) | (!output_table$BASECALLERS_AGREE_ON_VARIANT_CALL_NANOPOLISH)))
#
b3 <- gg_AF_data_extra[,c(3:5,7:9,11,12)]
b3_disagree <- b3[which((!((b3$nano_hac4>0 & b3$nano_rle>0 & b3$nano_flipflop>0) | (b3$nano_hac4==0 & b3$nano_rle==0 & b3$nano_flipflop==0))) | 
                          (!((b3$med_hac4>0 & b3$med_rle>0 & b3$med_flipflop>0) | (b3$med_hac4==0 & b3$med_rle==0 & b3$med_flipflop==0))) ),]
length(unique(b3_disagree$VAR))
#
##sum(output_table$ZERO_SUPPORTING_READS)
sum(output_table$IN_DEVESON_BLACKLIST)
sum(output_table$IN_GOLDMAN_MASK_LIST)
sum(output_table$IN_GOLDMAN_CAUTION_LIST)

#Read in homopolymer and primer positions
artic_primers <- read.table(artic_primers_path)
homopolymer_bed <- data.frame(BASE=NA, START=NA, STOP=NA, LENGTH=NA)
#For each position, check if next four positions are the same character
i=1
while(i <(length(cov_ref_vec)-4)){
  base <- cov_ref_vec[i]
  if(all(cov_ref_vec[(i+1):(i+4)] == rep(base, 4))){
    #Position is start of homopolymer, find homopolymer end and generate entry for homopolymer
    j <- (i+4)
    while(cov_ref_vec[i]==cov_ref_vec[j]){j <- (j+1)}
    j <- (j+1) #Return j to last index at which base still matched original base, this is the last position of the homopolymer (1-based)
    homopolymer_bed <- rbind(homopolymer_bed, c(base, i, j, (j-i+1)))
    i <- (j+1) #Update i value to look for next homopolymer
  }else{
    #Position is not homopolymer, update i value by one and check next position
    i <- (i+1)
  }
}
homopolymer_bed <- homopolymer_bed[-1,]
homopolymer_bed[,2] <- as.numeric(homopolymer_bed[,2])
homopolymer_bed[,3] <- as.numeric(homopolymer_bed[,3])
homopolymer_bed[,4] <- as.numeric(homopolymer_bed[,4])

#Read in qPCR primers from Japan and Wuhan qPCR sets
qPCR_bed <- read.table(qPCR_primers_path, skip = 4)
wuhan_qPCR_bed <- qPCR_bed[1:3,]
japan_qPCR_bed <- qPCR_bed[4:6,]
half_width <- 20
ymin1 <- 1.5
ymax1 <- 1.4
ymid <- (ymin1+ymax1)/2

#Get legend for top of plot with 6 types of blacklist + basecaller discrepancies
#Blacklist types are: Systematic bias >5% , DEL AF >20%, One VC only, 0% AF, Bull, De Maio
legend_plot <- ggplot() + ft_theme()+
  #geom_density(data=binned_data, mapping = aes(x=POS, y=BASECALLERDISPERPAT, fill=green), col=NA, stat = 'identity') +
  geom_rect(data = output_table[which(output_table$SYSTEMATIC_BIAS_OVER_5_PERCENT), ], mapping = aes(xmin=POS-half_width, xmax=POS+half_width, ymin=(ymin1-0.1), ymax=(ymax1-0.1), fill="Systematic bias > 5%")) +
  geom_rect(data = output_table[which(!is.na(output_table$HIGH_DEL_BASECALLERS)), ], mapping = aes(xmin=POS-half_width, xmax=POS+half_width, ymin=(ymin1-0.2), ymax=(ymax1-0.2), fill="DEL AF >20%")) +
  geom_rect(data = output_table[which(output_table$ONLY_CALLED_WITH_ONE_VARIANT_CALLER!=FALSE), ], mapping = aes(xmin=POS-half_width, xmax=POS+half_width, ymin=(ymin1-0.3), ymax=(ymax1-0.3), fill="One VC only")) +
  ##geom_rect(data = output_table[which(output_table$ZERO_SUPPORTING_READS), ], mapping = aes(xmin=POS-half_width, xmax=POS+half_width, ymin=(ymin1-0.4), ymax=(ymax1-0.4), fill="0% AF")) +
  geom_rect(data = output_table[which(output_table$IN_DEVESON_BLACKLIST), ], mapping = aes(xmin=POS-half_width, xmax=POS+half_width, ymin=(ymin1-0.4), ymax=(ymax1-0.4), fill="Bull")) +
  geom_rect(data = output_table[which(output_table$IN_GOLDMAN_MASK_LIST | output_table$IN_GOLDMAN_CAUTION_LIST), ], mapping = aes(xmin=POS-half_width, xmax=POS+half_width, ymin=(ymin1-0.5), ymax=(ymax1-0.5), fill="De Maio")) +
  scale_fill_manual(values = c("black", "dark red", aqua, yellow, green), name="Blacklist subset:",
                    labels=c("Bull et al.", "DEL AF >20%", "De Maio et al.", "One VC only", "Systematic bias > 5%")) +
                    guides(fill=guide_legend(nrow=2), byrow=TRUE) +
  geom_rect(data = artic_primers, mapping = aes(xmin=V2, xmax=V3, ymin=(ymin1-0.7), ymax=(ymax1-0.7), colour='Artic_primers')) +
  geom_rect(data = homopolymer_bed, mapping = aes(xmin=START, xmax=STOP, ymin=(ymin1-0.8), ymax=(ymax1-0.8), colour='Homopolymers')) +
  geom_rect(data = japan_qPCR_bed, mapping = aes(xmin=V2, xmax=V3, ymin=(ymin1-0.9), ymax=(ymax1-0.9), colour='Japan qPCR primers')) +
  geom_rect(data = wuhan_qPCR_bed, mapping = aes(xmin=V2, xmax=V3, ymin=(ymin1-1.0), ymax=(ymax1-1.0), colour='Wuhan qPCR primers')) +
  scale_color_manual(values = c(coral, navy, orange, purple), name="Annotation:",
                     labels=c("Artic primers","Homopolymers","Japan qPCR primers","Wuhan qPCR primers"))+
  guides(color = guide_legend(override.aes = list(fill = c(coral, navy, orange, purple)), nrow = 2)) +
  
  annotate(geom="rect",xmin=0,xmax=30000,ymin=ymin1,ymax=ymax1,color="grey80",fill="grey80")+
  annotate(geom="rect",xmin=0,xmax=200,ymin=ymin1,ymax=ymax1,color="grey30",fill="grey30")+
  annotate(geom="rect",xmin=200,xmax=13483,ymin=ymin1,ymax=ymax1,color="grey50",fill="grey50")+
  annotate(geom="text",x=7000,y=ymid,color="white",label="ORF1a",fontface='bold')+
  annotate(geom="rect",xmin=13483,xmax=30000,ymin=ymin1,ymax=ymax1,color="grey60",fill="grey60")+
  annotate(geom="text",x=17000,y=ymid,color="white",label="ORF1b",fontface='bold')+
  annotate(geom="rect",xmin=21532,xmax=30000,ymin=ymin1,ymax=ymax1,color=navy,fill=navy)+
  annotate(geom="text",x=23500,y=ymid,color="white",label="S",fontface='bold')+
  annotate(geom="rect",xmin=25361,xmax=30000,ymin=ymin1,ymax=ymax1,color=coral,fill=coral)+
  annotate(geom="text",x=25750,y=ymid,color="white",label="3a",fontface='bold')+
  annotate(geom="rect",xmin=26213,xmax=30000,ymin=ymin1,ymax=ymax1,color=green,fill=green)+
  annotate(geom="rect",xmin=26449,xmax=30000,ymin=ymin1,ymax=ymax1,color=purple,fill=purple)+
  annotate(geom="text",x=26700,y=ymid,color="white",label="M",fontface='bold')+
  annotate(geom="rect",xmin=27017,xmax=30000,ymin=ymin1,ymax=ymax1,color=yellow,fill=yellow)+
  annotate(geom="rect",xmin=27364,xmax=30000,ymin=ymin1,ymax=ymax1,color=orange,fill=orange)+
  annotate(geom="rect",xmin=27864,xmax=30000,ymin=ymin1,ymax=ymax1,color=grey,fill=grey)+
  annotate(geom="rect",xmin=28235,xmax=30000,ymin=ymin1,ymax=ymax1,color=aqua,fill=aqua)+
  annotate(geom="text",x=28800,y=ymid,color="white",label="N",fontface='bold')+
  annotate(geom="rect",xmin=29510,xmax=30000,ymin=ymin1,ymax=ymax1,color="grey80",fill="grey80")+
  annotate(geom="rect",xmin=0,xmax=30000,ymin=ymin1,ymax=ymax1,fill="white",alpha=0.3) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),#legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) + xlab(label = NULL) + coord_cartesian(xlim=c(1000, 29000))

disc_plot_legend <- get_legend(legend_plot)
disc_plot_nolegend <- legend_plot + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                                          axis.title.x=element_blank(),
                                          axis.title.y=element_blank(),legend.position="none",
                                          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                          panel.grid.minor=element_blank(), panel.spacing.y = unit(0,"null")) + xlab(label = NULL) +
  coord_cartesian(xlim=c(0, 30000), ylim = c(0.4,1.5), expand = FALSE)

#Get basecaller discrepancies density for bottom of plot
output_table_bdisc <- output_table[which((!output_table$BASECALLERS_AGREE_ON_VARIANT_CALL_MEDAKA) | (!output_table$BASECALLERS_AGREE_ON_VARIANT_CALL_NANOPOLISH)), ]

#Look at the distribution of all variants and recurrent variants, compared with the distribution of variants with basecaller discrepancies
all_unique_variant_positions <- sort(variants_ordered$POS)
ggplot() + 
  geom_histogram(data=output_table_bdisc, mapping = aes(x=POS), binwidth = 200)

ggplot() + 
  geom_histogram(data=variants_ordered, mapping = aes(x=POS), binwidth = 200)

non_bdisc_variants <- variants_ordered[which(!(rowglue(variants_ordered[,1:3]) %in% rowglue(output_table_bdisc[,1:3]))),]
ggplot() + 
  geom_histogram(data=non_bdisc_variants, mapping = aes(x=POS), binwidth = 200)

#The peaks in blacklisted variants also occur at the same positions for non-blacklisted variants, so these are positions where
#there are larger numbers of unique variants detected in the SARS CoV-2 genome in general. However, the peaks are far less pronounced for
#non-blacklisted variants, which show a much higher proportion of variants at non-peak positions.

nrow(output_table_bdisc)/(nrow(non_bdisc_variants) + nrow(output_table_bdisc))
#70.01% of all unique variants have at least one basecaller discrepancy, so unsurprisingly, the majority of variants found were
#in the blacklist. We recommend "caution" with these rather than outright masking, since many of these are likely real variants.

b_disc_plot <- ggplot() + #ft_theme()+
  #geom_density(data=binned_data, mapping = aes(x=POS, y=BASECALLERDISPERPAT, fill="black"), col=NA, stat = 'identity') +
  geom_histogram(data=output_table_bdisc, mapping = aes(x=POS), binwidth = 200)+
  theme(legend.position = "none", axis.line=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank()) +
  xlab("") + ylab("Number of variants with\n basecaller discrepancies") + coord_cartesian(xlim=c(1000, 29000))

#ggarrange(legend_plot, b_disc_plot, nrow=2, heights=c(1, 0.4), align = 'hv', labels = c('A', 'B'), label.y = c(1, 1.1))

b_disc_plot2 <- ggplot() + #ft_theme()+
  geom_histogram(data=non_bdisc_variants, mapping = aes(x=POS), binwidth = 200)+
  theme(legend.position = "none", axis.line=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  xlab("Locus") + ylab("Number of variants without\n basecaller discrepancies") + coord_cartesian(xlim=c(1000, 29000))

#ggarrange(legend_plot, b_disc_plot, b_disc_plot2, nrow=3, heights=c(1, 0.4, 0.4), align = 'hv', labels = c('A', 'B', 'C'), label.y = c(1, 1.1, 1.1))

#Generate similar plot, but showing proportion of variants with basecaller discrepancies as a fraction of total number of variants
#Within the 200bp bins of the histogram
get_histogram_heights <- function(ggframe, poscol, bin_width){
  poscolind <- which(colnames(ggframe)==poscol)
  maxpos <- max(ggframe[,poscolind])
  maxrange <- bin_width*ceiling(maxpos/bin_width)
  nbins <- maxrange/bin_width
  #Get the numbers of variants in each bin and return as vector of length nbins
  return_vec <- rep(0, nbins)
  for(i in 1:nbins){
    maxbin <- i*bin_width
    minbin <- (maxbin - bin_width + 1)
    return_vec[i] <- sum(ggframe[,poscolind] %in% minbin:maxbin)
  }
  return(return_vec)
}

#Get heights data for both data sets and divide by total to get fraction
bin_width = 200
heights_disc <- get_histogram_heights(output_table_bdisc, poscol='POS', bin_width)
heights_nodisc <- get_histogram_heights(non_bdisc_variants, poscol='POS', bin_width)
heights_disc_frac <- heights_disc/(heights_disc + heights_nodisc)
ggtemp <- data.frame(xv=seq(from=0, to=30000-1, by=bin_width), yv=heights_disc_frac)
b_disc_props <- ggplot(ggtemp, aes(x = xv, y = yv)) + geom_bar(stat='identity', width = bin_width)+
  theme(legend.position = "none", axis.line=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  xlab("Locus") + ylab("Proportion of variants with\n basecaller discrepancies") + coord_cartesian(xlim=c(1000, 29000), ylim = c(0,1))

b_disc_props_noxaxis <- ggplot(ggtemp, aes(x = xv, y = yv)) + geom_bar(stat='identity', width = bin_width)+
  theme(legend.position = "none", axis.line=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.spacing.y = unit(0,"null"), plot.margin=unit(c(5.5,5.5,-1,5.5), "points")) + xlab(NULL) +
  ylab("Proportion of variants with\n basecaller discrepancies") + coord_cartesian(xlim=c(0, 30000), ylim = c(0,1), expand = FALSE, )

#Plot ORF chart and this plot together:
summary_figure <- ggarrange(disc_plot_legend, plot_grid(b_disc_props_noxaxis,
                            disc_plot_nolegend+theme(plot.margin=unit(c(-1,5.5,5.5,5.5), "points")),
                            ncol=1, align = 'v', axis='lr', rel_heights = c(0.3, 0.7)), 
                            heights = c(0.1, 1), labels = c('', '', ''), label.y = c(1, 1, 1), ncol=1,
                            legend.grob = disc_plot_legend, legend = 'top') 

pdf(file=paste0(main_project_directory, "blacklist_summary_distribution_final.pdf"), width = 11, height = 8)
print(summary_figure)
dev.off()

png(file=paste0(main_project_directory, "blacklist_summary_distribution_final.png"), width = 800, height = 600)
print(summary_figure)
dev.off()

#Does it look better to plot the bar chart below the ORF chart instead of splitting up the legend and ORF chart?
#ggarrange(legend_plot, b_disc_props, nrow=2, heights=c(1, 0.4), align = 'v', labels = c('A', 'B'), label.y = c(1, 1.1))

##END OF SECTION 17

##SECTION 18: Carry out Fisher Exact tests to generate tables of odds ratios, showing any significant correlations
#between different variant blacklists, or between variant blacklists and specific genomic regions of note.

#Define function to get the contingency table for odds ratio calculation (OR) of two sets of values that are subsets of a larger group
get_OR_table <- function(subsetA, subsetB, fullset){
  AB <- unique(intersect(subsetA, subsetB))
  Ab <- unique(setdiff(subsetA, subsetB))
  aB <- unique(setdiff(subsetB, subsetA))
  ab <- unique(setdiff(fullset, unique(c(subsetA, subsetB))))
  c_table <- matrix(c(length(AB), length(Ab), length(aB), length(ab)), nrow = 2, dimnames = list(c("In B", "Not in B"),c("In A", "Not in A")))
  return(c_table)
}

#Get OR values for all major combinations of blacklisted variant features
#The full set in this case should be the list of all variants, rather than all positions
FULLSET <- rowglue(variants_ordered[,c(2,1,3)])

SYSBIAS5 <- rowglue(output_table[which(output_table$SYSTEMATIC_BIAS_OVER_5_PERCENT) ,c(2,1,3)])
DEL20 <- rowglue(output_table[which(!is.na(output_table$HIGH_DEL_BASECALLERS)) ,c(2,1,3)])
VC1 <- rowglue(output_table[which(output_table$ONLY_CALLED_WITH_ONE_VARIANT_CALLER!=FALSE) ,c(2,1,3)])
BCnot4 <- rowglue(output_table[which((!output_table$BASECALLERS_AGREE_ON_VARIANT_CALL_MEDAKA) | (!output_table$BASECALLERS_AGREE_ON_VARIANT_CALL_NANOPOLISH)) ,c(2,1,3)])
#ZAF <- rowglue(output_table[which(output_table$ZERO_SUPPORTING_READS) ,c(2,1,3)])
BULL <- rowglue(output_table[which(output_table$IN_DEVESON_BLACKLIST) ,c(2,1,3)])
DEMAIO1 <- rowglue(output_table[which(output_table$IN_GOLDMAN_MASK_LIST) ,c(2,1,3)])
DEMAIO2 <- rowglue(output_table[which(output_table$IN_GOLDMAN_CAUTION_LIST) ,c(2,1,3)])
DEMAIO <- c(DEMAIO1, DEMAIO2)

#Create empty table for blacklist x blacklist odds ratios
colnames1 <- c("Blacklist features (columns) / Region features (rows)", 
               "Systematic bias >5%", "DEL AF >20%", "Only called in one variant caller",
               "Variant call not consistent across all four basecallers",
               "Bull et al.", "De Maio et al.")
rownames1 <- colnames1[2:length(colnames1)]

blacklist_blacklist_OR_table <- data.frame("Blacklist features / Region features"=rownames1,
                                           "Systematic bias >5%"=NA, "DEL AF >20%"=NA, "Only called in one variant caller"=NA,
                                           "Variant call not consistent across all four basecallers"=NA,
                                           "Bull et al."=NA, "De Maio et al."=NA)

#Fill in table values using list
list_col <- list(SYSBIAS5, DEL20, VC1, BCnot4, BULL, DEMAIO)
list_rows <- list_col
for(i in 1:length(list_rows)){
  for(j in 1:length(list_col)){
    if(i>j){
      FT_result <- fisher.test(get_OR_table(list_rows[[i]], list_col[[j]], FULLSET))
      asterisks <- ""
      if(FT_result$p.value < 0.05){asterisks <- "*"}
      if(FT_result$p.value < 0.005){asterisks <- "**"}
      if(FT_result$p.value < 0.0005){asterisks <- "***"}
      blacklist_blacklist_OR_table[i,j+1] <- paste0(round(FT_result$estimate, digits = 3), asterisks)
    }
  }
}

write.table(blacklist_blacklist_OR_table, file = paste0(blacklist_directory, "OR_table_blacklists_in_blacklists.xlsx"),
            sep = '\t', col.names = colnames1, row.names = FALSE)

#Get OR values for enrichment of blacklisted variants at genomic regions of interest
#These genomic regions are loaded/calculated in lines 92-189 of fisher_plot_suspect.R and contain no repeated values

#Get BED file for 100bp homopolymer flanks
pre_homo_bed <- homopolymer_bed
post_homo_bed <- homopolymer_bed
pre_homo_bed$STOP <- (pre_homo_bed$START - 1)
pre_homo_bed$START <- (pre_homo_bed$STOP - 99)
post_homo_bed$START <- (post_homo_bed$STOP + 1)
post_homo_bed$STOP <- (post_homo_bed$START + 99)
homo_100flanks_bed <- rbind(pre_homo_bed, post_homo_bed)

#Do similar for short 3bp homopolymer flanks
pre_homo_bed2 <- homopolymer_bed
post_homo_bed2 <- homopolymer_bed
pre_homo_bed2$STOP <- (pre_homo_bed2$START - 1)
pre_homo_bed2$START <- (pre_homo_bed2$STOP - 2)
post_homo_bed2$START <- (post_homo_bed2$STOP + 1)
post_homo_bed2$STOP <- (post_homo_bed2$START + 2)
homo_3flanks_bed <- rbind(pre_homo_bed2, post_homo_bed2)

#For each position, get the GC% in the surrounding 100bp
GC100 <- c()
for(i in 51:(length(cov_ref_vec)-50)){
  GC100 <- c(GC100,length(grep(cov_ref_vec[(i-50):(i+50)][-51], pattern = "[GC]")))
}
#Get a vector for all genomic positions that have greater than the upper quartile GC%
p99 <- quantile(GC100, 0.99)
GCgt99_vec <- which(GC100 > p99)

#Convert primer and homopolymer BED files into vectors of positions covered
homo_vec <- c()
homoflanks_vec <- c()
homoflanks_vec_short <- c()
primer_vec <- c()
for(i in 1:nrow(homopolymer_bed)){
  homo_vec <- c(homo_vec, homopolymer_bed[i,2]:homopolymer_bed[i,3])
}
for(i in 1:nrow(artic_primers)){
  primer_vec <- c(primer_vec, artic_primers[i,2]:artic_primers[i,3])
}
for(i in 1:nrow(homo_100flanks_bed)){
  homoflanks_vec <- c(homoflanks_vec, homo_100flanks_bed[i,2]:homo_100flanks_bed[i,3])
}
for(i in 1:nrow(homo_3flanks_bed)){
  homoflanks_vec_short <- c(homoflanks_vec_short, homo_3flanks_bed[i,2]:homo_3flanks_bed[i,3])
}
homo_vec <- setdiff(unique(homo_vec), (length(cov_ref_vec)+1):(length(cov_ref_vec)+200))
homoflanks_vec <- setdiff(unique(homoflanks_vec), (length(cov_ref_vec)+1):(length(cov_ref_vec)+200))
homoflanks_vec_short <- setdiff(unique(homoflanks_vec_short), (length(cov_ref_vec)+1):(length(cov_ref_vec)+200))
primer_vec <- unique(primer_vec)

#Read in non-canonical sgRNAs from Matt Parker and convert to vector
novel_shef=read.csv(sheffield_sgRNA_path)
novel_shef = novel_shef %>% separate(orf,c("novel","pos"),sep="_")
sgRNA_vec <- as.numeric(unique(novel_shef$pos))

#homo_vec
#homoflanks_vec_short
#homoflanks_vec
#GCgt99_vec
#primer_vec
#sgRNA_vec

#Need to strip allele data from variant lists to only retain position values
numstrip <- function(stringvec){return(gsub("[^0-9.-]", "", stringvec))}
#Also need to substitute full range of genomic loci for "FULLSET" which only covers loci at which variants are present
FULL_LOCI <- 1:length(cov_ref_vec)
#Also compare with non-blacklisted loci to check that enrichment is not merely a function of a position having called variants
non_blacklisted_variants <- variants_ordered[which(!(rowglue(variants_ordered[,1:3]) %in% rowglue(output_table[,1:3]))),]

colnames2 <- c("Blacklist features (columns) / Region features (rows)", 
               "Systematic bias >5%", "DEL AF >20%", "Only called in one variant caller", "Variant call not consistent across all four basecallers",
               "Bull et al.", "De Maio et al.", "Non-blacklisted variants")
rownames2 <- c("Homopolymers", "Homopolymer 3bp flanks", "Homopolymer 100bp flanks", "GCgt99p", "Artic primers", "sgRNAs")

blacklist_region_OR_table <- data.frame("Blacklist features / Region features"=rownames2,
                                        "Systematic bias >5%"=NA, "DEL AF >20%"=NA, "Only called in one variant caller"=NA, "Variant call not consistent across all four basecallers"=NA,
                                        "Bull et al."=NA, "De Maio et al."=NA, "Non-blacklisted variants"=NA)

#Fill in table values using list
list_col <- list(numstrip(SYSBIAS5), numstrip(DEL20), numstrip(VC1), numstrip(BCnot4),
                 numstrip(BULL), numstrip(DEMAIO), non_blacklisted_variants$POS)
list_rows <- list(homo_vec, homoflanks_vec_short, homoflanks_vec, GCgt99_vec, primer_vec, sgRNA_vec)
for(i in 1:length(list_rows)){
  for(j in 1:length(list_col)){
    FT_result <- fisher.test(get_OR_table(as.numeric(list_rows[[i]]), as.numeric(list_col[[j]]), FULL_LOCI))
    asterisks <- ""
    if(FT_result$p.value < 0.05){asterisks <- "*"}
    if(FT_result$p.value < 0.005){asterisks <- "**"}
    if(FT_result$p.value < 0.0005){asterisks <- "***"}
    blacklist_region_OR_table[i,j+1] <- paste0(round(FT_result$estimate, digits = 3), asterisks)
  }
}

write.table(blacklist_region_OR_table, file = paste0(blacklist_directory, "OR_table_blacklists_in_regions.xlsx"),
            sep = '\t', col.names = colnames2, row.names = FALSE)

##END OF SECTION 18