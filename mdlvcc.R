# Load R packages
lapply(c("foreach", "doParallel", "expss", "rlang", "dplyr", "tidyr", "filesstrings"), 
       require, character.only = TRUE)

# Set the folder containing current R script as working directory
setwd(".")
print(paste("Current working dir:", getwd()))
if (dir.exists("original_files")==FALSE){
  dir.create("original_files")
}
if (dir.exists("cache")==FALSE){
  dir.create("cache")
}

# Set cpu cores for parallel computing
numCores <- detectCores(all.tests = FALSE, logical = TRUE)



# Convert .vcf to .txt files
vcfs <- list.files(pattern=".MDLVC.vcf")
if (is_empty(vcfs)==FALSE){
  registerDoParallel(numCores)
  foreach (vcf = vcfs) %dopar% {
    allvariants <- readLines(vcf)
    allvariants <- as.data.frame(allvariants)
    allvariants <-allvariants[-grep('#',allvariants$allvariants),]
    filename <- paste0(vcf, ".txt")
    write.table(allvariants, file = filename, sep = "\t", row.names = FALSE)
    file.move(vcf, "./original_files/a_vcf_files", overwrite=TRUE)
  }
  rm("vcfs")
}

# Convert .txt to .csv files
txts <- list.files(pattern=".MDLVC.vcf.txt")
if (is_empty(txts)==FALSE){
  registerDoParallel(numCores)
  foreach (txt = txts) %dopar% {
    allvariants <- read.table(txt,header=T,sep="\t");
    is.data.frame(allvariants)
    if(nrow(allvariants)>0){
      allvariants <- separate(allvariants, "x",
                              into=c("CHROM","POS","ID", "REF", "ALT", "QUAL", 
                                     "FILTER", "INFO"),sep="\t")
      allvariants$Sample <- rep(txt,nrow(allvariants))
      allvariants$Sample <- gsub('.MDLVC.vcf.txt','',allvariants$Sample)
      allvariants$Variant <- paste0(allvariants$CHROM, "-", allvariants$POS, "-", allvariants$REF, "-", allvariants$ALT)
    }
    file.move(txt, "./cache/a_txt_files", overwrite=TRUE)
    filename <- paste0(txt, ".csv")
    write.table(allvariants, file=filename, sep=",", row.names = FALSE)
  }
  rm("txts")
}

# Combine all .csv files and first filtering
csvs <- list.files(pattern='.MDLVC.vcf.txt.csv')
allvariants <- lapply(csvs, read.csv) 
allvariants <- do.call(rbind.data.frame, allvariants)

chr9 <- filter(allvariants, CHROM == "chr9")
jak2 <- filter(chr9, POS >= 4985033 & POS <= 5128183)
jak2$Gene <- rep("JAK2",nrow(jak2))

chr11 <- filter(allvariants, CHROM == "chr11")
atm <- filter(chr11, POS >= 108093211 & POS <= 108239829)
atm$Gene <- rep("ATM",nrow(atm))

variants <- rbind(jak2, atm)
variants <- variants[c("Sample", "Variant", "Gene", "QUAL", "FILTER", "INFO")]

file.move(csvs, "./cache/b_csv_files", overwrite=TRUE)
write.table(variants, file="variants_in_ATM&JAK2.csv", sep=",", row.names = FALSE)
rm("csvs", "allvariants", "chr9", "chr11", "jak2", "atm", "variants")



# Manually download GnomAD referece data

# Organize variant ref info
reffiles <- list.files(pattern='.ref.csv')
for (reffile in reffiles) {
  csv1 <- read.csv(reffile, header = TRUE, na.strings=c("","NA"))
  colnames(csv1) <- c("Chromosome","Position","rsIDs","Reference","Alternate","Source","Filters_exomes","Filters_genomes","Transcript","HGVS_Consequence",
                      "Protein_Consequence","Transcript_Consequence","VEP_Annotation","ClinVar_Clinical_Significance","ClinVar_Variation_ID","Flags",
                      "Allele_Count","Allele_Number","Allele_Frequency","Homozygote_Count","Hemizygote_Count",
                      "Allele_Count_African_or_African_American","Allele_Number_African_or_African_American","Homozygote_Count_African_or_African_American","Hemizygote_Count_African_or_African_American",
                      "Allele_Count_Latino_or_Admixed_American","Allele_Number_Latino_or_Admixed_American","Homozygote_Count_Latino_or_Admixed_American","Hemizygote_Count_Latino_or_Admixed_American",
                      "Allele_Count_Ashkenazi_Jewish","Allele_Number_Ashkenazi_Jewish","Homozygote_Count_Ashkenazi_Jewish","Hemizygote_Count_Ashkenazi_Jewish",
                      "Allele_Count_East_Asian","Allele_Number_East_Asian","Homozygote_Count_East_Asian","Hemizygote_Count_East_Asian",
                      "Allele_Count_European_Finnish","Allele_Number_European_Finnish","Homozygote_Count_European_Finnish","Hemizygote_Count_European_Finnish",
                      "Allele_Count_European_non_Finnish","Allele_Number_European_non_Finnish","Homozygote_Count_European_non_Finnish","Hemizygote_Count_European_non_Finnish",
                      "Allele_Count_Other","Allele_Number_Other","Homozygote_Count_Other","Hemizygote_Count_Other",
                      "Allele_Count_South_Asian","Allele_Number_South_Asian","Homozygote_Count_South_Asian","Hemizygote_Count_South_Asian")
  csv1$ENSG_ID <- rep(reffile,nrow(csv1))
  csv1$ENSG_ID <- gsub('.ref.csv', '', csv1$ENSG_ID)
  csv1$Variant <- paste0(csv1$Chromosome, "-", csv1$Position, "-", csv1$Reference, "-", csv1$Alternate)
  csv1 <- mutate(csv1, Global_AF=Allele_Count/Allele_Number)
  csv1 <- mutate(csv1, African_AF=Allele_Count_African_or_African_American/Allele_Number_African_or_African_American)
  csv1 <- mutate(csv1, Latino_AF=Allele_Count_Latino_or_Admixed_American/Allele_Number_Latino_or_Admixed_American)
  csv1 <- mutate(csv1, Ashkenazi_Jewish_AF=Allele_Count_Ashkenazi_Jewish/Allele_Number_Ashkenazi_Jewish)
  csv1 <- mutate(csv1, East_Asian_AF=Allele_Count_East_Asian/Allele_Number_East_Asian)
  csv1 <- mutate(csv1, European_Finnish_AF=Allele_Count_European_Finnish/Allele_Number_European_Finnish)
  csv1 <- mutate(csv1, European_Non_Finnish_AF=Allele_Count_European_non_Finnish/Allele_Number_European_non_Finnish)
  csv1 <- mutate(csv1, Other_AF=Allele_Count_Other/Allele_Number_Other)
  csv1 <- mutate(csv1, South_Asian_AF=Allele_Count_South_Asian/Allele_Number_South_Asian)
  csv1 <- csv1[c("Chromosome", "rsIDs", "Source", "Filters_exomes", "Filters_genomes", "HGVS_Consequence", "Protein_Consequence", 
                 "Transcript_Consequence", "VEP_Annotation", "ClinVar_Clinical_Significance", "Flags", "ENSG_ID", "Variant", "Global_AF", "African_AF", "Latino_AF", 
                 "Ashkenazi_Jewish_AF", "East_Asian_AF", "European_Finnish_AF", "European_Non_Finnish_AF", "South_Asian_AF", "Other_AF")]
  filename <- paste0(reffile, ".af.csv")
  write.table(csv1, file=filename, sep=",", row.names = FALSE)
  file.move(reffile, "./original_files/b_variant_ref", overwrite=TRUE)
}
csvs <- list.files(pattern='.ref.csv.af.csv')
refvariants <- lapply(csvs, read.csv) 
refvariants <- do.call(rbind.data.frame, refvariants)
write.table(refvariants, file="ref_variants_in_ATM&JAK2.csv", sep=",", row.names = FALSE)
file.move(csvs, "./original_files/b_variant_ref", overwrite=TRUE)
rm("csv1", "refvariants", "csvs", "filename", "reffile", "reffiles")csv2 <- filter(csv1, FILTER == "PASS" & Total_Depth >=50)

# Annotate variants found in sequencing
csv1 <- read.csv("variants_in_ATM&JAK2.csv")
csv1$Variant <- gsub('chr', '', csv1$Variant)
csv2 <- read.csv("ref_variants_in_ATM&JAK2.csv")
csv3 <- add_columns(csv1, csv2, by="Variant")
write.table(csv3, file="all_variants_in_ATM&JAK2.csv", sep=",", row.names = FALSE)
file.move(c("variants_in_ATM&JAK2.csv", "ref_variants_in_ATM&JAK2.csv"), "./cache/c_variant_formatting", overwrite=TRUE)
rm("csv1", "csv2", "csv3")

