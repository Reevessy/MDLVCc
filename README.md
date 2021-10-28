# MDLVC converter
*v1.0*

A simple tool to convert JHU sequencing center's in-house MDLVC files to .csv files as well as rare germline variants annotating and filtering.

##### Please install the below R packages for the running environment:
```
install.packages(c("foreach", "doParallel", "expss", "rlang", "dplyr", "tidyr", "filesstrings"))
```

##### Input file format:
> *.MDLVC.vcf files (JHU's special VCF files that can't be read by vcfR)    
