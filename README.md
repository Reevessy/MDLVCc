# MDLVCc
A simple tool for decoding JHU sequencing center's in-house MDLVC files as well as variants annotating.    
*Current version: v1.0.0*

##### Please install the below R packages for the running environment:
```
install.packages(c("foreach", "doParallel", "expss", "rlang", "dplyr", "tidyr", "filesstrings"))
```

##### Input file format:
> *.MDLVC.vcf files (JHU's special VCF files that can't be read by vcfR)    
