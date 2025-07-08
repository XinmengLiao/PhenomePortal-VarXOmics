## Conda environment configuration 

### Single gene query 
```bash
conda activate varxomics

# install r-base and relevant packages
conda install r-base=3.6 -y

# install g3viz by conda
conda install

# install epigraphdb and conflicted by CRAN
R
install.packages('epigraphdb')
install.packages('conflicted')
```

### Single gene query 
```bash
conda activate varxomics

# install r-base and relevant packages
conda install r-base=3.6 -y

# install g3viz by conda
conda install

# install maps by CRAN
R
install.packages('maps')
```

### For SVs and CNVs 
```bash
conda activate vep114

# AnnotSV and vcf2circos need to be installed by both conda and manually download 
conda install bioconda::annotsv
conda install bioconda::vcf2circos

# when running AnnotSV, AnnotSV comes from the manually downloaded ones, while vcf2circos comes from the conda one
```
