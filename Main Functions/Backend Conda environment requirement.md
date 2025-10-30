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

Additonally, some of the vcf2circos files should be changed. 
1) `vcf2circos/plotcategories/histogram.py` change to the one in `Scripts/vcf2circos/histogram.py`
2) `vcf2circos/config/json_config.py` change to the one in `Scripts/vcf2circos/json_config.py`
3) `vcf2circos/config/default_params.json` change to the one in `Scripts/vcf2circos/default_params.json`
4) `vcf2circos/config/__main__.py` change to the one in `Scripts/vcf2circos/__main__.py`
