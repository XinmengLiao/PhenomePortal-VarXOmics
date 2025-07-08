# vcf2circos need to be installed by both conda and manually download 
conda install bioconda::vcf2circos

# when running vcf2circos, using the one from the conda env path, so all the relevant files in conda env should be changed:
1) `vcf2circos/plotcategories/histogram.py` change to the one in `Scripts/vcf2circos/histogram.py`
2) `vcf2circos/config/json_config.py` change to the one in `Scripts/vcf2circos/json_config.py`
3) `vcf2circos/config/default_params.json` change to the one in `Scripts/vcf2circos/default_params.json`
4) `vcf2circos/config/__main__.py` change to the one in `Scripts/vcf2circos/__main__.py`
