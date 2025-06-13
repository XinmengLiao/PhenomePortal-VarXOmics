conda activate vep

zgrep -E "^#|PASS" sample.vcf.gz | bgzip > sample_PASS.vcf.gz 

echo "vep cache: 113, vep version: 113."

vep --cache -dir_cache vep_cache \
--offline \
--cache_version 113 \
--fork 48 \
--format vcf \
--dir_plugins vep_cache/VEP_plugins \
-i sample_PASS.vcf.gz \
-o sample_vep_annotated.vcf.gz \
--force_overwrite \
--compress_output bgzip \
--assembly GRCh38 \
--symbol --vcf --check_existing --variant_class \
--sift b --polyphen b \
--synonyms vep_cache/homo_sapiens_refseq/113_GRCh38/chr_synonyms.txt \
--hgvs --refseq \
--fasta vep_cache/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
--canonical \
--pick --pick_order mane_select,rank \
--af --af_gnomade --af_gnomadg --max_af \
--custom vep_cache/vep_custom/clinvar_20240611_PLPC_new_CPLP.vcf.gz,ClinVar,vcf,exact,0,ID,CLNSIG,CLNDN,CLNHGVS,CLNSIGINCL,CLNVC,GENEINFO,CLNDISDB,CLNSIGCONF,CLNREVSTAT,CLNDNINCL \
--custom vep_cache/vep_custom/Whole_out_sorted_2024.vcf.gz,Database,vcf,exact,0,Type,SZAID \
--plugin dbscSNV,vep_cache/dbscSNV1.1/dbscSNV1.1_GRCh38.txt.gz \
--plugin REVEL,file=vep_cache/REVEL/new_tabbed_revel_grch38.tsv.gz \
--plugin dbNSFP,vep_cache/dbNSFP4.7a/dbNSFP4.7a_grch38.gz,REVEL_score,REVEL_rankscore,BayesDel_addAF_score,BayesDel_addAF_rankscore,BayesDel_addAF_pred,BayesDel_noAF_score,BayesDel_noAF_rankscore,BayesDel_noAF_pred \
--plugin SpliceAI,snv=vep_cache/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=vep_cache/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \
--plugin AlphaMissense,file=vep_cache/AlphaMissense/AlphaMissense_hg38.tsv.gz,cols=all \
--plugin LoF,loftee_path:loftee/,human_ancestor_fa:loftee/human_ancestor.fa.gz,conservation_file:loftee/loftee.sql,gerp_bigwig:loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
--verbose