#%%Cell 1 Define the path for parameters and databases
# ==============================================================================
import sys
import pandas as pd
import numpy as np
import os
import gzip
import paramiko
import re
from scp import SCPClient
import time
from itertools import islice
from datetime import datetime
import pickle
from collections import defaultdict

all_start_time = time.time()
# ====== 1) setting up paths and running parameters ======
path = 'sysmed'  # local or sysmed

# User input filtration parameters
# user_af_clinvar = 1
# user_af_predict = 0.05
# user_ada_score = 0.6
# user_rf_score = 0.6
# user_revel_score=0.75
# user_spliceai_al=0.5
# user_spliceai_dg=0.5
# user_spliceai_dl=0.5
# user_spliceai_ag=0.5
# user_bayesdel_addaf_score=0.0692655
# user_bayesdel_noaf_score=-0.0570105
# user_am_classification=['likely_pathogenic','ambigous']
# user_am_pathogenicity=0.564
# user_gender= sys.argv[3] # Male, Female, blank
# user_clinvar= ["Pathogenic","Likely_pathogenic","Uncertain_significance","Conflicting_classifications_of_pathogenicity"]
# user_acmg_classification = ["Pathogenic","Likely_pathogenic","Uncertain_significance","Benign","Likely_benign"]

user_gender= sys.argv[3] # Male, Female, blank
user_af_clinvar = sys.argv[4]
user_af_clinvar = float(user_af_clinvar)
user_af_predict = sys.argv[5]
user_af_predict = float(user_af_predict)
user_ada_score = sys.argv[6]
user_ada_score = float(user_ada_score)
user_rf_score = sys.argv[7]
user_rf_score = float(user_rf_score)
user_revel_score = sys.argv[8]
user_revel_score = float(user_revel_score)
user_spliceai_al = sys.argv[9]
user_spliceai_al = float(user_spliceai_al)
user_spliceai_dg = sys.argv[10]
user_spliceai_dg = float(user_spliceai_dg)
user_spliceai_dl = sys.argv[11]
user_spliceai_dl = float(user_spliceai_dl)
user_spliceai_ag = sys.argv[12]
user_spliceai_ag = float(user_spliceai_ag)
user_bayesdel_addaf_score = sys.argv[13]
user_bayesdel_addaf_score = float(user_bayesdel_addaf_score)
user_bayesdel_noaf_score = sys.argv[14]
user_bayesdel_noaf_score = float(user_bayesdel_noaf_score)
user_am_classification = sys.argv[15]
user_am_pathogenicity = sys.argv[16]
user_am_pathogenicity = float(user_am_pathogenicity)
user_clinvar= sys.argv[17]
user_acmg_classification = sys.argv[18]

# config file of setting paths
config = {
    "local": {
        "fileName": sys.argv[1],
        "output_file": sys.argv[2],
        "clinvar_vcf_path": "/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/ClinVar/clinvar_20250504.vcf.gz",
        "GenCC_path" : "/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/GenCC/GenCCforGene.txt",
        "eqtl_catalog_file" : "/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/xQTL/eQTL_Catalogue/permuted_results/eqtl.catalogue.permuted.known.position.txt",
        "eqtl_gtex_file" : "/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/xQTL/GTEx_Analysis_v10_eQTL_updated/GTEx.all.tissue.known.gene.txt",
        "gwas_file" : "/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/gwas_catalog_v1.0-associations_e113_r2025-04-28.tsv",
        "pgx_file": "/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/PharmGKB/clinical_annotation_combined.txt",
        "otg_file": "/Users/xinmengliao/Documents/Project/20250516_Webserver/Datasets/xQTL/Open_Targets_Genetics/OTG.scores.txt"
        
    },
    "sysmed": {
        "fileName": sys.argv[1],
        "outFile": sys.argv[2],
        "clinvar_vcf_path": "/mnt/storage_pool/Genomics/VarXOmics/Databases/clinvar_20250504.vcf.gz",
        "GenCC_path" : "/mnt/storage_pool/Genomics/VarXOmics/Databases/GenCCforGene.txt",
        "eqtl_catalog_file" : "/mnt/storage_pool/Genomics/VarXOmics/Databases/eqtl.catalogue.permuted.known.position.txt",
        "eqtl_gtex_file" : "/mnt/storage_pool/Genomics/VarXOmics/Databases/GTEx.all.tissue.known.gene.txt",
        "gwas_file" : "/mnt/storage_pool/Genomics/VarXOmics/Databases/gwas_catalog_v1.0-associations_e113_r2025-04-28.tsv",
        "pgx_file": "/mnt/storage_pool/Genomics/VarXOmics/Databases/PharmGKB_clinical_annotation_combined.txt",
        "otg_file": "/mnt/storage_pool/Genomics/VarXOmics/Databases/OTG.scores.txt"
    }
}

if path not in config:
    raise ValueError(f"unknown path: {path}")

cfg = config[path]

# ====== 3) Reading necessary files ======

def read_db_file(filepath, encoding="ISO-8859-1", sep="\t", fillna_str="No info", drop_allna_cols=False):
    """help to formally read files"""
    df = pd.read_csv(filepath, sep=sep, encoding=encoding)
    if drop_allna_cols:
        df = df.dropna(axis=1, how='all')
    df = df.replace(np.nan, fillna_str)
    return df


eqtl_catalog = pd.read_csv(cfg["eqtl_catalog_file"], sep="\t", header=0) 
eqtl_catalog.columns = [f"eQTL_Catalog_{col}" for col in eqtl_catalog.columns]
eqtl_catalog_set = set(eqtl_catalog['eQTL_Catalog_variant'].to_list())

eqtl_gtex = pd.read_csv(cfg["eqtl_gtex_file"], sep="\t", header=0)
eqtl_gtex.columns = [f"eQTL_GTEx_{col}" for col in eqtl_gtex.columns]
eqtl_gtex_set = set(eqtl_gtex['eQTL_GTEx_variant_id'].to_list())

# ====== 4) Decision of whether running GWAS/Pharmaco======
gwas = pd.read_csv(cfg["gwas_file"], sep="\t", quotechar='"')
gwas['Risk.allele'] = gwas['STRONGEST SNP-RISK ALLELE'].str.split('-').str[1]
gwas_set = set(gwas['SNPS'].to_list())

pharmgkb_data = pd.read_csv(cfg["pgx_file"], sep="\t")
pgx_set = set(pharmgkb_data['Variant.Haplotypes'].to_list())

otg_database = pd.read_csv(cfg["otg_file"], sep="\t")

#%%Cell 2 start the analysis program for the vep annotated vcf files, generated report A B C and D
# ==============================================================================
start_time = time.time()
# col_map is used to build a mapping once colNames_CSQ is obtained
col_map = {}  # Will be populated after parsing the "ID=CSQ" line
headings = []

## For REVEL scores refinement from 0.460&0.460&.&. to 0.460
def extract_decimal_from_string(s):
    if not s or not isinstance(s, str):
        return None
    matches = re.findall(r"\d+\.\d+", s)
    if not matches:
        return None
    return matches[0]


# Reading vcf.gz file 
file = gzip.open(cfg["fileName"],'rt')
tLine = file.readline()
i = 0
reportA,reportgwas, reportpgx, reporteqtl = [], [], [], []

while tLine:
    # remove the newline character
    tLine = tLine.rstrip('\n')
    # split the current line
    iContent = tLine.split('\t')
    i += 1
    ##get the content from VCF annotation header
    if tLine.startswith('#'):
        if 'ID=CSQ' in tLine:
            annoText = iContent[0].split('Format: ')
            colNames_CSQ = annoText[1].replace('">','')
            colNames_CSQ = colNames_CSQ.split('|')
            # construct col_map for all use
            col_map = { name: idx for idx, name in enumerate(colNames_CSQ) }
        elif tLine.startswith('#CHROM'):
            headings = iContent
            print(headings)
        # directly goes into next line
        tLine = file.readline()
        #print(tLine)
        continue
    if not headings:
        tLine = file.readline()
        #print(tLine)
        continue
    
    iText = [s for s in iContent[headings.index('INFO')].split(';') if 'CSQ=' in s]
    iText = iText[0].replace('CSQ=','').split(',')
    
    saveFlag1, saveFlaggwas, saveFlagPGx, saveFlageqtl = False, False, False, False
    
    for j in range(0,len(iText)):
        jText = iText[j].split('|')
        # fixing REVEL score 
        if 'REVEL' in col_map and 'REVEL_score' in col_map:
            revel_idx = col_map['REVEL']
            revel_score_idx = col_map['REVEL_score']
            if jText[revel_idx] == '':
                parsed_value = extract_decimal_from_string(jText[revel_score_idx])
                if parsed_value is not None:
                    jText[revel_idx] = parsed_value
                            
        # 1) ClinVar and prediction filtered based on users' needs 
        if 'MAX_AF' in col_map:
            max_af_val = jText[col_map['MAX_AF']]
            try:
                max_af_numeric = float(max_af_val) if max_af_val != '' else 0.0
            except (ValueError, TypeError):
                max_af_numeric = 0.0
            jText[col_map['MAX_AF']] = max_af_numeric

        if jText[col_map['ClinVar_CLNSIG']] != "" and max_af_numeric < user_af_clinvar:
            saveFlag1 = "Keep"
        
        if max_af_numeric < user_af_predict:   
            predicted_impact = (
                    ('IMPACT' in col_map and jText[col_map['IMPACT']] == 'HIGH')
                    or ('ada_score' in col_map and jText[col_map['ada_score']] != '' and float(jText[col_map['ada_score']]) > user_ada_score)
                    or ('rf_score' in col_map and jText[col_map['rf_score']] != '' and float(jText[col_map['rf_score']]) > user_rf_score)
                    or ('REVEL' in col_map and jText[col_map['REVEL']] != '' and float(jText[col_map['REVEL']]) > user_revel_score)
                    or ('SpliceAI_pred_DS_AL' in col_map and jText[col_map['SpliceAI_pred_DS_AL']] != '' and float(jText[col_map['SpliceAI_pred_DS_AL']])>user_spliceai_al)
                    or ('SpliceAI_pred_DS_DG' in col_map and jText[col_map['SpliceAI_pred_DS_DG']] != '' and float(jText[col_map['SpliceAI_pred_DS_DG']])>user_spliceai_dg)
                    or ('SpliceAI_pred_DS_DL' in col_map and jText[col_map['SpliceAI_pred_DS_DL']] != '' and float(jText[col_map['SpliceAI_pred_DS_DL']])>user_spliceai_dl)
                    or ('SpliceAI_pred_DS_AG' in col_map and jText[col_map['SpliceAI_pred_DS_AG']] != '' and float(jText[col_map['SpliceAI_pred_DS_AG']])>user_spliceai_ag)
                    or ('BayesDel_addAF_score' in col_map and jText[col_map['BayesDel_addAF_score']] != '' and float(jText[col_map['BayesDel_addAF_score']])>user_bayesdel_addaf_score)
                    or ('BayesDel_noAF_score' in col_map and jText[col_map['BayesDel_noAF_score']] != '' and float(jText[col_map['BayesDel_noAF_score']])>user_bayesdel_noaf_score)
                    or ('am_class' in col_map and jText[col_map['am_class']] in user_am_classification
                        and 'am_pathogenicity' in col_map
                        and jText[col_map['am_pathogenicity']] != ''
                        and float(jText[col_map['am_pathogenicity']])>user_am_pathogenicity))
            if predicted_impact:
                saveFlag1 = "Keep"
        
        # 2) judge PGx, GWAS 
        if 'Existing_variation' in col_map:
            ex_variation_val = jText[col_map['Existing_variation']]
            ex_variants = ex_variation_val.split('&')
            if ex_variants and ex_variants[0] in pgx_set:
                saveFlagPGx = "PGx_var"
            if ex_variants and ex_variants[0] in gwas_set:
                saveFlaggwas = "GWAS_var"

        # 3) judge eqtl
        var_info = f'{iContent[0]}_{iContent[1]}_{iContent[3]}_{iContent[4]}'
        if var_info in eqtl_catalog_set or var_info in eqtl_gtex_set:
            saveFlageqtl = "eqtl_var"
    
    # after for j in range(len(iText)) loop, if saveFlag1/2/3/4/5 has value, append the line to respective report
    if saveFlag1:
        reportA.append(tLine)
    if saveFlagPGx:
        reportpgx.append(tLine)
    if saveFlageqtl:
        reporteqtl.append(tLine)
    if saveFlaggwas:
        reportgwas.append(tLine)


    # print progress every 1000000 lines
    if i % 1000000 == 0:
        print(f"{i} lines processed!")

    # read the next line
    tLine = file.readline()


file.close()
print('Cell 2 VEP annotated File processing done! Now start to map GeneDB and DiseaseDB')
end_time = time.time()
print("Total processing time: {:.2f} seconds".format(end_time - start_time))

# Manage text file into a dataframe
base_vcf_columns = headings
#base_vcf_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT','Sample_Info']

csq_columns = [''] * len(col_map)  
for col_name, col_index in col_map.items():
    csq_columns[col_index] = col_name

all_output_columns = base_vcf_columns + csq_columns

def process_vcf_line_to_expanded_format(vcf_line, col_map, csq_columns):
    expanded_rows = []

    fields = vcf_line.split('\t')
    info_field = fields[7]  # INFO
    csq_info = [s for s in info_field.split(';') if 'CSQ=' in s]
    
    if not csq_info:
        return expanded_rows
    
    csq_data = csq_info[0].replace('CSQ=', '').split(',')
    
    for transcript_annotation in csq_data:
        transcript_fields = transcript_annotation.split('|')
        
        while len(transcript_fields) < len(csq_columns):
            transcript_fields.append('')
        
        expanded_row = fields + transcript_fields[:len(csq_columns)]
        expanded_rows.append(expanded_row)
    
    return expanded_rows

expanded_reportA = []
for i, vcf_line in enumerate(reportA):
    expanded_rows = process_vcf_line_to_expanded_format(vcf_line, col_map, csq_columns)
    expanded_reportA.extend(expanded_rows)
if expanded_reportA:
    expanded_reportA = pd.DataFrame(expanded_reportA, columns=all_output_columns)

# Keeps the rows with Users defined prediction and ClinVar criteria
def contains_any_or_empty(x, keywords):
    if str(x).strip() == "":
        return True
    return any(kw in str(x) for kw in keywords)

expanded_reportA = expanded_reportA[expanded_reportA["ClinVar_CLNSIG"].apply(lambda x: contains_any_or_empty(x, user_clinvar))]



# Zygosity
# Get which colums contains the sampe information (start with 0/0 or 1/1)
flexible_pattern = re.compile(r'^[01\.][/|][01\.]')
matched_colnames = [] 
for col in expanded_reportA.columns:
    pattern = expanded_reportA[col].iloc[0]    
    if flexible_pattern.match(str(pattern)):
        matched_colnames.append(col) 

print("匹配的列名:", matched_colnames)


def parse_sample_info(sample_info, chrom, gender):

    gt_str = str(sample_info).split(':')[0]  

    if gender == 'Male' and chrom in ['X', 'chrX']:
        if gt_str in ['0/1', '1/0', '0|1', '1|0']:
            return 'Hemizygous'
        elif gt_str in ['1/1', '1|1', '0/0', '0|0']:
            return 'Homozygous' 
        else:   
            return 'Unknown'
        
    if gender == '' and chrom in ['X', 'chrX']:
        if gt_str in ['0/1', '1/0', '0|1', '1|0']:
            return 'Heterozygous or Hemizygous, Gender not defined.'
        else:
            return 'Unknown'

    if chrom not in ['X', 'chrX']:
        if gt_str in (('0/1', '1/0', '0|1', '1|0')):
            return 'Heterozygous'
        elif gt_str in (('1/1', '1|1', '0/0', '0|0')):
            return 'Homozygous'
        elif gt_str in ['./.', '.|.']:
            return 'Missing'
        else:
            return 'Unknown'

def parse_genotype_info(sample_info, ref, alt):

    gt_str = str(sample_info).split(':')[0]  
    gt_ref = str(ref)
    gt_alt = str(alt)

    if gt_str == '0/1':
        return f'{gt_ref}/{gt_alt}'
    elif gt_str == '1/0':
        return f'{gt_alt}/{gt_ref}'
    elif gt_str == '0|1':
        return f'{gt_ref}|{gt_alt}'
    elif gt_str == '1|0':
        return f'{gt_alt}|{gt_ref}'
    elif gt_str == '1|1':
        return f'{gt_alt}|{gt_alt}'
    elif gt_str == '1/1':
        return f'{gt_alt}/{gt_alt}'
    elif gt_str == '0|0':
        return f'{gt_ref}|{gt_ref}'
    elif gt_str == '0/0':
        return f'{gt_ref}/{gt_ref}'
    elif gt_str == './.':
        return 'Missing'
    elif gt_str == '.|.':
        return 'Missing'
    else:   
        return 'Unknown'


expanded_reportA['Zygosity'] = expanded_reportA.apply(lambda row: parse_sample_info(row[matched_colnames[0]], row['#CHROM'], user_gender), axis=1)
expanded_reportA['Genotype'] = expanded_reportA.apply(lambda row: parse_genotype_info(row[matched_colnames[0]], row['REF'], row['ALT']), axis=1)


#%% Cell 3 GeneBe ACMG classification ======================================
import genebe as gnb

small_df = expanded_reportA.loc[:,["#CHROM","POS","REF","ALT"]]
small_df = small_df.rename(columns={"#CHROM":"chr","POS":"pos","REF":"ref","ALT":"alt"})
unique_small_df = small_df.drop_duplicates()

try:
    annotated_df = gnb.annotate(
        unique_small_df,
        genome='hg38',
        use_ensembl=False,
        use_refseq=True,
        flatten_consequences=True,
        output_format="dataframe"
    )
except Exception as e:
    print(f"GeneBe annotation failed: {str(e)}")
    # create an empty DataFrame to keep the structure
    annotated_df = pd.DataFrame(columns=[
        'chr', 'pos', 'ref', 'alt', 'gene_symbol', 
        'acmg_score', 'acmg_classification', 'acmg_criteria'
    ])
annotated_df = annotated_df.rename(columns={"chr":"#CHROM","pos":"POS","ref":"REF","alt":"ALT"})

required_columns = ["#CHROM","POS","REF","ALT","gene_symbol",
                   "acmg_score","acmg_classification",'acmg_criteria']
for col in required_columns:
    if col not in annotated_df.columns:
        annotated_df[col] = None  # add missing columns

small_annotate_all = annotated_df[required_columns].rename(columns={'gene_symbol': 'SYMBOL'})
# merge automatically handles missing values
expanded_reportA = pd.merge(
    expanded_reportA, 
    small_annotate_all, 
    how="left", 
    on=["#CHROM","POS","REF","ALT","SYMBOL"]
)

#%% Cell 4 Output python managed file and extract unique genes ======================================

# Output expended_reportA
# expanded_reportA['MAX_AF'] = pd.to_numeric(expanded_reportA['MAX_AF'], errors='coerce')
# expanded_reportA = expanded_reportA[(expanded_reportA['MAX_AF'].isna()) | (expanded_reportA['MAX_AF'] < user_af_predict)]
expanded_reportA['rsID'] = expanded_reportA['Existing_variation'].str.extract(r'(rs\d+)', expand=False)
expanded_reportA['variant_info'] = expanded_reportA[['#CHROM', 'POS', 'REF', 'ALT']].fillna('').astype(str).agg('_'.join, axis=1)
expanded_reportA = expanded_reportA[expanded_reportA['acmg_classification'].isin(user_acmg_classification)]
expanded_reportA = expanded_reportA.drop_duplicates()

# Python output Statistic
print("\nPython output Statistic")
print(f"Original rows: {i}")
print(f"Filtered and output rows: {len(expanded_reportA)}")

expanded_reportA.to_csv(cfg["output_file"], sep="\t", index=False, quoting=3)

#extract the genes
gene_list = expanded_reportA['SYMBOL'].unique()
formatted_genes = ",".join([f"'{gene.strip()}'" for gene in gene_list])
with open(cfg["output_file"].replace('.txt', '_genes.txt'), 'w') as f:
    f.write(formatted_genes)


#%% Cell 6 GWAS
print("Now generating GWAS results")
csq_data = []
gwas_data = [line.split('\t') for line in reportgwas]
gwas_df = pd.DataFrame(gwas_data, columns=base_vcf_columns)
for info in gwas_df['INFO']:
    csq_entry = [s for s in info.split(';') if 'CSQ=' in s]
    if csq_entry:
        csq_content = csq_entry[0].replace('CSQ=', '')
        csq_records = csq_content.split(',') 
        
        for csq_rec in csq_records:
            csq_split = csq_rec.split('|')
            csq_data.append(csq_split)
    else:
        csq_split = [''] * len(base_vcf_columns)
        csq_data.append(csq_split)

csq_df = pd.DataFrame(csq_data, columns=csq_columns)
gwas_merged_df = pd.concat([gwas_df.reset_index(drop=True), csq_df.reset_index(drop=True)], axis=1)

gwas_merged_df['combine1'] = gwas_merged_df['Existing_variation'].str.split('&').str[0] + '_' + gwas_merged_df['REF']
gwas_merged_df['combine2'] = gwas_merged_df['Existing_variation'].str.split('&').str[0] + '_' + gwas_merged_df['ALT']

gwas['combine'] = gwas['SNPS'] + '_' + gwas['Risk.allele']
gwas_combine_set1 = set(gwas_merged_df['combine1'])
gwas_combine_set2 = set(gwas_merged_df['combine2'])
mask = (gwas['combine'].isin(gwas_combine_set1)) | (gwas['combine'].isin(gwas_combine_set2))
gwas_filtered_df = gwas[mask].reset_index(drop=True)

variant_set = set(expanded_reportA['Existing_variation'].str.split('&').str[0])
gwas_filtered_user = gwas_filtered_df[gwas_filtered_df['SNPS'].isin(variant_set)].reset_index(drop=True)

gwas_output_file = cfg["output_file"].replace("txt", "gwas.txt")
gwas_output_file_user = cfg["output_file"].replace("txt", "gwas_filtered.txt")
gwas_filtered_df.to_csv(gwas_output_file, sep="\t", index=False, quoting=3)
gwas_filtered_user.to_csv(gwas_output_file_user, sep="\t", index=False, quoting=3)


#%% Cell 7 PGx
print("Now generating PGx results")
csq_data = []
pgx_data = [line.split('\t') for line in reportpgx]
pgx_df = pd.DataFrame(pgx_data, columns=base_vcf_columns)
for info in pgx_df['INFO']:
    csq_entry = [s for s in info.split(';') if 'CSQ=' in s]
    if csq_entry:
        csq_content = csq_entry[0].replace('CSQ=', '')
        csq_records = csq_content.split(',') 
        
        for csq_rec in csq_records:
            csq_split = csq_rec.split('|')
            csq_data.append(csq_split)
    else:
        csq_split = [''] * len(base_vcf_columns)
        csq_data.append(csq_split)

csq_df = pd.DataFrame(csq_data, columns=csq_columns)
pgx_merged_df = pd.concat([pgx_df.reset_index(drop=True), csq_df.reset_index(drop=True)], axis=1)


pgx_merged_df['Genotype'] = pgx_merged_df.apply(lambda row: parse_genotype_info(row[matched_colnames[0]], row['REF'], row['ALT']), axis=1)

pgx_merged_df['sorted_genotype'] = (
    pgx_merged_df['Genotype']
    .str.replace(r'[\|/]', '', regex=True) 
    .apply(lambda x: ''.join(sorted(x))) 
)

pgx_merged_df['combine'] = pgx_merged_df['Existing_variation'] + '_' + pgx_merged_df['sorted_genotype']

pharmgkb_data['combine'] = pharmgkb_data['Variant.Haplotypes'] + '_' + pharmgkb_data['Genotype.Allele']
pgx_combine_set = set(pgx_merged_df['combine'])
mask = (pharmgkb_data['combine'].isin(pgx_combine_set))
pgx_filtered_df = pharmgkb_data[mask].reset_index(drop=True)

# Only the ones overlapped with expended reportA
variant_set = set(expanded_reportA['Existing_variation'].str.split('&').str[0])
pgx_filtered_user = pgx_filtered_df[pgx_filtered_df['Variant.Haplotypes'].isin(variant_set)].reset_index(drop=True)

pgx_output_file = cfg["output_file"].replace("txt", "pgx.txt")
pgx_output_file_user = cfg["output_file"].replace("txt", "pgx_filtered.txt")
pgx_filtered_df.to_csv(pgx_output_file, sep="\t", index=False, quoting=3)
pgx_filtered_user.to_csv(pgx_output_file_user, sep="\t", index=False, quoting=3)


#%% Cell 8 eqtl
print("Now generating eqtl results")
csq_data = []
eqtl_data = [line.split('\t') for line in reporteqtl]
eqtl_df = pd.DataFrame(eqtl_data, columns=base_vcf_columns)
for info in eqtl_df['INFO']:
    csq_entry = [s for s in info.split(';') if 'CSQ=' in s]
    if csq_entry:
        csq_content = csq_entry[0].replace('CSQ=', '')
        csq_records = csq_content.split(',') 
        
        for csq_rec in csq_records:
            csq_split = csq_rec.split('|')
            csq_data.append(csq_split)
    else:
        csq_split = [''] * len(base_vcf_columns)
        csq_data.append(csq_split)

csq_df = pd.DataFrame(csq_data, columns=csq_columns)
eqtl_merged_df = pd.concat([eqtl_df.reset_index(drop=True), csq_df.reset_index(drop=True)], axis=1)

eqtl_merged_df['variant_info'] = (eqtl_merged_df["#CHROM"] + '_' + eqtl_merged_df["POS"] + '_' + eqtl_merged_df["REF"] + '_' + eqtl_merged_df["ALT"])
eqtl_set = set(eqtl_merged_df['variant_info'])
mask_eqtl_catalog = (eqtl_catalog['eQTL_Catalog_variant'].isin(eqtl_set))
mask_eqtl_gtex = (eqtl_gtex['eQTL_GTEx_variant_id'].isin(eqtl_set))

eqtl_catalog_filtered_df = eqtl_catalog[mask_eqtl_catalog].reset_index(drop=True)
eqtl_gtex_filtered_df = eqtl_gtex[mask_eqtl_gtex].reset_index(drop=True)

eqtl_catalog_output_file = cfg["output_file"].replace("txt", "eqtl_catalog.txt")
eqtl_gtex_output_file = cfg["output_file"].replace("txt", "eqtl_gtex.txt")
eqtl_catalog_output_file_user = cfg["output_file"].replace("txt", "eqtl_catalog_filtered.txt")
eqtl_gtex_output_file_user = cfg["output_file"].replace("txt", "eqtl_gtex_filtered.txt")

eqtl_catalog_filtered_df.to_csv(eqtl_catalog_output_file, sep="\t", index=False, quoting=3)
eqtl_gtex_filtered_df.to_csv(eqtl_gtex_output_file, sep="\t", index=False, quoting=3)

expanded_reportA['variant_info'] = (expanded_reportA["#CHROM"] + '_' + expanded_reportA["POS"] + '_' + expanded_reportA["REF"] + '_' + expanded_reportA["ALT"])
variant_set = set(expanded_reportA['variant_info'])
eqtl_catalog_user = eqtl_catalog_filtered_df[eqtl_catalog_filtered_df['eQTL_Catalog_variant'].isin(variant_set)].reset_index(drop=True)
etl_gtex_user = eqtl_gtex_filtered_df[eqtl_gtex_filtered_df['eQTL_GTEx_variant_id'].isin(variant_set)].reset_index(drop=True)

eqtl_catalog_user.to_csv(eqtl_catalog_output_file_user, sep="\t", index=False, quoting=3)
etl_gtex_user.to_csv(eqtl_gtex_output_file_user, sep="\t", index=False, quoting=3)


#%% Cell 9 pqtl-Open Targets Genetics
print("Now generating pqtl-OTG results")
csq_data = []
otg_data = [line.split('\t') for line in reporteqtl]
otg_df = pd.DataFrame(otg_data, columns=base_vcf_columns)
for info in otg_df['INFO']:
    csq_entry = [s for s in info.split(';') if 'CSQ=' in s]
    if csq_entry:
        csq_content = csq_entry[0].replace('CSQ=', '')
        csq_records = csq_content.split(',') 
        
        for csq_rec in csq_records:
            csq_split = csq_rec.split('|')
            csq_data.append(csq_split)
    else:
        csq_split = [''] * len(base_vcf_columns)
        csq_data.append(csq_split)

csq_df = pd.DataFrame(csq_data, columns=csq_columns)
otg_merged_df = pd.concat([otg_df.reset_index(drop=True), csq_df.reset_index(drop=True)], axis=1)

otg_merged_df['variant_info'] = (otg_merged_df["#CHROM"] + '_' + otg_merged_df["POS"] + '_' + otg_merged_df["REF"] + '_' + otg_merged_df["ALT"])
otg_set = set(otg_merged_df['variant_info'])
mask = (otg_database['variant_info'].isin(otg_set))

otg_filtered_df = otg_database[mask].reset_index(drop=True)

otg_output_file = cfg["output_file"].replace("txt", "pqtl_otg.txt")
otg_output_file_user = cfg["output_file"].replace("txt", "pqtl_otg_filtered.txt")
otg_filtered_df.to_csv(otg_output_file, sep="\t", index=False, quoting=3)

expanded_reportA['variant_info'] = (expanded_reportA["#CHROM"] + '_' + expanded_reportA["POS"] + '_' + expanded_reportA["REF"] + '_' + expanded_reportA["ALT"])
variant_set = set(expanded_reportA['variant_info'])
otg_filtered_user = otg_filtered_df[otg_filtered_df['variant_info'].isin(variant_set)].reset_index(drop=True)

otg_filtered_user.to_csv(otg_output_file_user, sep="\t", index=False, quoting=3)
