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
path = 'local' 

run_traits = True

# User input filtration parameters
user_af = 0.05
user_ada_score = 0.6
user_rf_score = 0.6
user_revel_score=0.75
user_spliceai_al=0.5
user_spliceai_dg=0.5
user_spliceai_dl=0.5
user_spliceai_ag=0.5
user_bayesdel_addaf_score=0.0692655
user_bayesdel_noaf_score=-0.0570105
user_am_classification=['likely_pathogenic','ambigous']
user_am_pathogenicity=0.564
user_gender= "Male"
user_clinvar = ["Pathogenic","Likely_pathogenic","Uncertain_significance","Benign","Likely_benign","Conflicting_classifications_of_pathogenicity","Conflicting_interpretations_of_pathogenicity"]
user_acmg_classification = ["Pathogenic","Likely_pathogenic","Uncertain_significance","Benign","Likely_benign"]


# config file of setting paths
# shoudld make a new Pharma_df! ________________________________--
config = {
    "local": {
        "fileName": "sample_vep_annotated.vcf.gz",
        "OMIM_inheritance_DBfile": "pheno_OMIM_all.txt",
        "Pharma_dbfile": "clinical_annotation_combined.txt",
        "Trait_dbfile": "Reports_genome_databases_traits_merged_2.txt",
        "clinvar_vcf_path": "clinvar_20250504.vcf.gz",
        "GenCC_path" : "enCCforGene.txt",
        "eqtl_catalog_file" : "eqtl.catalogue.permuted.known.position.txt",
        "eqtl_gtex_file" : "GTEx.all.tissue.known.gene.txt",
        "output_file": "sample.nodup.txt"
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

OMIM_Inheritance_DB = read_db_file(cfg["OMIM_inheritance_DBfile"])
OMIM_Inheritance_DB['phenotypeMimNumber'] = OMIM_Inheritance_DB['phenotypeMimNumber'].astype(str).str.strip().str.replace(r'\.0$', '', regex=True)
OMIM_Inheritance_DB['inheritances'] = OMIM_Inheritance_DB['inheritances'].replace(np.nan,"Inheritance Not provided by OMIM")
geneBaseFile = cfg["geneBaseFile"]
eqtl_catalog = pd.read_csv(cfg["eqtl_catalog_file"], sep="\t", header=0) 
eqtl_catalog.columns = [f"eQTL_Catalog_{col}" for col in eqtl_catalog.columns]
eqtl_gtex = pd.read_csv(cfg["eqtl_gtex_file"], sep="\t", header=0)
eqtl_gtex.columns = [f"eQTL_GTEx_{col}" for col in eqtl_gtex.columns]

# ====== 4) Decision of whether running GWAS/Pharmaco/Trait ======
if run_traits:
    Trait_db = read_db_file(cfg["Trait_dbfile"], drop_allna_cols=True)
    trait_list = Trait_db['variants'].to_list()
    print("Traits matching will run")

#%% ClinVar DB dynamic updates: parse into clinvar_cache
def parse_clinvar_vcf(clinvar_vcf_path):
    print(f"Parsing ClinVar VCF from {clinvar_vcf_path}...")
    clinvar_dict = {}      
    with gzip.open(clinvar_vcf_path, 'rt') as file:
        for line in file:
            if line.startswith("#"):
                continue  
            cols = line.strip().split("\t")
            chrom, pos, id_, ref, alt, qual, filter_, info = cols[:8]
            # get needed clinvar info
            clinvar_info = {}
            for key in ["CLNSIG", "CLNDN", "CLNHGVS", "CLNSIGINCL", "CLNVC", "GENEINFO","CLNSIGCONF", "CLNREVSTAT", "CLNDNINCL"]:
                match = re.search(fr"{key}=([^;]+)", info)
                if match:
                    value = match.group(1)
                    ## note this is very important that all csq single value should separate not by |
                    # if isinstance(value, list):
                    #     value = "&".join(value)  # 如果 ClinVar 解析出来是 list，则连接
                    # **CLNDISDB 里的 `,` 也要替换成 `&`**
                    if key == "CLNDISDB" or key == "CLNDN":
                        value = value.replace(",", "&").replace("|", "&")
                else:
                    value = None
                clinvar_info[key] = value
            # build the key
            key = (chrom, pos, ref, alt)
            clinvar_dict[key] = clinvar_info
    return clinvar_dict

#%%Cell 2 start the analysis program for the vep annotated vcf files, generated report A B C and D || E for the traits
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
file = gzip.open(fileName,'rt')
tLine = file.readline()
i = 0
reportA,reportE = [], []

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
        # directly goes into next line
        tLine = file.readline()
        #print(tLine)
        continue
    #start processing real data rows
    if not headings:
        tLine = file.readline()
        #print(tLine)
        continue
    
    iText = [s for s in iContent[headings.index('INFO')].split(';') if 'CSQ=' in s]
    iText = iText[0].replace('CSQ=','').split(',')
    
    
    #### =================ClinVar version updates ##### =================
    saveFlag1, saveFlag2,saveFlag3,saveFlag4, saveFlag5 = False, False, False, False, False
    
    ## 1.26 2024: In the conference I checked rs6025, which should be SZAvar360812 a famous PGx genes for F5 gene. However it only appears on nodup4files, as a clinVar genes,
    ## not in PGx reports. so I realized this elif here should be all changed to "if", because elif ignores other cases if the first one satisfied the criteria
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
                    #print(parsed_value)
                    #print(jText[col_map['REVEL']])
                            
        # 1) ClinVar and prediction filtered based on users' needs 
        if 'MAX_AF' in col_map:
            max_af_val = jText[col_map['MAX_AF']]
            try:
                max_af_numeric = float(max_af_val) if max_af_val != '' else 0.0
            except (ValueError, TypeError):
                max_af_numeric = 0.0
            jText[col_map['MAX_AF']] = max_af_numeric
        
        if max_af_numeric < user_af:
            
            if any(clinvar_term in jText[col_map['ClinVar_CLNSIG']] for clinvar_term in user_clinvar if 'ClinVar_CLNSIG' in col_map):
                saveFlag1 = "Keep"

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
        
        # 3) judge Trait, need run_traits=True and trait_list
        if run_traits and 'Existing_variation' in col_map:
            ex_variation_val = jText[col_map['Existing_variation']]
            ex_variants = ex_variation_val.split('&')
            if ex_variants and ex_variants[0] in trait_list:
                saveFlag5 = "Trait_var"
    # after for j in range(len(iText)) loop, if saveFlag1/2/3/4/5 has value, append the line to respective report
    if saveFlag1:
        reportA.append(tLine)
    if saveFlag5:
        reportE.append(tLine)

    # print progress every 1000000 lines
    if i % 1000000 == 0:
        print(f"{i} lines processed!")

    # read the next line
    tLine = file.readline()


file.close()
print('Cell 2 VEP annotated File processing done! Now start to map GeneDB and DiseaseDB')
end_time = time.time()
print("Total processing time: {:.2f} seconds".format(end_time - start_time))

base_vcf_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT','Sample_Info']

# Get the column names
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

#%% Cell 4 Scoring System
# ClinVar and ACMG scoring system
expanded_reportA = expanded_reportA[expanded_reportA['acmg_classification'].isin(user_acmg_classification)]

acmg_conditions = [
    expanded_reportA['acmg_classification'] == 'Pathogenic',
    expanded_reportA['acmg_classification'] == 'Likely_pathogenic',
    expanded_reportA['acmg_classification'] == 'Benign',
    expanded_reportA['acmg_classification'] == 'Likely_benign',
    expanded_reportA['acmg_classification'] == 'Uncertain_significance'
]
acmg_choices = [1, 0.8, 0.2, 0.2, 0.5]

clinvar_conditions = [
    expanded_reportA['ClinVar_CLNSIG'].str.contains('Pathogenic', na=False, case=True),
    expanded_reportA['ClinVar_CLNSIG'].str.contains('Likely_pathogenic', na=False, case=True),
    expanded_reportA['ClinVar_CLNSIG'].str.contains('Benign', na=False, case=True), 
    expanded_reportA['ClinVar_CLNSIG'].str.contains('Likely_benign', na=False, case=True),
    expanded_reportA['ClinVar_CLNSIG'].str.contains('Uncertain_significance', na=False, case=True),
    expanded_reportA['ClinVar_CLNSIG'].str.contains('Conflicting_classifications_of_pathogenicity', na=False, case=True),
    expanded_reportA['ClinVar_CLNSIG'].str.contains('Conflicting_interpretations_of_pathogenicity', na=False, case=True)
]
clinvar_choices = [1, 0.8, 0.2, 0.2, 0.5, 0.5, 0.5]

expanded_reportA['ACMG_score'] = np.select(acmg_conditions, acmg_choices, default=0)
expanded_reportA['ClinVar_score'] = np.select(clinvar_conditions, clinvar_choices, default=0)
expanded_reportA['Pathogenicity_score'] = expanded_reportA['ACMG_score'] * expanded_reportA['ClinVar_score']

# AF Score
af_values = pd.to_numeric(expanded_reportA['MAX_AF'], errors='coerce')
min_af = af_values[af_values > 0].min()  
af_filled = af_values.fillna(min_af)
af_filled = af_filled.replace(0, min_af)
af_score_raw = -np.log10(af_filled)
min_score = af_score_raw.min()
max_score = af_score_raw.max()

if max_score > min_score: 
    expanded_reportA['AF_score'] = (af_score_raw - min_score) / (max_score - min_score)
else:
    expanded_reportA['AF_score']

# Prediction Score
spliceai_cols = ['SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL']
for col in spliceai_cols:
    expanded_reportA[col] = pd.to_numeric(expanded_reportA[col], errors='coerce').fillna(0)
expanded_reportA['spliceai_max'] = expanded_reportA[spliceai_cols].max(axis=1)

prediction_conditions = [
    expanded_reportA['Consequence'].str.startswith('splice_', na=False),
    (expanded_reportA['IMPACT'] == 'HIGH') & (~expanded_reportA['Consequence'].str.startswith('splice_', na=False)),
    expanded_reportA['Consequence'].str.startswith('missense_variant', na=False),
    (expanded_reportA['IMPACT'] != 'LOW') & (~expanded_reportA['Consequence'].str.startswith('missense_', na=False))
]

expanded_reportA['am_pathogenicity'] = pd.to_numeric(expanded_reportA['am_pathogenicity'], errors='coerce').fillna(0)
prediction_choices = [expanded_reportA['spliceai_max'], 0.99, 0.66 * expanded_reportA['am_pathogenicity'], 0.33]
expanded_reportA['Prediction_score'] = np.select(prediction_conditions, prediction_choices, default=0)

# eQTL
expanded_reportA['variant_info'] = expanded_reportA[['#CHROM', 'POS', 'REF', 'ALT']].fillna('').astype(str).agg('_'.join, axis=1)

def eqtl_score_fun(data: pd.DataFrame, eqtl_catalog: pd.DataFrame, eqtl_gtex: pd.DataFrame) -> pd.DataFrame:
    
    final_data = data.merge(
        eqtl_catalog, 
        left_on='variant_info', 
        right_on='eQTL_Catalog_variant', 
        how='left'
    ).drop_duplicates()
    
    final_data = final_data.merge(
        eqtl_gtex,
        left_on='variant_info',
        right_on='eQTL_GTEx_variant_id',
        how='left'
    ).drop_duplicates()
    
    return final_data

expanded_reportA = eqtl_score_fun(expanded_reportA, eqtl_catalog, eqtl_gtex)
expanded_reportA = expanded_reportA.drop_duplicates()

expanded_reportA.to_csv(cfg["output_file"], sep='\t', index=False)


#%% Cell 5 Traits
## note: one SNP (rsID) could match to more than one trait phenotypes. such as rs1815739. In the result file we only output one db matching but we do the multiple match in the updated_traits.csv
if run_traits == True:
    from datetime import datetime
    now = datetime.now()
    newLineHeadings = ['Traits name', 'category', 'genes', 'variants', 'Description','Patient genotypes','Genotype Description','Zygosity'];
    newLine = '\t'.join(newLineHeadings)+'\n'
    with open(cfg["output_file"].replace('.txt','_traits.txt'),'w') as f:
        f.write(newLine)
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    print('processing ' + str(len(reportE)) + ' lines of traits' )
    #this irs_id_set is used to record the trait variants this individual have, and later write the one don't have to a "wild type"
    irs_id_set = {s.split('|')[col_map['Existing_variation']].split('&')[0] for s in reportE}
    for i in range(0,len(reportE)):
        current_time = now.strftime("%H:%M:%S")
        ##iText means all transcripts
        iText = [s for s in reportE[i].split('\t')[7].split(';') if 'CSQ=' in s][0].replace('CSQ=','').split(',')
        irs_id = iText[0].split("|")[col_map['Existing_variation']].split('&')[0]
        iREF = reportE[i].split('\t')[3].split(',')
        iALT = reportE[i].split('\t')[4].split(',')
        iGenoTypeList = iREF+iALT
        if "." not in reportE[i].split('\t')[9].split(':')[0]:
            iGenotypInd1 = int(reportE[i].split('\t')[9][0])
            iGenotypInd2 = int(reportE[i].split('\t')[9][2])
            try:
                iGenotype = iGenoTypeList[0]+'/'+ iGenoTypeList[0] + '>' + iGenoTypeList[iGenotypInd1] +'/'+iGenoTypeList[iGenotypInd2]
            except:
                iGenotype = iREF[0]+'>'+iALT[0]+':'+reportE[i].split('\t')[9][0]+reportE[i].split('\t')[9][1]+reportE[i].split('\t')[9][2]
            iTemp = reportE[i].replace('\n','').split('\t')
            iTemp[headings.index('INFO')] = iTemp[headings.index('INFO')].split(';CSQ=')[0]
            if iGenotypInd1 == iGenotypInd2:
                iZygo = 'Homozygous'
            elif iGenotypInd1 == 0 or iGenotypInd2 == 0:
                iZygo = 'Heterozygous'
            else:
                iZygo = 'Compound heterozygous'
            iTraits= Trait_db[Trait_db['variants']== irs_id]['Traits name'].tolist()[0]
            iCategory = Trait_db[Trait_db['variants']== irs_id]['category'].tolist()[0]
            iGenes= Trait_db[Trait_db['variants']== irs_id]['genes'].tolist()[0]
            iDescription = Trait_db[Trait_db['variants']== irs_id]['Description'].tolist()[0].strip()
            iGenoType_description =Trait_db[Trait_db['variants']== irs_id]['Genotype_Description'].tolist()[0].strip()
            newLine_temp ='\t'.join([iTraits,iCategory,iGenes,irs_id,iDescription,iGenotype,iGenoType_description,iZygo])+'\n'
            with open(cfg["output_file"].replace('.txt','_traits.txt'),'a') as f:
                f.write(newLine_temp)
    # Iterate through Trait_db and write entries not in reportE
    for index, row in Trait_db.iterrows():
        if row['variants'] not in irs_id_set:
            print(row['variants'] + "This patient carries wildtype genotype (keep as homo_ref)")
            newLine_temp = '\t'.join([row['Traits name'], row['category'], row['genes'], row['variants'], row['Description'], row['REF']+'/'+row['REF'], row['Genotype_Description'], 'homo_ref']) + '\n'
            with open(cfg["output_file"].replace('.txt', '_traits.txt'), 'a') as f:
                f.write(newLine_temp)
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time) 

# APOE haploytype based on generated results file
def get_apoe_isoform_and_risk(rs429358_genotype, rs7412_genotype):
    """
    based on rs429358_genotype and rs7412_genotype
    return (APOE Isoform, AD risk description) tuple.

    reference table:
    rs429358(Ref:T)  rs7412(Ref:C)   Isoform   AD risk
**  TT        TT       ε2/ε2    Reduced Risk (protective)
?   TT        CT       ε2/ε3    Neutral to Reduced Risk
?   CC        CT       ε2/ε3    Neutral to Reduced Risk
**  TT        CC       ε3/ε3    Neutral (most common)  
?   TC        CC       ε3/ε4    Increased Risk
**  CC        CC       ε4/ε4    Highest Risk
    """
    genotype_map = {
        ("TT", "TT"): ("e2/e2", "Reduced Risk (protective)"),
        ("TT", "CT"): ("e2/e3", "Neutral to Reduced Risk"),
        ("TT", "CC"): ("e3/e3", "Neutral (most common)"),
        ("TC", "CC"): ("e3/e4", "Increased Risk"),
        ("CC", "CC"): ("e4/e4", "Highest Risk"),
    }
    return genotype_map.get((rs429358_genotype, rs7412_genotype), ("Unknown", "Unknown"))

    
rows = []
trait_outfile = cfg["output_file"].replace('.txt','_traits.txt')


with open(trait_outfile, 'r', encoding='utf-8') as f_in:
    lines = f_in.readlines()

rs429358_gt = None
rs7412_gt = None

for line in lines:
    parts = line.strip().split('\t')
    if len(parts) < 6:
        continue  
    
    variant = parts[3]
    genotype_str = parts[5]
    if '>' in genotype_str:
        genotype = genotype_str.split('>')[1]  
    else:
        genotype = genotype_str 
    genotype = genotype.replace('/', '')  
    if variant == 'rs429358':
        rs429358_gt = genotype
        print(parts[5])
        print("rs429358_gt",rs429358_gt)
    elif variant == 'rs7412':
        rs7412_gt = genotype
        print("rs7412",rs7412_gt)

# 3) if both are got, then calculate Isoform & Risk
if rs429358_gt and rs7412_gt:
    print("get APOE genotype from original trait file")
    isoform, risk = get_apoe_isoform_and_risk(rs429358_gt, rs7412_gt)
    newLine_temp_parts = [
        "Alzheimer's disease risk",                    
        "Neurogenic and Cognitive functions",         
        "APOE",                                   
        "rs429358+rs7412",                             
        "Combined APOE genotype from rs429358+rs7412", 
        isoform,                 
        risk,                                     
        ""                                  
    ]
    newLine_temp = '\t'.join(newLine_temp_parts) + '\n'

    with open(trait_outfile, 'a', encoding='utf-8') as f_out:
        f_out.write(newLine_temp)
    print(f"Done. Check: {trait_outfile}")
else:
    isoform, risk = ("Unknown", "Unknown")
    print("did not find two genptypes of rs429358 and rs7412")


