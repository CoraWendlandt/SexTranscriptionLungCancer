import pandas as pd
import numpy as np

#note* these are genes enriched in males relative to females
#USER
#	parse_data.py and annotated_volcano.py must be run on every dataset (ALL, CPTAC3, LUAD, LUSC) first
#	run from the home directory
#		python code/tumor_enriched_genes.py


# collect upregulated genes
CPTAC3_TUMOR_up = set(pd.read_csv('data/CPTAC3_TUMOR_up.csv')['gene_name'].values)
CPTAC3_NORMAL_up = set(pd.read_csv('data/CPTAC3_NORMAL_up.csv')['gene_name'].values)

LUAD_TUMOR_up = set(pd.read_csv('data/LUAD_TUMOR_up.csv')['gene_name'].values)
LUAD_NORMAL_up = set(pd.read_csv('data/LUAD_NORMAL_up.csv')['gene_name'].values)

LUSC_TUMOR_up = set(pd.read_csv('data/LUSC_TUMOR_up.csv')['gene_name'].values)
LUSC_NORMAL_up = set(pd.read_csv('data/LUSC_NORMAL_up.csv')['gene_name'].values)

ALL_TUMOR_up = set(pd.read_csv('data/ALL_TUMOR_up.csv')['gene_name'].values)
ALL_NORMAL_up = set(pd.read_csv('data/ALL_NORMAL_up.csv')['gene_name'].values)

# get genes only upregulated in tumors
CPTAC3_TUMOR_only = CPTAC3_TUMOR_up-CPTAC3_NORMAL_up
LUAD_TUMOR_only = LUAD_TUMOR_up-LUAD_NORMAL_up
LUSC_TUMOR_only = LUSC_TUMOR_up-LUSC_NORMAL_up
ALL_TUMOR_only = ALL_TUMOR_up-ALL_NORMAL_up

# look at intersection of projects
combined_projects = CPTAC3_TUMOR_only.intersection(LUAD_TUMOR_only,LUSC_TUMOR_only)
all_combined_projects = combined_projects.intersection(ALL_TUMOR_only)

# save results
with open('data/tumor_enriched_genes.txt', 'w') as f:
	f.write('CPTAC3 Tumor Genes: '+ str(CPTAC3_TUMOR_only) + '\n')
	f.write('LUAD Tumor Genes: '+ str(LUAD_TUMOR_only) + '\n')
	f.write('LUSC Tumor Genes: ' + str(LUSC_TUMOR_only) + '\n')
	f.write('ALL Tumor Genes: ' + str(ALL_TUMOR_only) + '\n')
	f.write('Combined Projects Tumor Genes: ' + str(combined_projects) + '\n')
	f.write('ALL + Combined Tumor Genes: ' + str(all_combined_projects) + '\n')
f.close()
