import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import sys

#USER
#	parse_data.py and annotated_volcano.py must be run on every dataset (ALL, CPTAC3, LUAD, LUSC) first
#	run from the home directory
#		python code/PCA.py CPTAC3 2

proj = sys.argv[1]
n_components = int(sys.argv[2])

# read in data (samples x genes)
gene_list = pd.read_csv('data/'+proj+'_TUMOR_combined.csv').iloc[:,0].values

FEMALE_TUMOR = pd.read_csv('data/'+proj+'_TUMOR_fem_collected_expression.csv',header=None).T
FEMALE_TUMOR.set_axis(gene_list, axis=1, inplace=True)
FEMALE_TUMOR['Label'] = 'FEMALE_TUMOR'

FEMALE_NORMAL = pd.read_csv('data/'+proj+'_NORMAL_fem_collected_expression.csv',header=None).T
FEMALE_NORMAL.set_axis(gene_list, axis=1, inplace=True)
FEMALE_NORMAL['Label'] = 'FEMALE_NORMAL'

MALE_TUMOR = pd.read_csv('data/'+proj+'_TUMOR_male_collected_expression.csv',header=None).T
MALE_TUMOR.set_axis(gene_list, axis=1, inplace=True)
MALE_TUMOR['Label'] = 'MALE_TUMOR'

MALE_NORMAL = pd.read_csv('data/'+proj+'_NORMAL_male_collected_expression.csv',header=None).T
MALE_NORMAL.set_axis(gene_list, axis=1, inplace=True)
MALE_NORMAL['Label'] = 'MALE_NORMAL'

# combine data
combined_df = pd.concat([FEMALE_TUMOR,FEMALE_NORMAL,MALE_TUMOR,MALE_NORMAL], ignore_index=True)

# rescale data
x = combined_df.loc[:, gene_list].values
x = StandardScaler().fit_transform(x)
y = combined_df.loc[:,['Label']].values

# PCA
pca = PCA(n_components=n_components)

principalComponents = pca.fit_transform(x)
pcColumns = ['PC'+str(i+1) for i in range(n_components)]

principalDf = pd.DataFrame(data = principalComponents
             , columns = pcColumns)

finalDf = pd.concat([principalDf, combined_df[['Label']]], axis = 1)

# plot
fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('PC1', fontsize = 15)
ax.set_ylabel('PC2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)

targets = ['FEMALE_TUMOR','FEMALE_NORMAL','MALE_TUMOR','MALE_NORMAL']
colors = ['r', 'g', 'b','y']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['Label'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'PC1']
               , finalDf.loc[indicesToKeep, 'PC2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()

plt.savefig('figures/'+proj+'_PCA.png')