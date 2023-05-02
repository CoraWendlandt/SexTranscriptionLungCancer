import pandas as pd
import numpy as np
import sys 

#USER (*make sure to unzip the tar files b/c only the compressed ones fit on git)
# run from home directory
#	python code/deseqs_meta_data.py CPTAC3

extraTitle = str(sys.argv[1])
femNorm = pd.read_csv('data/'+extraTitle+'_NORMAL_fem_counts.csv', index_col=0)
femTum = pd.read_csv('data/'+extraTitle+'_TUMOR_fem_counts.csv',index_col=0)
maleNorm = pd.read_csv('data/'+extraTitle+'_NORMAL_male_counts.csv', index_col=0)
maleTum = pd.read_csv('data/'+extraTitle+'_TUMOR_male_counts.csv', index_col=0)

indexNorm = femNorm.columns.values
indexNorm = np.append(indexNorm, maleNorm.columns.values)
conditionNorm = ['Female']*len(femNorm.columns.values)
conditionNorm = np.append(conditionNorm,['Male']*len(maleNorm.columns.values))
normDf = pd.DataFrame(conditionNorm, index=indexNorm, columns = ['gender'])
normDf.to_csv('data/'+extraTitle+'_normal_meta.csv')
normCounts = pd.concat([femNorm, maleNorm], axis=1)
normCounts.to_csv('data/'+extraTitle+'_NORMAL_deseqs_counts.csv')

indexTum = femTum.columns.values
indexTum = np.append(indexTum, maleTum.columns.values)
conditionTum = ['Female']*len(femTum.columns.values)
conditionTum = np.append(conditionTum,['Male']*len(maleTum.columns.values))
tumDf = pd.DataFrame(conditionTum, index=indexTum, columns = ['gender'])
tumDf.to_csv('data/'+extraTitle+'_tumor_meta.csv')
tumCounts = pd.concat([femTum, maleTum], axis=1)
tumCounts.to_csv('data/'+extraTitle+'_TUMOR_deseqs_counts.csv')
