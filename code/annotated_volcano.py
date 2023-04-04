import pandas as pd 
import sys
import matplotlib.pylab as plt
import numpy as np

#USER
# run from home directory
#	python code/annotated_volcano.py data/ALLcombined.csv

# work inspired from: https://hemtools.readthedocs.io/en/latest/content/Bioinformatics_Core_Competencies/Volcanoplot.html

# read in combined expression data
combinedCSV = sys.argv[1]
df = pd.read_csv(combinedCSV)

# scatter everytime as nonspecific
plt.scatter(x=df['Log Fold Change'],y=df['Adjusted P Value'].apply(lambda x:-np.log10(x)),s=1,label="Not significant")

# highlight down- and up- regulated genes
down = df[(df['Log Fold Change']<=-2)&(df['Adjusted P Value']<=0.01)]
up = df[(df['Log Fold Change']>=2)&(df['Adjusted P Value']<=0.01)]

# plot down- and up- regulated genes
plt.scatter(x=down['Log Fold Change'],y=down['Adjusted P Value'].apply(lambda x:-np.log10(x)),s=3,label="Down-regulated",color="blue")
plt.scatter(x=up['Log Fold Change'],y=up['Adjusted P Value'].apply(lambda x:-np.log10(x)),s=3,label="Up-regulated",color="red")

# label plot
plt.xlabel("Log Fold Change")
plt.ylabel("-Log Adjusted p value")
plt.axvline(-2,color="grey",linestyle="--")
plt.axvline(2,color="grey",linestyle="--")
plt.axhline(2,color="grey",linestyle="--")
plt.title("Sex Differences in Lung Cancer")
plt.legend()

plt.savefig('figures/annotated_volcano.png')
