import os
import sys
import pandas as pd
import numpy as np
from scipy import stats
from operator import add
import matplotlib.pyplot as plt
import math

#USER (*make sure to unzip the tar files b/c only the compressed ones fit on git)
# run from home directory
#	python code/parse_data.py data/ALL_FEMALE_TUMOR/ data/ALL_MALE_TUMOR/ ALL
# OR if running separated data
#	python code/parse_data.py data/CPTAC3_FEMALE_TUMOR/ data/CPTAC3_MALE_TUMOR/ CPTAC3
# expect about 6ish min to run

def get_avg_expression(root):
	# takes in the root directory holding expression data
	rootVals = [] # will collect the sum of expression value to make the average
	geneList = [] # will collect the list of genes recorded
	ttestVals = np.array([]) # collect a matrix (#genes x #samples) recording each expression level
	first = True # first entry needs to set up variables

	for subdir, dirs, files in os.walk(root):
		for file in files:
			if file[-3:] == 'tsv':
				# get current expression file
				curTSV = os.path.join(subdir, file)
				curDf = pd.read_csv(curTSV, delimiter ='\t', skiprows=1)
				if first:
					# set up the variables
					geneList = list(curDf['gene_name'][4:])
					rootVals = list(curDf['tpm_unstranded'][4:])
					ttestVals = np.zeros((len(rootVals),1))
					ttestVals[:,0] = np.array(curDf['tpm_unstranded'][4:]).T
					first = False
				else:
					# add to the variables
					rootVals = list( map(add, rootVals, list(curDf['tpm_unstranded'][4:])) )
					nextCol = np.array(curDf['tpm_unstranded'][4:])
					nextCol = np.reshape(nextCol, (len(nextCol),1))
					ttestVals = np.hstack((ttestVals, nextCol))
					tempGeneList = list(curDf['gene_name'][4:])
					if geneList != tempGeneList:
						# ensure genes are the same
						print("DIFFERENT GENES!")
	rootVals[:] = list(map(lambda x: x/len(rootVals), rootVals)) #get average
	return geneList, rootVals, ttestVals

# get female and male results
femRoot = sys.argv[1]
maleRoot = sys.argv[2]

extraTitle = ''
if len(sys.argv) == 4:
	extraTitle = str(sys.argv[3])

femGenes, femVals, femttest = get_avg_expression(femRoot)
np.savetxt('data/'+extraTitle+'fem_collected_expression.csv', femttest,  delimiter=',')
maleGenes, maleVals, malettest = get_avg_expression(maleRoot)
np.savetxt('data/'+extraTitle+'male_collected_expression.csv', malettest,  delimiter=',')

if femGenes != maleGenes:
	# make sure they have the same genes
	print("GENDER DIFFERENT GENES!")

logFoldChanges = []
def fold_change(i, j):	
	# takes in two values and calculates the log fold change
	if j != 0:
		val = i / j
		if val > 0: 
			return math.log2(val)
		elif val < 0:
			return -1*math.log2(-1*val)
		else:
			return np.nan
	else:
		return 0

# get log full change of variables
logFoldChanges = list(map(lambda i,j: fold_change(i,j), maleVals, femVals))

tvals = []
pvals = []
for i in range(len(femVals)):
	# calculate t statistic and p value for changes in genes from female to male
	curttestMale = malettest[i,:]
	curttestFem = femttest[i,:]
	tval, p = stats.ttest_ind(curttestMale, curttestFem)
	tvals.append(tval)
	pvals.append(p)

# multiple hypothesis correction on p values
adjPvals = [x * len(pvals) for x in pvals]

# collect and save all the data
combinedDf = pd.DataFrame(list(zip(femVals, maleVals, logFoldChanges, tvals, pvals, adjPvals)),
               index = femGenes, columns =['Female', 'Male', 'Log Fold Change', 't statistic', 'P Value', 'Adjusted P Value'])
combinedDf.to_csv('data/'+extraTitle+'combined.csv')

# volcano plot of changes from female to male
plt.scatter(x=combinedDf['Log Fold Change'],y=combinedDf['Adjusted P Value'].apply(lambda x:-np.log10(x)),s=1)
plt.xlabel("Log Fold Changes")
plt.ylabel("-Log Adjusted p value")
plt.title('Sex Differences in Lung Cancer'+extraTitle)
plt.savefig('figures/'+extraTitle+'SexLungVolcano.png')
