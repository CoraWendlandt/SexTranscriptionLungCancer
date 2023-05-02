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
#	python code/prep_deseqs.py data/ALL_FEMALE_TUMOR/ data/ALL_MALE_TUMOR/ ALL_TUMOR
# OR if running separated data
#	python code/prep_deseqs.py data/CPTAC3_FEMALE_TUMOR/ data/CPTAC3_MALE_TUMOR/ CPTAC3_TUMOR

def get_counts(root, extraTitle):
	# takes in the root directory holding expression data
	rootVals = [] # will collect the sum of expression value to make the average
	geneList = [] # will collect the list of genes recorded
	ttestVals = np.array([]) # collect a matrix (#genes x #samples) recording each expression level
	first = True # first entry needs to set up variables
	i = 0
	cols = []
	for subdir, dirs, files in os.walk(root):
		for file in files:
			if file[-3:] == 'tsv':
				# get current expression file
				curTSV = os.path.join(subdir, file)
				curDf = pd.read_csv(curTSV, delimiter ='\t', skiprows=1)
				if first:
					# set up the variables
					geneList = list(curDf['gene_id'][4:])
					rootVals = list(curDf['unstranded'][4:])
					ttestVals = np.zeros((len(rootVals),1))
					ttestVals[:,0] = np.array(curDf['unstranded'][4:]).T
					first = False
				else:
					# add to the variables
					rootVals = list( map(add, rootVals, list(curDf['unstranded'][4:])) )
					nextCol = np.array(curDf['unstranded'][4:])
					nextCol = np.reshape(nextCol, (len(nextCol),1))
					ttestVals = np.hstack((ttestVals, nextCol))
					tempGeneList = list(curDf['gene_id'][4:])
					if geneList != tempGeneList:
						# ensure genes are the same
						print("DIFFERENT GENES!")
				cols.append(extraTitle+str(i))
				i += 1
	return geneList, cols, ttestVals

# get female and male results
femRoot = sys.argv[1]
maleRoot = sys.argv[2]

extraTitle = str(sys.argv[3])

femGenes, femCols, femttest = get_counts(femRoot, extraTitle+'Fem')
femDf = pd.DataFrame(femttest,columns=femCols,index=femGenes)
femDf.to_csv('data/'+extraTitle+'_fem_counts.csv')
maleGenes, maleCols, malettest = get_counts(maleRoot, extraTitle+'Male')
maleDf = pd.DataFrame(malettest,columns=maleCols,index=maleGenes)
maleDf.to_csv('data/'+extraTitle+'_male_counts.csv')
