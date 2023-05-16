import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cross_decomposition import PLSRegression
import matplotlib.pyplot as plt
import sys

#USER
#	parse_data.py and annotated_volcano.py must be run on every dataset (ALL, CPTAC3, LUAD, LUSC) first
#	run from the home directory
#		python code/PLSDA.py LUAD 2

proj = sys.argv[1]
n_components = int(sys.argv[2])

# read in data (samples x genes)
gene_list = pd.read_csv('data/'+proj+'_TUMOR_combined.csv').iloc[:,0].values

FEMALE_TUMOR = pd.read_csv('data/'+proj+'_TUMOR_fem_collected_expression.csv',header=None).T # transposed, columns are genes, rows are patients
FEMALE_TUMOR.set_axis(gene_list, axis=1, inplace=True)
FEMALE_TUMOR['Label'] = -1

MALE_TUMOR = pd.read_csv('data/'+proj+'_TUMOR_male_collected_expression.csv',header=None).T
MALE_TUMOR.set_axis(gene_list, axis=1, inplace=True)
MALE_TUMOR['Label'] = 1

FEMALE_TEST = FEMALE_TUMOR.tail(20)
FEMALE_TEST_y = FEMALE_TEST['Label']
FEMALE_TEST.drop('Label', axis=1)

FEMALE_TRAIN = FEMALE_TUMOR[:-20]
FEMALE_TRAIN_y = FEMALE_TRAIN['Label']
FEMALE_TRAIN.drop('Label', axis=1)


MALE_TEST = MALE_TUMOR.tail(20)
MALE_TEST_y = MALE_TEST['Label']
MALE_TEST.drop('Label', axis=1)

MALE_TRAIN = MALE_TUMOR[:-20]
MALE_TRAIN_y = MALE_TRAIN['Label']
MALE_TRAIN.drop('Label', axis=1)

# combine data
train_x = pd.concat([FEMALE_TRAIN,MALE_TRAIN], ignore_index=True)
train_y = pd.concat([FEMALE_TRAIN_y,MALE_TRAIN_y], ignore_index=True)

test_x = pd.concat([FEMALE_TEST,MALE_TEST], ignore_index=True)
test_y = pd.concat([FEMALE_TEST_y,MALE_TEST_y], ignore_index=True)

train_x = train_x.loc[:, gene_list].values
train_x = StandardScaler().fit_transform(train_x)
test_x = test_x.loc[:, gene_list].values
test_x = StandardScaler().fit_transform(test_x)


def plsda_accuracy(n_comp,train_x,train_y,test_x):

	plsda = None
	plsda = PLSRegression(n_components=n_comp)
	plsda.fit(train_x,train_y)

	Y_pred = 0
	Y_pred = plsda.predict(test_x)
	wrong = 0
	right = 0
	for i in range(len(Y_pred)):
		if Y_pred[i] > 0 and i < 20:
			wrong += 1
		elif Y_pred[i] < 0 and i >= 20:
			wrong +=1
		else:
			right +=1

	return right/(wrong+right)

comps_to_test = [2,3,4,5,6,8,12,24] # for CPTAC3: 2,4,5,12,24 for LUSC: 2,5,6,7,8,12,24 LUAD: 2,3,4,5,6,8,12,24
accuracy = []
for n in comps_to_test:
	temp_acc = plsda_accuracy(n,train_x,train_y,test_x)
	accuracy.append(temp_acc*100)
	print(temp_acc)


def plsda_plot(train_x,train_y,test_x):
	plsda = None
	plsda = PLSRegression(n_components=24)
	plsda.fit(train_x,train_y)

	PCspace = plsda.transform(test_x)

	variance_in_x = np.var(plsda.x_scores_, axis = 0) 
	explained_variance = variance_in_x / np.sum(variance_in_x)

	sx_colors = ['hotpink','dodgerblue'] # colors for each classification for plot
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1) 
	ax.set_xlabel('PC1 (explained variance: '+str(round(explained_variance[0]*100,2))+'%)')
	ax.set_ylabel('PC2 (explained variance: '+str(round(explained_variance[1]*100,2))+'%)')
	ax.set_title('First 2 Components of PLS-DA for '+proj)

	for i in range(PCspace.shape[0]): 
  		if i < 20: # for fem data
  			ax.scatter(PCspace[i,0],PCspace[i,1], c = sx_colors[0]) # plot it pink
  		else: # for male data
  			ax.scatter(PCspace[i,0],PCspace[i,1], c = sx_colors[1]) # plot it blue
	ax.legend(["female","male"]) # set legend
	plt.savefig('figures/'+proj+'_PLSDA_Plot.png')

#plsda_plot(train_x,train_y,test_x)

plt.plot(comps_to_test,accuracy)
plt.xlabel('Number of Components')
plt.ylabel('Prediction Accuracy (%)')
plt.title('PLS-DA Prediction Accuracy for '+proj)
plt.ylim([40,80])
plt.savefig('figures/'+proj+'_PLSDA_Accuracy.pdf')