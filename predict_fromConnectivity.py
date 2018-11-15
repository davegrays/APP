#!/Users/dgrayson1/anaconda/bin/python2.7

import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
from sklearn import linear_model, preprocessing
from glmnet import LogitNet
from sklearn.model_selection import train_test_split, cross_val_score

#load brain connectivity data
msubs = pd.read_table('test/SUBS_164_123maleids.txt',header=None,names=['connSubjects']) #ordered list of subject IDs in mat file
mconns = pd.read_csv('test/connection_names_scale33.csv',header=None) #ordered list of connections in mat file
mat = loadmat('test/123connectome_scale33_COM_ranknormal.mat') # load mat-file
mdata = mat['subject_array_2D'] # variable in mat file

#load large demographic/behavior dataframe
df=pd.read_table('test/subject_info_164.txt')

#intersect connectivity subjects with larger dataframe
df2 = msubs.merge(df,left_on='connSubjects', right_on='SubjectID', how='inner')

#define function to regress out nuisance covariates from a feature
def regress(outcome, predictors):
    model = sm.OLS(outcome, predictors).fit()
    return model.resid

#regress age and TCV out from each connection
ageandTCV = sm.add_constant(df2.loc[:,['age', 'TCV']])
for col in range(0, len(mdata[1])):
    mdata[:,col] = regress(mdata[:,col],ageandTCV)

#prep variables for logistic regression
y = df2.diagnosis.copy()
y[y < 0] = 0 #ensure diagnosis is either 0 or 1, not -1 or 1
#X = preprocessing.StandardScaler().fit_transform(mdata) #scale features
X = preprocessing.MinMaxScaler().fit_transform(mdata) #scale features

#regularized logistic regression with cross-validation
#try to find optimal or semi-optimal lambdas (inverse of C in sklearn)
numcrossvals = 10
Rpar = np.array((50, 10, 5, 1, 0.5, 0.1, 0.05))
regpenalty = "l1"

def plot_logr_cvmetrics(X,y,numcrossvals,Rpar,regpenalty):
    scores=np.zeros((len(Rpar),numcrossvals))
    numcoeffs=np.zeros(len(Rpar))
    for i, C in enumerate(Rpar):
        clf = linear_model.LogisticRegression(C=C,class_weight="balanced",penalty=regpenalty)
        scores[i,:] = cross_val_score(clf, X, y, cv=numcrossvals)
        logm = clf.fit(X, y)
        numcoeffs[i] = sum(logm.coef_.ravel()>0)
    for i in range(0, numcrossvals):
        plt.plot(np.log(1/Rpar),scores[:,i])
    plt.title('Accuracy vs ' + regpenalty + ' regularization')
    plt.ylabel('Accuracy');plt.xlabel('Lambda')
    plt.show()
    plt.plot(np.log(1/Rpar),numcoeffs)
    plt.title('Coefficients kept vs ' + regpenalty + ' regularization')
    plt.ylabel('#coeffs');plt.xlabel('Lambda')
    plt.show()

plot_logr_cvmetrics(X,y,numcrossvals,Rpar,regpenalty)

#building the model with all the data for future predictions
logm = linear_model.LogisticRegression(C=1,class_weight="balanced",penalty="l1")
logm = logm.fit(X, y)
#predict labels
p = logm.predict(X)
#or probability estimates
p = logm.predict_proba(X)
#show the connections implicated
print mconns.T[logm.coef_.ravel()>0]







#See also
#sklearn.linear_model.SGDClassifier
#sklearn.svm.LinearSVC

#glmnet stuff
#logm = LogitNet(alpha=0.75, n_splits=3)
#logm = logm.fit(X, y)
# predict labels
#p = logm.predict(X)
# or probability estimates
#p = logm.predict_proba(X)