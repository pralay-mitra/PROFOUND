#!/usr/bin/env python
# coding: utf-8

#The demo_code.py file is a dummy python file illustrating how the PU learning based classifiers were constructed and how the reported results 
#were generated in PROFOUND. Please keep in mind that it is a dummy code and is only provided for information.

import os
import numpy as np
from sklearn.metrics import recall_score
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier as rfc


# class that will assign weight to each instance.
class PUClassifier(object):
    def __init__(self, trad_clf=None, n_folds=2):
        self.trad_clf = trad_clf
        self.n_folds = n_folds

    def fit(self, X, s):    #Function that calculates value of c.
        if self.trad_clf is None:
            self.trad_clf=rfc(n_estimators=1500, class_weight="balanced", n_jobs=4)
        c = np.zeros(self.n_folds)
        skf=StratifiedKFold(n_splits=self.n_folds, shuffle=True)
        for i, (itr, ite) in enumerate(skf.split(X,s)):
            self.trad_clf.fit(X[itr], s[itr])
            c[i] = self.trad_clf.predict_proba(X[ite][s[ite]==1])[:,1].mean()
        self.c = c.mean()
        return self

    def sample(self, X, s):     #Function that assigns weight to each instance.
        if not hasattr(self, "c"):
            self.fit(X, s)
        X_positive = X[s==1]
        X_unlabeled = X[s==0]
        n_positive = X_positive.shape[0]
        n_unlabeled = X_unlabeled.shape[0]
        X_train = np.r_[X_positive, X_unlabeled, X_unlabeled]
        y_train = np.concatenate([np.repeat(1, n_positive), np.repeat(1, n_unlabeled), np.repeat(0, n_unlabeled)])
        self.trad_clf.fit(X, s)
        coun1=coun2=0
        p_unlabeled = self.trad_clf.predict_proba(X_unlabeled)[:,1]
        w_positive = ((1 - self.c) / self.c) * (p_unlabeled / (1 - p_unlabeled))
        for n in w_positive:
            coun1=coun1+n
        coun2=n_unlabeled-coun1
        coun1=coun1+n_positive
        zero_rat=coun1/(coun2+coun1)
        one_rat=1-zero_rat
        w_negative = 1 - w_positive
        sample_weight = np.concatenate([np.repeat(1.0, n_positive), w_positive, w_negative])
        return X_train, y_train, sample_weight, one_rat, zero_rat
 


os.system("mkdir result") #creating result folder


k=10     #Number of fold

f1=open("known_positive_sample.txt","r")    #File that contain feature set of known positive instance.
f2=open("unknown_sample.txt","r")       #File that contain feature set of unknown positive instance.
f3=open("result/"+str(k)+"fold_positive.txt","w")	#File that contain number of times each instances is predicated as label '1'
f5=open("result/"+str(k)+"fold_result.txt","w")	#File that contain recall and fall out rate for each prediction.
f6=open("result/"+str(k)+"fold_feature_importance.txt","w")	#File that contain contribution of each feature in each prediction.


lines1=f1.readlines()
lines2=f2.readlines()
size1=len(lines1)
size2=len(lines2)
fsize=len(lines1[0].split("\t"))

data_P=[]
data_U=[]

for line in lines1:
    line=line[:-1]
    line.strip()
    l=line.split("\t")
    l = list(map(float, l))
    data_P+=l

for line in lines2:
    line=line[:-1]
    line.strip()
    l=line.split("\t")
    l = list(map(float, l))
    data_U+=l

data_P=np.array(data_P)
data_P=data_P.reshape(size1,fsize)

data_U=np.array(data_U)
data_U=data_U.reshape(size2,fsize)

data_P_label=np.ones(len(data_P))
data_U_label=np.zeros(len(data_P))

data_T = np.concatenate((data_P, data_U), axis=0)
data_T_label= np.concatenate([np.repeat(1, len(data_P)), np.repeat(0, len(data_U))])



res=np.zeros(len(data_T_label))
result_all_list=[]

for itra in range(10):

    skf=StratifiedKFold(n_splits=k, shuffle=True)

    final_array=[]
    result=[]
    result1=[]
    recall_score_mat=[]    
    fall_out_mat=[]
    result_array_list=[]

    f5.write(str(itra+1)+"\n")

    for itr, ite in skf.split(data_T, data_T_label):
            data_T_xtr=data_T[itr]
            data_T_label_ytr=data_T_label[itr]
            data_T_xtt=data_T[ite]
            data_T_label_ytt=data_T_label[ite]

            n = len(data_T_xtr)
            X = data_T_xtr
            y = data_T_label_ytr
            s = y.copy()
            pu = PUClassifier(n_folds=5)

            X_train, y_train, sample_weight, one_rat, zero_rat = pu.sample(X, s)    #Assigning weight to each instance. 
 
            ypr1=np.zeros(len(data_T_xtt))

            clf = rfc(n_estimators=1500,criterion='gini',max_features=0.7, class_weight={0:zero_rat, 1:one_rat}, n_jobs=4).fit(X_train, y_train, sample_weight=sample_weight)       #Training random forest model that predict the probability of native folding after deletion
            ypr = clf.predict(data_T_xtt)       #Predicting the probability of folding after deletion
            
            rec_score=recall_score(data_T_label_ytt, ypr) # calculating recall score
            importances=clf.feature_importances_
            out = confusion_matrix(data_T_label_ytt, ypr)
            fall_out_rate=float(out[0][1])/float(out[0][0]+out[0][1])  # calculating fall out rate
            ypr1=ypr
            importances_list = list(importances)
            recall_score_mat.append(rec_score)
            fall_out_mat.append(fall_out_rate)


            result_array = [out[0][0], out[0][1], out[1][0], out[1][1], rec_score, fall_out_rate]
            print (result_array)

            result_array_list.append(result_array)
            result_all_list.append(result_array)

            string6 = "\t".join(str(a) for a in importances_list)
            f6.write(string6+"\n")


            string5 = "\t".join(str(a) for a in result_array)
            f5.write(string5+"\n")

            index=0
            for count in ite:						
                res[count]=res[count]+ypr1[index]
                index=index+1

    f5.write("\n\n")

    ncols = len(result_array_list[0])
    nrows = len(result_array_list)
    mean_results = ncols*[0] 
    for col in range(ncols):
        for row in range(nrows):
            mean_results[col] += result_array_list[row][col]
    nelem = float(nrows)
    mean_results = [val/nelem for val in mean_results]
    string5 = "\t".join(str(val) for val in mean_results)
    f5.write(string5+"\n\n\n")



for i in res:
    f3.write(str(i)+"\n")
f3.close()




ncols = len(result_all_list[0])
nrows = len(result_all_list)
mean_results = ncols*[0] 
for col in range(ncols):
    for row in range(nrows):
        mean_results[col] += result_all_list[row][col]
nelem = float(nrows)
mean_results = [val/nelem for val in mean_results]
string5 = "\t".join(str(val) for val in mean_results)
f5.write(string5+"\n\n\n")

f5.close()
f6.close()