#### Load the necessary modules###

import pandas as pd
import numpy as np
import pickle

import seaborn as sns

sns.set(style='whitegrid', rc={'axes.grid': False})

from itertools import cycle

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler




from sklearn import model_selection
from sklearn.model_selection import train_test_split
from collections import Counter




from sklearn.metrics import roc_curve, auc, precision_recall_curve, f1_score, roc_auc_score
 
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import confusion_matrix


from sklearn.pipeline import make_pipeline

from collections import Counter
from imblearn.over_sampling import SMOTE



import warnings
warnings.filterwarnings("ignore")

class binary_classification:
    
    def __init__(self):
        self.clf_name = None
        self.clf_model = None
        self.X_train = None
        self.y_train = None
        self.X_test = None
        self.y_test = None
        self.to_drop = None
        self.dirp = None
        self.rf_imp_fea = [None] * 2
        #self.imp.ft.df = [None] * 2
        #self.df = None

    def set_classifier(self,clf_name,clf_model):
        self.clf_name=clf_name
        self.clf_model=clf_model
        
    def set_data(self,dirp,fileName,test_decimal,cutoff,cl1,cl2):
        self.dirp = dirp
        df = pd.read_excel(rf'{self.dirp}{fileName}')
        df=df.drop(['friendlyName'], axis=1)
        #df=df.drop(['Unnamed: 0'], axis=1)
        df['cluster']=df['cluster'].map({cl1:0,cl2:1})
        X = df[df.columns[:-1]]
        y = df.cluster
        X_train, self.X_test, y_train, self.y_test = train_test_split(X, y, test_size = test_decimal, shuffle=True,random_state = 0)
        
        if(len(y_train) > 20.0):
            corr_matrix = X_train.corr().abs()
            upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(np.bool))
            self.to_drop = [column for column in upper.columns if any(upper[column] > cutoff)]
            print(len(self.to_drop))
            df = df.drop(columns=self.to_drop)
            X = df[df.columns[:-1]]
            y = df.cluster 
            X_train,self.X_test,y_train, self.y_test=train_test_split(X,y,test_size=test_decimal, shuffle=True, random_state=0)
        
        counter = Counter(y_train)
        print('Before',counter)

        ###Random oversample with SMOTE
        oversample = SMOTE()
        X_train_sm, y_train_sm = oversample.fit_resample(X_train, y_train)

        counter = Counter(y_train_sm)
        print('After',counter)

        self.X_train = X_train_sm
        self.y_train = y_train_sm
        print(len(self.y_train))

    def classify(self): 
        self.clf_model.fit(self.X_train, self.y_train)

        with open(rf'{self.dirp}{self.clf_name}.sav','wb') as f:
             pickle.dump(self.clf_model,f)
        #  
        pred_probs = self.clf_model.predict_proba(self.X_test)
        preds = self.clf_model.predict(self.X_test)
        
        #print(pred_probs)
        #print(preds)

        if (self.clf_name == 'RF' and (len(self.X_test.columns) >= 10)):
            importances = self.clf_model.feature_importances_
            sorted_indices = np.argsort(importances)[::-1]
            not_imp_features_10 = self.X_test.columns[sorted_indices][10:len(self.X_test.columns)]
            not_imp_features_20 = self.X_test.columns[sorted_indices][20:len(self.X_test.columns)]
            self.rf_imp_fea[0] = not_imp_features_10
            self.rf_imp_fea[1] = not_imp_features_20
            
        
        
        # Confusion Matrix
        TN,FP,FN,TP = confusion_matrix(self.y_test, preds).ravel()
        
        
        print(confusion_matrix(self.y_test, preds))
        
        # Cohen Kappa score
        c_score = round(cohen_kappa_score(self.y_test,preds),2)
        
        # Recall / Sensitivity
        SN = round(TP/(TP+FN),2)
        # Specificity
        SP = round(TN/(TN+FP),2)
        # Precision
        PRE = round(TP/(TP+FP),2)
        # Accuracy
        ACC = round((TP+TN)/(TP+TN+FP+FN),2)
        # F1 score
        f1 = round(2*(PRE * SN)/(PRE + SN),2)
        

        # ROC curve
        fpr, tpr, thresh = roc_curve(self.y_test, pred_probs[:,1], pos_label=1)
        auc_score = roc_auc_score(self.y_test, pred_probs[:,1])
        
        
        # Precision recall curve
        P, R, _ = precision_recall_curve(self.y_test, pred_probs[:,1], pos_label=1)
        auc_pre = auc(R, P)
        
        from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score
        from sklearn.model_selection import RepeatedKFold

        scoring = {'accuracy' : make_scorer(accuracy_score), 
            'sensitivity' : make_scorer(recall_score),
            'specificity': make_scorer(recall_score,pos_label=0),
            'precision' : make_scorer(precision_score),
            'f1_score' : make_scorer(f1_score)} 
        
        # K fold cross validation
        kfold = model_selection.KFold(n_splits=10, shuffle=True,random_state=0)
        
        pipeline = make_pipeline(StandardScaler(), self.clf_model)
        
        with open(rf'{self.dirp}{self.clf_name}_cv.sav','wb') as f:
             pickle.dump(pipeline,f)
        
        
        
        #cv_results = model_selection.cross_val_score(model, X_train, y_train, cv=kfold, scoring="accuracy")
        cv_results = model_selection.cross_validate(pipeline, self.X_train, self.y_train, cv=kfold, scoring=scoring)
        cv_results = pd.DataFrame(cv_results)
        cv_results = cv_results.round(decimals = 2)
        cv_results.loc['mean'] = cv_results.mean()


        df_res = pd.DataFrame({'Classifier':[self.clf_name], 'Sensitivity' : [SN], "Specificity" : [SP], "Cohen_Kappa_Score":[c_score], "Precision": [PRE], "f-1 score":[f1], "Accuracy" : [ACC]})
        result_table = pd.DataFrame({"Classifier":[self.clf_name], 'fpr':[fpr], 'tpr':[tpr],'auc':[auc_score]})
        table_pre_recal = pd.DataFrame({"Classifier":[self.clf_name], 'Precision':[P], 'Recall':[R],'auc':[auc_pre]})
        return [df_res,result_table,table_pre_recal,cv_results]
    
    
        
        
