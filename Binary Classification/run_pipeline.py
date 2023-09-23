import pandas as pd
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn import svm
#from sklearn.svm import LinearSVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from binary_classification import binary_classification
import sys

dir_path = sys.argv[1]
file_name = sys.argv[2]
g1 = sys.argv[3]
g2 = sys.argv[4]
fn = sys.argv[5]
OutFile1 = sys.argv[6]
OutFile2 = sys.argv[7]

seed = 0
models = []
models.append(('ADB', AdaBoostClassifier(random_state=seed)))
models.append(('DTREE', DecisionTreeClassifier(random_state=seed)))
models.append(('GB', GradientBoostingClassifier(random_state=seed)))
models.append(('KNN', KNeighborsClassifier(n_neighbors=5)))
models.append(('LDA', LinearDiscriminantAnalysis()))
models.append(('LR', LogisticRegression(random_state=seed)))
models.append(('NB', GaussianNB()))
models.append(('RF', RandomForestClassifier(random_state=seed)))
models.append(('SVM', svm.SVC(probability=True,random_state=seed)))

classifier = binary_classification()
classifier.set_data(dirp = dir_path,
                    fileName = file_name,
                    test_decimal=0.30,cutoff = 0.80, cl1 = g1, cl2 = g2)
all_results={}

for name, model in models:
    classifier.set_classifier(name,model)    
    all_results[name]=classifier.classify()


all_df_res = []
all_roc_res = []
all_pre_res = []
#all_cv_res = []

for v in all_results.values():
    all_df_res.append(v[0])
    all_roc_res.append(v[1])
    all_pre_res.append(v[2])
    #all_cv_res.append(v[3])

all_df_res_dframe = pd.concat(all_df_res,ignore_index=True)
all_roc_res_dframe = pd.concat(all_roc_res,ignore_index=True)
all_pre_res_dframe = pd.concat(all_pre_res,ignore_index=True)
#all_cv_res_dframe = pd.concat(all_cv_res)

acc_df = pd.DataFrame(list(zip(all_results['ADB'][3]['test_accuracy'], all_results['DTREE'][3]['test_accuracy'],
                              all_results['GB'][3]['test_accuracy'], all_results['KNN'][3]['test_accuracy'],
                              all_results['LDA'][3]['test_accuracy'], all_results['LR'][3]['test_accuracy'],
                              all_results['NB'][3]['test_accuracy'], all_results['RF'][3]['test_accuracy'],
                              all_results['SVM'][3]['test_accuracy'])),
               columns =['ADB', 'DTREE','GB','KNN','LDA','LR','NB','RF','SVM'])

from plots import plots

visualization = plots()
visualization.plot_metrics(all_df_res_dframe,
                           dirp = dir_path,
                           filename=fn)
visualization.plot_roc_curve(all_roc_res_dframe,
                             dirp = dir_path,
                             filename=fn)
visualization.plot_pre_recall_curve(all_pre_res_dframe,
                                    dirp = dir_path,
                                    filename=fn)
visualization.cv_boxplot(acc_df,
                         dirp = dir_path,
                         filename=fn)

df = pd.read_excel(rf'{dir_path}{file_name}')
df = df.drop(columns=classifier.to_drop)
df_10 = df.drop(columns=classifier.rf_imp_fea[0])
#df_10.to_excel(r'{dir_path}{OutFile1}',index=False)

out_path = f"{dir_path}{OutFile1}"
with pd.ExcelWriter(out_path) as writer:
     df_10.to_excel(writer,index=False)

df_20 = df.drop(columns=classifier.rf_imp_fea[1])
#df_20.to_excel(r'{dir_path}{OutFile2}',index=False)

out_path2 = f"{dir_path}{OutFile2}"
with pd.ExcelWriter(out_path2) as writer:
     df_20.to_excel(writer,index=False)

