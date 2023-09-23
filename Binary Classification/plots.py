import dataframe_image as dfi
import pandas as pd
import numpy as np

import matplotlib
from matplotlib import pyplot
import matplotlib.pyplot as plt
import seaborn as sns


feat1 = "ROC curve"
feat2 = "PRE-RECAL curve"

class plots:
    def __init__(self):
        self.clf_name = None
        self.dirp = None
    
    
    def plot_metrics(self,results_df,dirp,filename):
        self.dirp = dirp
        results_styled = results_df.style.background_gradient()
        dfi.export(results_df, f'{self.dirp}{filename}_test_Performance_metrics.png')



    def cv_boxplot(self,acc_df,dirp,filename):
        self.dirp = dirp
        fig = plt.figure()
        fig.suptitle('Algorithm Comparison: 10 fold cross validation')
        ax = fig.add_subplot(1,1,1)
        plt.boxplot(acc_df)
        ax.set_xticklabels(acc_df.columns)
        #plt.show()
        plt.savefig(f"{self.dirp}{filename}_10-fold_cv-accuracy_boxplot.png",bbox_inches='tight',dpi=100)

    def plot_roc_curve(self,roc_table,dirp,filename):
        self.dirp = dirp
        fig = plt.figure(figsize=(8,6))
        import numpy as np
        for i in roc_table.index:
            plt.plot(roc_table.loc[i]['fpr'], 
                    roc_table.loc[i]['tpr'], 
                    label="{}, AUC={:.3f}".format(i, roc_table.loc[i]['auc']))
            
        plt.plot([0,1], [0,1], color='blue', linestyle='--')

        plt.xticks(np.arange(0.0, 1.1, step=0.1))
        plt.xlabel("False Positive Rate", fontsize=15)

        plt.yticks(np.arange(0.0, 1.1, step=0.1))
        plt.ylabel("True Positive Rate", fontsize=15)

        plt.title(feat1, fontweight='bold', fontsize=15)
        plt.legend(prop={'size':11}, loc='lower right')

        #plt.show()

        plt.savefig(f"{self.dirp}{filename}_ROC_curve_test-set_ForPaper.png",bbox_inches='tight',dpi=100)

    
    def plot_pre_recall_curve(self,pre_recal_table,dirp,filename):
        self.dirp = dirp

        fig = plt.figure(figsize=(8,6))

        for i in pre_recal_table.index:
            plt.plot(pre_recal_table.loc[i]['Recall'], 
                    pre_recal_table.loc[i]['Precision'],
                    label="{}, AUC={:.3f}".format(i, pre_recal_table.loc[i]['auc']))


        #plt.plot([0,0], color='blue', linestyle='--')

        plt.xticks(np.arange(0.0, 1.1, step=0.2))
        plt.xlabel("Recall", fontsize=15)

        plt.yticks(np.arange(0.0, 1.1, step=0.2))
        plt.ylabel("Precision", fontsize=15)

        plt.title(feat2, fontweight='bold', fontsize=15)
        plt.legend(prop={'size':9}, loc='lower left')

        #plt.show()

        plt.savefig(f"{self.dirp}{filename}_PRE_curve_test-set_ForPaper.png",bbox_inches='tight',dpi=100)

                
