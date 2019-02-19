import uproot
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt
import math

print("Ciao Chiara!")

def create_roc_curve(labels, scores, positive_label):
    fpr, tpr, thresholds = metrics.roc_curve(labels, scores, pos_label=positive_label)
    roc_auc = metrics.auc(fpr, tpr)
    
    plt.title('Receiver Operating Characteristic')
    #plt.plot(fpr, tpr, 'b', label='AUC = %0.2f'% roc_auc)
    plt.semilogy(tpr, 1./fpr, 'b', label='AUC = %0.2f'% roc_auc)
    plt.legend(loc='lower right')
    plt.plot([0,1],[0,1],'r--')
    plt.xlim([0.3,1.0])
    plt.ylim([1,100000])
    plt.ylabel('light rejection')
    plt.xlabel('b efficiency')
    plt.grid(True)
    plt.show()


f = uproot.open("/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Initial/file.root")
t = f["bTag_AntiKt4EMTopoJets"]

mv1 = t.array("jet_mv1").flatten()
flav = t.array("jet_truthflav").flatten()
eta = t.array("jet_eta").flatten()
pt = t.array("jet_pt").flatten()

index_list = []
for ij in range(len(pt)):
    if pt[ij]>20 and math.fabs(eta[ij]) > 2.4 and math.fabs(eta[ij]) < 4.0 and (not math.fabs(flav[ij])==4): 
        index_list.append(ij)
        
create_roc_curve(np.take(flav, index_list), np.take(mv1, index_list), 5)
