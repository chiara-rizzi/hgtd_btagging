import uproot
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
import argparse

print("Ciao Chiara!")

parser = argparse.ArgumentParser(description='Make ROC curves.')
parser.add_argument('--tagger', type=str, default="MV1")
parser.add_argument('--largeEta', type=int, default=1)

args = parser.parse_args()

tagger=args.tagger
largeEta=args.largeEta>0

def getVariable( tagger ):
  if tagger=="MV1":    return "jet_mv1"
  if tagger=="MV1c":   return "jet_mv1c"
  if tagger=="MV2c00": return "jet_mv2c00"
  if tagger=="MV2c10": return "jet_mv2c10"
  if tagger=="MV2c20": return "jet_mv2c20"
  if tagger=="MVb":    return "jet_mvb"
  if tagger=="IP3D":   return  "jet_ip3d_llr"
  if tagger=="IP3D+SV1":  return  "jet_sv1ip3d"
  if tagger=="SV1":       return "jet_sv1_m" 
  if tagger=="JetFitter": return  "jet_jf_m"
  return "0"


keys = ["ITK", "Initial", "Int1", "Int2", "Final"]
#keys = ["ITK","Int1"]

condition = {
    "ITK":{
        "file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/ITK/file.root",
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/output/Intermediate_post_Replacement/user.crizzi.mc15_14TeV.117050.PowhegPythia_s3348_s3347_r10900_r11003.btag_HGTDemul_Intermediate_post_Replacement_v1_Akt4EMTo/user.crizzi.17121712.Akt4EMTo._001967.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": "ITK-only",
        "color": "k"
        },
    "Int1":{
        "file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Intermediate_pre_Replacement_trkEffselfTag/file.root",
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/output/Intermediate_post_Replacement/user.crizzi.mc15_14TeV.117050.PowhegPythia_s3348_s3347_r10900_r11003.btag_HGTDemul_Intermediate_post_Replacement_v1_Akt4EMTo/user.crizzi.17121712.Akt4EMTo._001966.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": "Int. pre-repl",
        "color": "y"
        },
    "Int2":{
        "file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Intermediate_post_Replacement_trkEffselfTag/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": "Int. post-repl",
        "color": "b"
        },
    "Initial":{
        "file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Initial_trkEffselfTag/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": "Initial",
        "color": "g"
        },
    "Final":{
        "file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Final_trkEffselfTag/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": "Final",
        "color": "r"
        }
    }


def input_arrays(t):
    mv1 = t.array(getVariable(tagger)).flatten()
    flav = t.array("jet_LabDr_HadF").flatten()
    eta = t.array("jet_eta").flatten()
    pt = t.array("jet_pt").flatten()
    dRminToB = t.array("jet_dRminToB").flatten()
    dRminToC = t.array("jet_dRminToC").flatten()
    dRminToT = t.array("jet_dRminToT").flatten()
    isPU = t.array("jet_isPU").flatten()
    truthMatch = t.array("jet_truthMatch").flatten()
    PVz_appo = t.array("PVz").flatten()
    truth_PVz_appo = t.array("truth_PVz").flatten()

    truth_PVz=[]
    PVz=[]
    mv1_appo = t.array("jet_pt")
    for i in range(len(mv1_appo)):
        for j in range(len(mv1_appo[i])):
            PVz.append(PVz_appo[i])
            truth_PVz.append(truth_PVz_appo[i])
        
    index_list = []

    for ij in range(len(pt)):
        #if pt[ij]>20 and math.fabs(eta[ij]) > 2.4 and math.fabs(eta[ij]) < 4.0 and (not math.fabs(flav[ij])==4): 
        if pt[ij]>20  and math.fabs(eta[ij]) < 3.6 and truthMatch[ij]==1 and isPU[ij]==0 and math.fabs(PVz[ij] - truth_PVz[ij])<0.1:
            if math.fabs(flav[ij]) == 5 or ((not math.fabs(flav[ij])==4) and (not math.fabs(flav[ij])==15) and dRminToB[ij]>0.8 and dRminToC[ij]>0.8 and dRminToT[ij]>0.8):        
                if (not largeEta) or math.fabs(eta[ij]) > 2.4:
                    index_list.append(ij)
    return np.take(flav, index_list), np.take(mv1, index_list)

#fig = plt.figure(figsize=(5,6))
fig = plt.figure()
gs = gridspec.GridSpec(1, 1,
                       height_ratios=[1],
                       hspace=0.05)
ax1 = plt.subplot(gs[0])
plt.title(tagger+' ROC curve') #
fpr_all={}
tpr_all={}
for key in keys:
    f = uproot.open(condition[key]["file"])
    t = f[condition[key]["tree_name"]]
    labels, scores = input_arrays(t)
    fpr, tpr, thresholds = metrics.roc_curve(labels, scores, 5) 
    fpr_all[key]=fpr
    tpr_all[key]=tpr
    roc_auc = metrics.auc(fpr, tpr)
    print(condition[key]['label']) 
    print(roc_auc)
    plt.semilogy(tpr, 1./fpr, label=condition[key]['label']+' AUC = %0.4f'% roc_auc, color=condition[key]['color'])
    
plt.legend(loc='best') 
plt.xlim([0.3,1.0])
plt.ylim([1,10000])
plt.ylabel('light-jet rejection')
plt.xlabel('b-jet efficiency')
plt.text(0.05, 0.25, 'ATLAS Simulation Internal', size='large',transform=ax1.transAxes)
plt.text(0.05, 0.15, '$t\overline{t}$ simulation', size='medium',transform=ax1.transAxes)
if largeEta:
    plt.text(0.05, 0.10, 'jet p$_{T}>$ 20 GeV, $|\eta|>$2.4', size='medium',transform=ax1.transAxes)
else:
    plt.text(0.05, 0.10, 'jet p$_{T}>$ 20 GeV', size='medium',transform=ax1.transAxes)
#ax.show()
#plt.gca().axes.get_xaxis().set_visible(False)
#ax1.set_xticklabels( () )
plt.grid(True)
#ax2 = plt.subplot(gs[1])

#plt.xlabel('b efficiency')
#plt.ylabel('Ratio to ITK')
#for key in keys:
#    plt.semilogy(tpr_all[key], fpr_all["ITK"]/fpr_all[key], color=condition[key]['color'])
#    print(tpr_all[key])


if largeEta:
    plt.savefig("roc_"+getVariable(tagger)+"_largeEta_trkEffselfTag.pdf")   
else:
    plt.savefig("roc_"+getVariable(tagger)+".pdf")   

