import uproot
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
import argparse
from scipy import interpolate
import collections

print("Ciao Chiara!")

parser = argparse.ArgumentParser(description='Make ROC curves.')
parser.add_argument('--tagger', type=str, default="MV1")
parser.add_argument('--largeEta', type=int, default=1)
parser.add_argument('--twoTrk', type=int, default=1)

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

def getVariableNtrk( tagger ):
  if tagger=="MV1":    return "jet_ip3d_ntrk"
  if tagger=="MV1c":   return "jet_ip3d_ntrk"
  if tagger=="MV2c00": return "jet_ip3d_ntrk"
  if tagger=="MV2c10": return "jet_ip3d_ntrk"
  if tagger=="MV2c20": return "jet_ip3d_ntrk"
  if tagger=="MVb":    return "jet_ip3d_ntrk"
  if tagger=="IP3D":   return  "jet_ip3d_ntrk"
  if tagger=="IP3D+SV1":  return  "jet_sv1_ntrk"
  if tagger=="SV1":       return "jet_sv1_ntrk" 
  if tagger=="JetFitter": return  "jet_ip3d_ntrk"
  return "0"


keys = ["ITK", "Initial", "Int1", "Int2", "Final"]
#keys = ["ITK","Initial"]

condition = {
    "ITK":{
    #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/ITK/file.root",
         "file":"../input/ITK/file.root",
         "tree_name":"bTag_AntiKt4EMTopoJets",
         "label": "ITK-only",
         "color": "k"
    },
    "Int1":{
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Intermediate_pre_Replacement_trkEffselfTag/file.root",
        "file":"../input/Intermediate_pre_Replacement_trkEffselfTag/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": "Int. pre-repl",
        "color": "y"
        },
    "Int2":{
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Intermediate_post_Replacement_trkEffselfTag/file.root",
        "file":"../input/Intermediate_post_Replacement_trkEffselfTag/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": "Int. post-repl",
        "color": "b"
        },
    "Initial":{
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Initial_trkEffselfTag/file.root",
        "file":"../input/Initial_trkEffselfTag/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": "Initial",
        "color": "g"
        },
    "Final":{
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Final_trkEffselfTag/file.root",
        "file":"../input/Final_trkEffselfTag/file.root",
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
#fig = plt.figure()
gs = gridspec.GridSpec(2, 1,
                       height_ratios=[3,1],
                       hspace=0.05)
ax1 = plt.subplot(gs[0])
plt.title(tagger+' ROC curve') #
fpr_all={}
tpr_all={}
for i,key in enumerate(keys):
    print(key)
    f = uproot.open(condition[key]["file"])
    t = f[condition[key]["tree_name"]]
    labels, scores = input_arrays(t)
    fpr, tpr, thresholds = metrics.roc_curve(labels, scores, 5) 
    x = tpr
    y = 1/fpr
    points = list(zip(x, y))
    #print("   Removing duplicates...")
    x1 = []
    y1 = []
    for ip,p in enumerate(points):
      if p[0] > 0.29999:
        if p[0] > points[ip-1][0]:
          x1.append(p[0])
          y1.append(p[1])
    x = np.array(x1)
    y = np.array(y1)
    if i==0:
      print("   setting reference_beff")
      xnew = x      
    #print("   x: ",x[0:20])
    #print("   y: ",y[0:20])
    #f2 = interpolate.interp1d(x, y, kind='cubic',fill_value="extrapolate")
    #ynew=f2(xnew)
    tck = interpolate.splrep(x, y)
    ynew = interpolate.splev(xnew, tck, der=0)

    fpr_all[key]=ynew
    tpr_all[key]=xnew
    roc_auc = metrics.auc(fpr, tpr)
    print(condition[key]['label']) 
    print(roc_auc)
    plt.semilogy(xnew, ynew, label=condition[key]['label']+' AUC = %0.4f'% roc_auc, color=condition[key]['color'])
    
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
plt.grid(True)

ax2 = plt.subplot(gs[1])
for i,key in enumerate(keys):
    print("\n")
    print(key)
    y=fpr_all[key] # light rejection
    x=tpr_all[key]
    #x=np.nan_to_num(x)
    #y=np.nan_to_num(y)
    #print("   x: ",x[0:20])
    #print("   y: ",y[0:20])
    #tck = interpolate.splrep(x, y, s=0)
    '''
    points = list(zip(x, y))
    #print("   Removing duplicates...")
    x1 = []
    y1 = []
    for ip,p in enumerate(points):
      if p[0] > 0.29999:
        if p[0] > points[ip-1][0]:
          x1.append(p[0])
          y1.append(p[1])
    x = np.array(x1)
    y = np.array(y1)
    '''
    #print("   Done eemoving duplicates")        
    '''
    if i==0:
      print("   setting reference_beff")
      xnew = x      
    #print("   x: ",x[0:20])
    #print("   y: ",y[0:20])
    #ynew = interpolate.splev(xnew, tck, der=0)
      
    f2 = interpolate.interp1d(x, y, kind='cubic',fill_value="extrapolate")
    ynew=f2(xnew)
    #print("   ynew: ",ynew[0:20])
    '''
    if i==0:
      reference_ynew = y
    #print("   reference_ynew: ",reference_ynew[0:20])
    ljr_ratio = np.divide(y, reference_ynew)
    #print("   ljr_ratio: ",ljr_ratio[0:20])
    print(condition[key]['label'])
    plt.plot(x, ljr_ratio, color=condition[key]['color'])

#plt.legend(loc='best') 
plt.xlim([0.3,1.0])
#plt.ylim([0.1,2])
plt.ylabel('Ratio to ITK')
plt.xlabel('b-jet efficiency')

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

