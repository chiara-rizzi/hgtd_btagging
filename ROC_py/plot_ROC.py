import uproot
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
import argparse
from scipy import interpolate
import collections
import pandas as pd
import sys

print("Ciao Chiara!")

parser = argparse.ArgumentParser(description='Make ROC curves.')
parser.add_argument('--tagger', type=str, default="MV1")
parser.add_argument('--etaRange', type=float , nargs="+", default=[2.4, 4])
parser.add_argument('--twoTrk', type=int, default=1)
parser.add_argument('--onlyCommonEvents', type=int, default=1)

args = parser.parse_args()

tagger=args.tagger
twoTrk = args.twoTrk>0
onlyCommonEvents = args.onlyCommonEvents > 0

print("etaRange:",args.etaRange)

def getVariable( tagger ):
  if tagger=="MV1":    return "jet_mv1"
  if tagger=="MV1c":   return "jet_mv1c"
  if tagger=="MV2c00": return "jet_mv2c00"
  if tagger=="MV2c10": return "jet_mv2c10"
  if tagger=="MV2c20": return "jet_mv2c20"
  if tagger=="MVb":    return "jet_mvb"
  if tagger=="IP3D":   return  "jet_ip3d_llr"
  if tagger=="IP3D+SV1":  return  "jet_my_sv1ip3d"
  if tagger=="IP3D+SV1old":  return  "jet_sv1ip3d"
  if tagger=="SV1":       return "jet_sv1_llr" 
  if tagger=="JetFitter": return  "jet_jf_llr"
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
  if tagger=="IP3D+SV1old":  return  "jet_sv1_ntrk"
  if tagger=="SV1":       return "jet_sv1_ntrk" 
  if tagger=="JetFitter": return  "jet_ip3d_ntrk"
  return "0"


keys = ["ITK", "Initial", "Int1", "Int2", "Final"]
#keys = ["ITK","Initial"]

condition = {
    "ITK":{
    #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/ITK/file.root",
         "file":"../input/AC_ITKonly_trkEff1hit/file.root",
         #"file":"../input/compare/Chiara_ntuple.root",
         "tree_name":"bTag_AntiKt4EMTopoJets",
         "label": "ITK",
         "color": "k"
    },
    "Int1":{
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Intermediate_pre_Replacement_trkEffselfTag/file.root",
        "file":"../input/AC_Intermediate_pre_Replacement_trkEff1hit/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": r'ITK+HGTD (2000 fb$^{-1}$)',
        "color": "y"
        },
    "Int2":{
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Intermediate_post_Replacement_trkEffselfTag/file.root",
        "file":"../input/AC_Intermediate_post_Replacement_trkEff1hit/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": r'ITK+HGTD (2001 fb$^{-1}$)',
        "color": "b"
        },
    "Initial":{
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Initial_trkEffselfTag/file.root",
        #"file":"../input/compare/Chiara_ntuple.root",
        "file":"../input/AC_Initial_trkEff1hit/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": r'ITK+HGTD (0 fb$^{-1}$)',
        "color": "g"
        },
    "Final":{
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Final_trkEffselfTag/file.root",
        "file":"../input/AC_Final_trkEff1hit/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": r'ITK+HGTD (4000 fb$^{-1}$)',
        "color": "r"
        }
    }


def compute_ip3d_sv1(row):
    if  row["jet_ip3d_pu"] >0 and row["jet_ip3d_pb"] >0:
        if  row["jet_ip3d_pu"] != 1 or  row["jet_ip3d_pu"] != 1.e9:
          IP3DPlusSV1w= math.log ( row["jet_ip3d_pb"] / row["jet_ip3d_pu"] ) # logarithm neperien
    else:
        IP3DPlusSV1w = 0
    if row["jet_sv1_ntrkv"]>0: #//SV1 found.
        if  row["jet_sv1_pu"] >0 and row["jet_sv1_pb"] >0 :
           IP3DPlusSV1w += math.log ( row["jet_sv1_pb"] / row["jet_sv1_pu"] );
    else:
        IP3DPlusSV1w += -1.55 # log(row["jet_sv1_pb"]/row["jet_sv1_pu"])  = log((1-eff_b)/(1-eff_u)) from the endpoint of the SV1 ROC curve
    return IP3DPlusSV1w


def get_df(t):  
    if tagger == "IP3D+SV1":
        var_needed_for_ip3d_sv1=["jet_sv1_pb",
                                 "jet_sv1_pu",
                                 "jet_sv1_ntrk",
                                 "jet_ip3d_ntrk",
                                 "jet_sv1_ntrkv",
                                 "jet_ip3d_pb",
                                 "jet_ip3d_pu"]
        df = t.pandas.df(["jet_pt","jet_eta","eventnb","jet_LabDr_HadF", "jet_dRminToB", 
                          "jet_dRminToC", "jet_dRminToT", "jet_isPU", "jet_truthMatch", "PVz", 
                          "truth_PVz", getVariableNtrk(tagger)]+var_needed_for_ip3d_sv1) # add all the variables needed to compute IP3D+SV1
        #df["jet_my_sv1ip3d"] = df.apply (lambda row: compute_ip3d_sv1(row), axis=1) # sloooow
    else:
        df = t.pandas.df(["jet_pt","jet_eta","eventnb","jet_LabDr_HadF", "jet_dRminToB", 
                          "jet_dRminToC", "jet_dRminToT", "jet_isPU", "jet_truthMatch", "PVz", 
                          "truth_PVz", getVariableNtrk(tagger), getVariable(tagger)])

    df=df.rename_axis(index=['entry', 'jet_pos'])
    df=df.reset_index(level=['jet_pos'])
    df = df.set_index(["eventnb","jet_pos"])    
    CutBase=" jet_pt>20000 & jet_truthMatch==1 & jet_isPU==0 & abs(PVz-truth_PVz)<0.1  & abs(jet_eta)>"+str(args.etaRange[0])
    CutBase = CutBase+" & abs(jet_eta)<"+str(args.etaRange[1])
    CutBase = CutBase+" & ( (jet_LabDr_HadF!=4 & jet_LabDr_HadF!=5 & jet_LabDr_HadF!=15 & jet_dRminToB>0.8 & jet_dRminToC>0.8 & jet_dRminToT>0.8)"
    CutBase = CutBase+" | (jet_LabDr_HadF==5))"

    df = df.query(CutBase)

    print(df.shape)
    return df

def common_index(df_list):
    idx = df_list[0].index
    for df in df_list[1:]:
        idx = idx.intersection(df.index)
    return idx

def filter_df(df1, ITKcomp=0, twoTrk=0, df2=0):
    print("shape df1 input: ", df1.shape)
    # va messo dopo che ho sistemato le cose 
    if ITKcomp:
        # index in common with ITK-only
        idx = sorted(df1.index.intersection(df2.index)) 
        df1 = df1.loc[idx] 
        df2 = df2.loc[idx]
        if twoTrk:
            # < 2 trk: use ITK, >=2 trk: use HGTS

            if tagger == "IP3D+SV1":
                df2.columns = [str(col) + '_ITK' for col in df2.columns]
                df3 = pd.concat([df1, df2], axis=1)
              #print(df3.columns)
                df3["jet_ip3d_pu"] = np.where(df3['jet_ip3d_ntrk_ITK']>1, df3["jet_ip3d_pu"], df3["jet_ip3d_pu_ITK"])
                df3["jet_ip3d_pb"] = np.where(df3['jet_ip3d_ntrk_ITK']>1, df3["jet_ip3d_pb"], df3["jet_ip3d_pb_ITK"])
                df3["jet_sv1_pu"] = np.where(df3['jet_sv1_ntrk_ITK']>1, df3["jet_sv1_pu"], df3["jet_sv1_pu_ITK"])
                df3["jet_sv1_pb"] = np.where(df3['jet_sv1_ntrk_ITK']>1, df3["jet_sv1_pb"], df3["jet_sv1_pb_ITK"])
                df3["jet_sv1_ntrkv"] = np.where(df3['jet_sv1_ntrk_ITK']>1, df3["jet_sv1_ntrkv"], df3["jet_sv1_ntrkv_ITK"])
                df3["jet_my_sv1ip3d"] = np.where((df3["jet_ip3d_pu"]>0) & (df3["jet_ip3d_pb"] >0) & (df3["jet_sv1_pu"]>0) & (df3["jet_sv1_pb"]>0), np.log ( df3["jet_ip3d_pb"] / df3["jet_ip3d_pu"] )+np.log ( df3["jet_sv1_pb"] / df3["jet_sv1_pu"] ),
                                                 np.where( (df3["jet_ip3d_pu"]>0) & (df3["jet_ip3d_pb"] >0) & (df3["jet_sv1_ntrkv"]<1), np.log ( df3["jet_ip3d_pb"] / df3["jet_ip3d_pu"] ) -1.55,
                                                           np.where( df3["jet_sv1_ntrkv"]<0, -1.55, 
                                                                     np.where( ((df3["jet_sv1_pu"]>0) & (df3["jet_sv1_pb"]>0)), np.log ( df3["jet_sv1_pb"] / df3["jet_sv1_pu"] ), 0)                                                                 
                                                                     )
                                                           )
                                                 ) # HERE
                itk_col = [col for col in df3.columns if "ITK" in col]
                df3=df3.drop(itk_col, axis=1)
                return df3
            else:             # previously was done also for IP3D+SV1 looking at the number of tracks for SV1

                df2_sel = df2.query(getVariableNtrk(tagger)+" < 2") 
                df1_sel = df1.loc[df1.index.difference(df2_sel.index)]
                df3 = pd.concat([df1_sel,df2_sel])   
                return df3
        else:
            if tagger == "IP3D+SV1":
                df1["jet_my_sv1ip3d"] = np.where((df1["jet_ip3d_pu"]>0) & (df1["jet_ip3d_pb"] >0) & (df1["jet_sv1_pu"]>0) & (df1["jet_sv1_pb"]>0), np.log ( df1["jet_ip3d_pb"] / df1["jet_ip3d_pu"] )+np.log ( df1["jet_sv1_pb"] / df1["jet_sv1_pu"] ),
                                                 np.where( (df1["jet_ip3d_pu"]>0) & (df1["jet_ip3d_pb"] >0) & (df1["jet_sv1_ntrkv"]<1), np.log ( df1["jet_ip3d_pb"] / df1["jet_ip3d_pu"] ) -1.55,
                                                           np.where( df1["jet_sv1_ntrkv"]<0, -1.55, 
                                                                     np.where( ((df1["jet_sv1_pu"]>0) & (df1["jet_sv1_pb"]>0)), np.log ( df1["jet_sv1_pb"] / df1["jet_sv1_pu"] ), 0))))
            return df1
    else:      
        if tagger == "IP3D+SV1":
            df1["jet_my_sv1ip3d"] = np.where((df1["jet_ip3d_pu"]>0) & (df1["jet_ip3d_pb"] >0) & (df1["jet_sv1_pu"]>0) & (df1["jet_sv1_pb"]>0), np.log ( df1["jet_ip3d_pb"] / df1["jet_ip3d_pu"] )+np.log ( df1["jet_sv1_pb"] / df1["jet_sv1_pu"] ),
                                             np.where( (df1["jet_ip3d_pu"]>0) & (df1["jet_ip3d_pb"] >0) & (df1["jet_sv1_ntrkv"]<1), np.log ( df1["jet_ip3d_pb"] / df1["jet_ip3d_pu"] ) -1.55,
                                                       np.where( df1["jet_sv1_ntrkv"]<0, -1.55, 
                                                                 np.where( ((df1["jet_sv1_pu"]>0) & (df1["jet_sv1_pb"]>0)), np.log ( df1["jet_sv1_pb"] / df1["jet_sv1_pu"] ), 0))))

        return df1

def input_arrays(df):
  return df.as_matrix(columns=['jet_LabDr_HadF']), df.as_matrix(columns=[getVariable(tagger)])
  

gs = gridspec.GridSpec(2, 1,
                       height_ratios=[3,1],
                       hspace=0.10)
ax1 = plt.subplot(gs[0])
plt.title(tagger) #
fpr_all={}
tpr_all={}
df_all={}
for i,key in enumerate(keys):
    print(key)
    f = uproot.open(condition[key]["file"]) # open file
    t = f[condition[key]["tree_name"]] #
    df = get_df(t)
    if int(i)==0:
        df_ref = df.copy()
        df = filter_df(df, 1, 0, df_ref)
    else:
        df = filter_df(df, 1, twoTrk, df_ref)

    df_all[key] = df

if onlyCommonEvents:
    df_list=[]
    for i,key in enumerate(keys):
        df_list.append(df_all[key]) 
    idx = common_index(df_list)
    for i,key in enumerate(keys):
        df_all[key] =  df_all[key].loc[idx]

min_max_x = 1
x_all={}
y_all={}
#roc_auc_all={}
for i,key in enumerate(keys):
    df = df_all[key]
    print(key)
    print(df.shape)
    labels, scores = input_arrays(df)
    fpr, tpr, thresholds = metrics.roc_curve(labels, scores, 5) 
    #roc_auc = metrics.auc(tpr, fpr)
    #print(roc_auc)
    #roc_auc_all[key] = roc_auc
    x = tpr
    y = 1/fpr
    if x[-2]<min_max_x: min_max_x = x[-2] # for SV1 and the other ones with endpoint
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
    x_all[key]=x
    y_all[key]=y

xnew = xnew[xnew<=min_max_x]

for i,key in enumerate(keys):
    x = x_all[key]
    y = y_all[key]
    tck = interpolate.splrep(x, y)
    ynew = interpolate.splev(xnew, tck, der=0)

    fpr_all[key]=ynew
    tpr_all[key]=xnew

    print(condition[key]['label'])     
    #plt.semilogy(xnew, ynew, label=condition[key]['label']+' AUC = %0.4f'% roc_auc_all[key], color=condition[key]['color'])
    plt.semilogy(xnew, ynew, label=condition[key]['label'], color=condition[key]['color'])
    
plt.legend(loc='best') 
plt.xlim([0.5,1.0])
if args.etaRange[0] > 2:
  plt.ylim([1,600])
else:
  plt.ylim([1,3000])
plt.gca().set_xticklabels(['']*10)
#plt.xticks([])
plt.ylabel('light-jet rejection')
plt.xlabel('')
plt.text(0.02, 0.25, r'$\mathbf{ATLAS}$ Simulation Internal', size='large',transform=ax1.transAxes)
plt.text(0.02, 0.15, '$t\overline{t}$ simulation', size='medium',transform=ax1.transAxes)
plt.text(0.02, 0.10, 'jet p$_{T}>$ 20 GeV, '+ str(args.etaRange[0])+ '<$|\eta|<$'+str(args.etaRange[1]), size='medium',transform=ax1.transAxes)
plt.grid(True, which="both", linestyle='--')

ax2 = plt.subplot(gs[1])
for i,key in enumerate(keys):
    print("\n")
    print(key)
    y=fpr_all[key] # light rejection
    x=tpr_all[key]
    if i==0:
      reference_ynew = y
    #print("   reference_ynew: ",reference_ynew[0:20])
    ljr_ratio = np.divide(y, reference_ynew)
    #print("   ljr_ratio: ",ljr_ratio[0:20])
    print(condition[key]['label'])
    plt.plot(x, ljr_ratio, color=condition[key]['color'])

plt.xlim([0.5,1.0])
plt.ylim([0.85,1.7])
plt.ylabel('Ratio to ITK')
plt.xlabel('b-jet efficiency')

plt.grid(True, linestyle='--')

name_plot = "roc_"+getVariable(tagger)+"_eta_"+str(args.etaRange[0])+"_"+str(args.etaRange[1])
if twoTrk:
    name_plot+="_twoTrk"
name_plot+="_trkEffselfTag"
plt.savefig(name_plot+".pdf")   

