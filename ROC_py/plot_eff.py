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

print("Ciao Chiara!")

parser = argparse.ArgumentParser(description='Make efficiency plots.')
parser.add_argument('--tagger', type=str, default="SV1")
parser.add_argument('--twoTrk', type=int, default=0)
parser.add_argument('--bEff', type=int, default=1)

args = parser.parse_args()

tagger=args.tagger
twoTrk = args.twoTrk>0
bEff = args.bEff>0


def getCut( tagger ):
  if tagger=="MV1":  return  0.945487725
  if tagger=="IP3D":      return  2.582 # 70: 2.431553, 85: -0.149126
  if tagger=="IP3D+SV1":  return  2.369 # 70: 2.133121, 85: -1.572315
  if tagger=="SV1":       return  -1.129 # 70: -1.358764, 85: -98
  if tagger=="JetFitter": return  -1.612 
  if tagger=="JetFitter": return  -1.6125 
  return 0


def getVariable( tagger ):
  if tagger=="MV1":    return "jet_mv1"
  if tagger=="MV1c":   return "jet_mv1c"
  if tagger=="MV2c00": return "jet_mv2c00"
  if tagger=="MV2c10": return "jet_mv2c10"
  if tagger=="MV2c20": return "jet_mv2c20"
  if tagger=="MVb":    return "jet_mvb"
  if tagger=="IP3D":   return  "jet_ip3d_llr"
  if tagger=="IP3D+SV1":  return  "jet_my_sv1ip3d"
  if tagger=="SV1":       return "jet_sv1_llr" 
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
         "file":"../input/AC_ITKonly_trkEff1hit/file.root",
         #"file":"../input/compare/Chiara_ntuple.root",
         "tree_name":"bTag_AntiKt4EMTopoJets",
         "label": "ITK-only",
         "color": "k"
    },
    "Int1":{
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Intermediate_pre_Replacement_trkEffselfTag/file.root",
        "file":"../input/AC_Intermediate_pre_Replacement_trkEff1hit/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": "Int. pre-repl",
        "color": "y"
        },
    "Int2":{
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Intermediate_post_Replacement_trkEffselfTag/file.root",
        "file":"../input/AC_Intermediate_post_Replacement_trkEff1hit/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": "Int. post-repl",
        "color": "b"
        },
    "Initial":{
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Initial_trkEffselfTag/file.root",
        #"file":"../input/compare/Chiara_ntuple.root",
        "file":"../input/AC_Initial_trkEff1hit/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": "Initial",
        "color": "g"
        },
    "Final":{
        #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/Final_trkEffselfTag/file.root",
        "file":"../input/AC_Final_trkEff1hit/file.root",
        "tree_name":"bTag_AntiKt4EMTopoJets",
        "label": "Final",
        "color": "r"
        }
    }


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
    CutBase=" jet_pt>20000 & jet_truthMatch==1 & jet_isPU==0 & abs(PVz-truth_PVz)<0.1"
    CutBase = CutBase+" & abs(jet_eta)<4"
    CutBase = CutBase+" & ( (jet_LabDr_HadF!=4 & jet_LabDr_HadF!=5 & jet_LabDr_HadF!=15 & jet_dRminToB>0.8 & jet_dRminToC>0.8 & jet_dRminToT>0.8)"
    CutBase = CutBase+" | (jet_LabDr_HadF==5))"

    df = df.query(CutBase)

    print(df.shape)
    return df

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
                       hspace=0.05)
ax1 = plt.subplot(gs[0])
plt.title(tagger) #
eff_all={}
tpr_all={}
fpr_all={}
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

    print("So far so good!")

    
    df['jet_eta_bin'] = pd.cut(df['jet_eta'].abs(), bins=[0,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0])
    #bins = df.groupby('jet_eta_bin').apply(lambda g: (g.shape[0], g.loc[g[getVariable(tagger)] > getCut(tagger), :].shape[0]))
    #print(bins)

    def get_ratio(g, invert=False):
      num = float(g.loc[g[getVariable(tagger)] > getCut(tagger), :].shape[0])
      den = float(g.shape[0])
      if invert:
        num,den = den,num
      try:
        return num/den
      except ZeroDivisionError:
        return np.nan

    if bEff:
      eff = df.loc[df["jet_LabDr_HadF"]==5].groupby('jet_eta_bin').apply(get_ratio)
    else:
      eff = df.loc[df["jet_LabDr_HadF"]==0].groupby('jet_eta_bin').apply(get_ratio, invert=True)

    eff = eff.to_frame()
    eff.columns = ['eff']
    eff["jet_eta_bin"] = eff.index
    eff["jet_eta_plot"] = eff["jet_eta_bin"].apply(lambda g: (g.right+g.left)/2)
    #print(df)
    
    # here
    x = eff.as_matrix(columns=['jet_eta_plot']) 
    y = eff.as_matrix(columns=['eff'])
    """
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
    """
    #if i==0:
    #  print("   setting reference_beff")
    #  xnew = x      
    #print("   x: ",x[0:20])
    #print("   y: ",y[0:20])
    #tck = interpolate.splrep(x, y)
    #ynew = interpolate.splev(xnew, tck, der=0)
    tpr_all[key]=x
    fpr_all[key]=y
    #roc_auc = metrics.auc(fpr, tpr)
    print(condition[key]['label']) 
    #print(roc_auc)
    print(y)
    plt.plot(x, y, label=condition[key]['label'], color=condition[key]['color'])
    
plt.legend(loc='best') 
plt.xlim([0.,4.])
if bEff:
  plt.ylim([0.,1.0])
  plt.ylabel('b-jet efficiency')
else:  
  max_ljr=max(max(fpr_all.values(), key = lambda g: max(g)))
  plt.ylim([0.,1.3*max_ljr])
  plt.ylabel('light rejection')
plt.xlabel('jet eta')
plt.text(0.02, 0.25, 'ATLAS Simulation Internal', size='large',transform=ax1.transAxes)
plt.text(0.02, 0.15, '$t\overline{t}$ simulation', size='medium',transform=ax1.transAxes)
plt.text(0.02, 0.10, 'jet p$_{T}>$ 20 GeV', size='medium',transform=ax1.transAxes)
plt.grid(True)

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

plt.xlim([0.,4.])
plt.ylim([0.8,3])
plt.ylabel('Ratio to ITK')
plt.xlabel('jet eta')

plt.grid(True)

if bEff:
  name_plot = "eff_"+getVariable(tagger)
else:
  name_plot = "mistag_"+getVariable(tagger)
if twoTrk:
  name_plot+="_twoTrk"
name_plot+="_trkEffselfTag"
plt.savefig(name_plot+".pdf")   

