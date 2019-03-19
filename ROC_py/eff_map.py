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
parser.add_argument('--largeEta', type=int, default=0)
parser.add_argument('--twoTrk', type=int, default=0)
parser.add_argument('--bEff', type=int, default=1)

args = parser.parse_args()

tagger=args.tagger
largeEta=args.largeEta>0
smallEta=args.largeEta==-1
twoTrk = args.twoTrk>0
bEff = args.bEff>0

print("largeEta:",largeEta)
print("smallEta:",smallEta)

def getCut( tagger ):
  if tagger=="MV1":  return  0.945487725
  if tagger=="MV1c": return  0.779833333333
  if tagger=="MV2c00": return  0.0308333333333
  if tagger=="MV2c10": return  -0.00416666666667
  if tagger=="MV2c20" : return  -0.0215
  if tagger=="IP3D":      return  2.007
  if tagger=="IP3D+SV1":  return  4.3625
  if tagger=="MVb":       return  -0.120991666667
  if tagger=="SV1":       return  -97
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
  if tagger=="IP3D+SV1":  return  "jet_sv1ip3d"
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


#keys = ["ITK", "Initial", "Int1", "Int2", "Final"]
keys = ["ITK"]

condition= {
    "ITK":{
    #"file":"/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/ITK/file.root",
         "file":"../input/ITK/file.root",
         #"file":"../input/compare/Chiara_ntuple.root",
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


def get_df(t):  
    df = t.pandas.df(["jet_pt","jet_eta","eventnb","jet_LabDr_HadF", "jet_dRminToB", 
                      "jet_dRminToC", "jet_dRminToT", "jet_isPU", "jet_truthMatch", "PVz", 
                      "truth_PVz", getVariableNtrk(tagger), getVariable(tagger)])
    df["jet_pt"] = df["jet_pt"]/1000.
    df=df.rename_axis(index=['entry', 'jet_pos'])
    df=df.reset_index(level=['jet_pos'])
    df = df.set_index(["eventnb","jet_pos"])    
    CutBase=" jet_pt>20 & jet_truthMatch==1 & jet_isPU==0 & abs(PVz-truth_PVz)<0.1  & abs(jet_eta)<4"
    if largeEta:
      CutBase = CutBase+" & abs(jet_eta)>2.4"
      print("add large eta requirement")
    if smallEta:
      CutBase = CutBase+" & abs(jet_eta)<2.4"
      print("add small eta requirement")
    CutBase = CutBase+" & ( (jet_LabDr_HadF!=4 & jet_LabDr_HadF!=5 & jet_LabDr_HadF!=15 & jet_dRminToB>0.8 & jet_dRminToC>0.8 & jet_dRminToT>0.8)"
    CutBase = CutBase+" | (jet_LabDr_HadF==5))"

    df = df.query(CutBase)

    print(df.shape)
    return df

def filter_df(df1, ITKcomp=0, twoTrk=0, df2=0):
    print("shape df1 input: ", df1.shape)
    if ITKcomp:
        # index in common with ITK-only
        idx = sorted(df1.index.intersection(df2.index)) 
        df1 = df1.loc[idx] 
        df2 = df2.loc[idx]
        if twoTrk:
            # < 2 trk: use ITK, >=2 trk: use HGTS
            df2_sel = df2.query(getVariableNtrk(tagger)+" < 2")
            df1_sel = df1.loc[df1.index.difference(df2_sel.index)]
            df3 = pd.concat([df1_sel,df2_sel])
            return df3
        else:
            return df1
    else:      
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
    else:
      df = filter_df(df, 1, twoTrk, df_ref)

    print("So far so good!")

    
    df['jet_eta_bin'] = pd.cut(df['jet_eta'].abs(), bins=[0, 0.5, 1.0, 1.5, 2.0, 2.4, 2.8, 3.2, 4.0])
    df['jet_pt_bin'] =  pd.cut(df['jet_pt'],        bins=[0, 20, 30, 40, 60, 85, 110, 150, 250])
    

    def get_ratio(g, invert=False):
      num = float(g.loc[g[getVariable(tagger)] > getCut(tagger), :].shape[0])
      den = float(g.shape[0])
      if invert:
        num,den = den,num
      try:
        return num/den
      except ZeroDivisionError:
        return np.nan

    eff = df.loc[df["jet_LabDr_HadF"]==5].groupby(['jet_eta_bin','jet_pt_bin']).apply(get_ratio)

    eff = eff.to_frame()
    eff.columns = ['eff']
    #entries_tagged = entries_tagged.to_frame()
    #print(eff.to_string())

    def print_tf2(eff):
      tf2_string = "Double_t func(Double_t *val, Double_t *par) { \n"
      tf2_string += "Float_t pt = val[0]; \n"
      tf2_string += "Float_t eta = val[1]; \n"
      tf2_string += "Double_t eff =0; \n"
      for index, row in eff.iterrows():
        pt_low = index[1].left
        pt_high = index[1].right
        eta_low = index[0].left
        eta_high = index[0].right
        value = row["eff"]
        tf2_string += 'if (pt > {} && pt <= {} && eta > {} && eta <= {})  eff = {}; \n'.format(pt_low, pt_high, eta_low, eta_high, value)
      tf2_string += "\n"
      tf2_string += "return eff; \n}"

      return tf2_string
    
    print(print_tf2(eff))


"""    
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
if largeEta:
    plt.text(0.02, 0.10, 'jet p$_{T}>$ 20 GeV, $|\eta|>$2.4', size='medium',transform=ax1.transAxes)
elif smallEta:
  plt.text(0.02, 0.10, 'jet p$_{T}>$ 20 GeV, $|\eta|<$2.4', size='medium',transform=ax1.transAxes)
else:
    plt.text(0.02, 0.10, 'jet p$_{T}>$ 20 GeV', size='medium',transform=ax1.transAxes)
plt.grid(True)

ax2 = plt.subplot(gs[1])
for i,key in enumerate(keys):
    print('\n')
    print(key)
    y=fpr_all[key] # light rejection
    x=tpr_all[key]
    if i==0:
      reference_ynew = y
    #print('   reference_ynew: ',reference_ynew[0:20])
    ljr_ratio = np.divide(y, reference_ynew)
    #print('   ljr_ratio: ',ljr_ratio[0:20])
    print(condition[key]['label'])
    plt.plot(x, ljr_ratio, color=condition[key]['color'])

plt.xlim([0.,4.])
plt.ylim([0.8,3])
plt.ylabel('Ratio to ITK')
plt.xlabel('jet eta')

plt.grid(True)

if bEff:
  name_plot = 'eff_'+getVariable(tagger)
else:
  name_plot = 'mistag_'+getVariable(tagger)
if largeEta:
  name_plot+='_largeEta'
if smallEta:
  name_plot+='_smallEta'
if twoTrk:
  name_plot+='_twoTrk'
name_plot+='_trkEffselfTag'
plt.savefig(name_plot+'.pdf')
"""   

