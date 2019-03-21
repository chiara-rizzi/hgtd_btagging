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
parser.add_argument('--twoTrk', type=int, default=1)
parser.add_argument('--op', type=int, default=70)

args = parser.parse_args()

tagger=args.tagger
twoTrk = args.twoTrk>0


def getCut70( tagger ):
  if tagger=="MV1":  return  0.945487725
  if tagger=="IP3D":      return  2.582 # 70: 2.431553, 85: -0.149126
  if tagger=="IP3D+SV1":  return  2.369 # 70: 2.133121, 85: -1.572315
  if tagger=="SV1":       return  -1.129 # 70: -1.358764, 85: -98
  if tagger=="JetFitter": return  -1.612 
  return 0

def getCut85( tagger ):
  if tagger=="MV1":  return  0.945487725
  if tagger=="IP3D":      return  -0.008 # 70: 2.431553, 85: -0.149126
  if tagger=="IP3D+SV1":  return  -1.391 # 70: 2.133121, 85: -1.572315
  if tagger=="SV1":       return  -97 # 70: -1.358764, 85: -98
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
  if tagger=="IP3D+SV1old":  return  "jet_sv1ip3d"
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

def print_tf2(eff, key, num):
  tf2_string = "Double_t fitfun_mu200_{}_{}_{}_{}(Double_t *val, Double_t *par)".format(tagger.replace('+','').lower(), args.op, num, key)
  tf2_string+= "{ \n"
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
  tf2_string += "return eff; \n}\n\n\n"
  return tf2_string

def get_ratio(g, invert=False):
  if args.op == 70:
    num = float(g.loc[g[getVariable(tagger)] > getCut70(tagger), :].shape[0])
  elif args.op == 85:
    num = float(g.loc[g[getVariable(tagger)] > getCut85(tagger), :].shape[0])
  else: 
    print("OP not supported")
  den = float(g.shape[0])
  if invert:
    num,den = den,num
  try:
    return num/den
  except ZeroDivisionError:
    return np.nan

#keys = ["ITK", "Initial", "Int1", "Int2", "Final"]
keys = ["ITK","Initial"]

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


def compute_ip3d_sv1(row):
    if  row["jet_ip3d_pu"] >0 and row["jet_ip3d_pb"] >0:
        if  row["jet_ip3d_pu"] != 1 or  row["jet_ip3d_pu"] != 1.e9:
          IP3DPlusSV1w= np.log ( row["jet_ip3d_pb"] / row["jet_ip3d_pu"] ) # logarithm neperien
    else:
        IP3DPlusSV1w = 0
    if row["jet_sv1_ntrkv"]>0: #//SV1 found.
        if  row["jet_sv1_pu"] >0 and row["jet_sv1_pb"] >0 :
           IP3DPlusSV1w += np.log ( row["jet_sv1_pb"] / row["jet_sv1_pu"] );
    else:
        IP3DPlusSV1w += -1.55 # log(row["jet_sv1_pb"]/row["jet_sv1_pu"])  = log((1-eff_b)/(1-eff_u)) from the endpoint of the SV1 ROC curve
    return IP3DPlusSV1w


def get_df(t):  
    if tagger == "IP3D+SV1":
        var_needed_for_ip3d_sv1=["jet_sv1_pb",
                                 "jet_sv1_pu",
                                 "jet_sv1_ntrkv",
                                 "jet_sv1_ntrk",
                                 "jet_ip3d_ntrk",
                                 "jet_ip3d_pb",
                                 "jet_ip3d_pu"]
        df = t.pandas.df(["jet_pt","jet_eta","eventnb","jet_LabDr_HadF", "jet_dRminToB", 
                          "jet_dRminToC", "jet_dRminToT", "jet_isPU", "jet_truthMatch", "PVz", 
                          "truth_PVz", getVariableNtrk(tagger)]+var_needed_for_ip3d_sv1) # add all the variables needed to compute IP3D+SV1
        #df["jet_my_sv1ip3d"] = df.apply (lambda row: compute_ip3d_sv1(row), axis=1)
    else:
        df = t.pandas.df(["jet_pt","jet_eta","eventnb","jet_LabDr_HadF", "jet_dRminToB", 
                          "jet_dRminToC", "jet_dRminToT", "jet_isPU", "jet_truthMatch", "PVz", 
                          "truth_PVz", getVariableNtrk(tagger), getVariable(tagger)])

    df["jet_pt"] = df["jet_pt"]/1000.
    df=df.rename_axis(index=['entry', 'jet_pos'])
    df=df.reset_index(level=['jet_pos'])
    df = df.set_index(["eventnb","jet_pos"])    
    CutBase=" jet_pt>20 & abs(jet_eta)<4 & (jet_truthMatch==1 & jet_isPU==0 & abs(PVz-truth_PVz)<0.1"
    CutBase = CutBase+" & ( (jet_LabDr_HadF!=4 & jet_LabDr_HadF!=5 & jet_LabDr_HadF!=15 & jet_dRminToB>0.8 & jet_dRminToC>0.8 & jet_dRminToT>0.8)"
    CutBase = CutBase+" | (jet_LabDr_HadF==5) | (jet_LabDr_HadF==4))) | jet_isPU>0"

    df['flav'] = np.where(df['jet_isPU']>0, 'PU',
                          np.where( df['jet_LabDr_HadF'] ==5, 'B',
                                    np.where( df['jet_LabDr_HadF'] ==4, 'C', 'L'                                              
                                              )
                                    )
                          )

    df = df.query(CutBase)

    print(df.shape)
    #print(df)
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
  

num_flav = dict()
num_flav['L'] = '0'
num_flav['C'] = '1'
num_flav['B'] = '2'
num_flav['PU'] = '3'

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
    
    df['jet_eta_bin'] = pd.cut(df['jet_eta'].abs(), bins=[0, 0.5, 1.0, 1.5, 2.0, 2.4, 2.8, 3.1, 4.0])
    df['jet_pt_bin'] =  pd.cut(df['jet_pt'],        bins=[0, 20, 30, 40, 60, 85, 110, 150, 250])

    f= open('eff_files/efficiency_{}_{}_{}.C'.format(key, tagger.replace('+','').lower(), args.op),"w")
    for flav in num_flav:
        print(flav)
        eff = df.loc[df["flav"]==flav].groupby(['jet_eta_bin','jet_pt_bin']).apply(get_ratio)
        eff = eff.to_frame()
        eff.columns = ['eff']
        f.write(print_tf2(eff, key, num_flav[flav])) 
    f.close()



