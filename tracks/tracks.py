import ROOT

infile = "/afs/cern.ch/user/c/crizzi/myeos/HGTD/btagging/input_eff_plot/ITK/file.root"
f = ROOT.TFile.Open(infile, "READ")
t = f.Get("bTag_AntiKt4EMTopoJets")

ROOT.gStyle.SetOptStat(0)
#ip3d: 100, -100, 100
#flav: 17, -0.5, 16.5
#sv1: 100, -100, 5500
#ntrk: 31, -0.5, 30.5

# jet_btag_ntrk:jet_sv1_ntrk:jet_ip3d_ntrk
# jet_ip3d_llr jet_sv1_m

selections = {
    "all":" jet_pt>20e3 && jet_truthMatch==1 && jet_isPU==0 && abs(PVz-truth_PVz)<0.1 && abs(jet_eta)<4",  
    "all_eta_larger2p4":" jet_pt>20e3 && jet_truthMatch==1 && jet_isPU==0 && abs(PVz-truth_PVz)<0.1 && abs(jet_eta)<4 && abs(jet_eta)>2.4",  
    "real_b":" jet_pt>20e3 && jet_truthMatch==1 && jet_isPU==0 && abs(PVz-truth_PVz)<0.1 && abs(jet_eta)<4 && jet_LabDr_HadF==5",  
    "real_b_eta_larger2p4":" jet_pt>20e3 && jet_truthMatch==1 && jet_isPU==0 && abs(PVz-truth_PVz)<0.1 && abs(jet_eta)<4 && abs(jet_eta)>2.4 && jet_LabDr_HadF==5",  
    "real_light":" jet_pt>20e3 && jet_truthMatch==1 && jet_isPU==0 && abs(PVz-truth_PVz)<0.1 && abs(jet_eta)<4  && (jet_LabDr_HadF!=4 && jet_LabDr_HadF!=5 && jet_LabDr_HadF!=15) && jet_dRminToB>0.8 && jet_dRminToC>0.8 && jet_dRminToT>0.8",  
    "real_light_eta_larger2p4":" jet_pt>20e3 && jet_truthMatch==1 && jet_isPU==0 && abs(PVz-truth_PVz)<0.1 && abs(jet_eta)<4  && (jet_LabDr_HadF!=4 && jet_LabDr_HadF!=5 && jet_LabDr_HadF!=15) && jet_dRminToB>0.8 && jet_dRminToC>0.8 && jet_dRminToT>0.8 && abs(jet_eta)>2.4"
    }

for sel in selections:
    print sel
    # trk_sv1 vs sv1 
    h_jet_sv1_ntrk_jet_sv1_m = ROOT.TH2F("h_jet_sv1_ntrk_jet_sv1_m","h_jet_sv1_ntrk_jet_sv1_m",100, -100, 5000, 31, -0.5, 30.5)
    # trk_ip3d vs ip3d
    h_jet_ip3d_ntrk_jet_ip3d_llr = ROOT.TH2F("h_jet_ip3d_ntrk_jet_ip3d_llr","h_jet_ip3d_ntrk_jet_ip3d_llr",100, -100, 100, 31, -0.5, 30.5)
    # trk_sv1 vs trk_ip3d
    h_jet_ip3d_ntrk_jet_sv1_ntrk = ROOT.TH2F("h_jet_ip3d_ntrk_jet_sv1_ntrk","h_jet_ip3d_ntrk_jet_sv1_ntrk", 31, -0.5, 30.5, 31, -0.5, 30.5)
    # trk_btagging vs trk_ip3d
    h_jet_ip3d_ntrk_jet_btag_ntrk = ROOT.TH2F("h_jet_ip3d_ntrk_jet_btag_ntrk","h_jet_ip3d_ntrk_jet_btag_ntrk", 31, -0.5, 30.5, 31, -0.5, 30.5)

    c_jet_ip3d_ntrk_jet_ip3d_llr = ROOT.TCanvas();
    t.Draw("jet_ip3d_ntrk:jet_ip3d_llr>>h_jet_ip3d_ntrk_jet_ip3d_llr",selections[sel],"goff")
    c_jet_ip3d_ntrk_jet_ip3d_llr.cd()
    h_jet_ip3d_ntrk_jet_ip3d_llr.GetXaxis().SetTitle("jet_ip3d_llr")
    h_jet_ip3d_ntrk_jet_ip3d_llr.GetYaxis().SetTitle("jet_ip3d_ntrk")
    h_jet_ip3d_ntrk_jet_ip3d_llr.SetTitle(sel)
    h_jet_ip3d_ntrk_jet_ip3d_llr.Draw("COLZ")
    c_jet_ip3d_ntrk_jet_ip3d_llr.SaveAs("jet_ip3d_ntrk_jet_ip3d_llr_"+sel+".pdf")

    c_jet_ip3d_ntrk_jet_sv1_ntrk = ROOT.TCanvas();
    t.Draw("jet_ip3d_ntrk:jet_sv1_ntrk>>h_jet_ip3d_ntrk_jet_sv1_ntrk",selections[sel],"goff")
    c_jet_ip3d_ntrk_jet_sv1_ntrk.cd()
    h_jet_ip3d_ntrk_jet_sv1_ntrk.GetXaxis().SetTitle("jet_sv1_ntrk")
    h_jet_ip3d_ntrk_jet_sv1_ntrk.GetYaxis().SetTitle("jet_ip3d_ntrk")
    h_jet_ip3d_ntrk_jet_sv1_ntrk.SetTitle(sel)
    h_jet_ip3d_ntrk_jet_sv1_ntrk.Draw("COLZ")
    c_jet_ip3d_ntrk_jet_sv1_ntrk.SaveAs("jet_ip3d_ntrk_jet_sv1_ntrk_"+sel+".pdf")

    c_jet_ip3d_ntrk_jet_btag_ntrk = ROOT.TCanvas();
    t.Draw("jet_ip3d_ntrk:jet_btag_ntrk>>h_jet_ip3d_ntrk_jet_btag_ntrk",selections[sel],"goff")
    c_jet_ip3d_ntrk_jet_btag_ntrk.cd()
    h_jet_ip3d_ntrk_jet_btag_ntrk.GetXaxis().SetTitle("jet_btag_ntrk")
    h_jet_ip3d_ntrk_jet_btag_ntrk.GetYaxis().SetTitle("jet_ip3d_ntrk")
    h_jet_ip3d_ntrk_jet_btag_ntrk.SetTitle(sel)
    h_jet_ip3d_ntrk_jet_btag_ntrk.Draw("COLZ")
    c_jet_ip3d_ntrk_jet_btag_ntrk.SaveAs("jet_ip3d_ntrk_jet_btag_ntrk_"+sel+".pdf")
        
    c_jet_sv1_ntrk_jet_sv1_m = ROOT.TCanvas();
    t.Draw("jet_sv1_ntrk:jet_sv1_m>>h_jet_sv1_ntrk_jet_sv1_m",selections[sel],"goff")
    c_jet_sv1_ntrk_jet_sv1_m.cd()
    h_jet_sv1_ntrk_jet_sv1_m.GetXaxis().SetTitle("jet_sv1_m")
    h_jet_sv1_ntrk_jet_sv1_m.GetYaxis().SetTitle("jet_sv1_ntrk")
    h_jet_sv1_ntrk_jet_sv1_m.SetTitle(sel)
    h_jet_sv1_ntrk_jet_sv1_m.Draw("COLZ")
    c_jet_sv1_ntrk_jet_sv1_m.SaveAs("jet_sv1_ntrk_jet_sv1_m_"+sel+".pdf")
    
