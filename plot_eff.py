import pickle, sys
import glob, os
import re
import ROOT
from array import array
import time
import json
import math
import types


ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

def make_leg(n, labels, x1=0.7, y1=0.6, x2=0.876, y2=0.87, hdata=None, textSize=0):
    if len(n)<4:
        y1 = y1 + (y2-y1)*(4.0-len(n))/4.0
    print "y1",y1
    leg=ROOT.TLegend(x1,y1,x2,y2)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetLineWidth(0)
    if textSize>0:
        leg.SetTextFont(43)
        leg.SetTextSize(textSize)
    i=len(n)
    if not hdata is None:
        leg.AddEntry(hdata, "Data")
    for h in reversed(n):
        print "hey"
        i-=1
        if "hh" in labels[i]:
            continue
        if h.GetFillColor() or "ulti" in labels[i] :
            leg.AddEntry(h,labels[i],"f")
        else:
            leg.AddEntry(h,labels[i].replace("GGM_","").replace("_"," "),"l")
    i=len(n)
    for h in reversed(n):
        i-=1
        if not "hh" in labels[i]:
            continue
        if h.GetFillColor() :
            leg.AddEntry(h,labels[i],"f")
        else:
            leg.AddEntry(h,labels[i].replace("GGM_","").replace("_"," "),"l")

    return leg

if __name__ == "__main__":
    print "Ciao Chiara"

    CutBase=" jet_pt>20e3 && jet_truthMatch==1 && jet_isPU==0 && abs(PVz-truth_PVz)<0.1";
    
    CutFlav=dict()
    CutFlav["B"]=" && jet_LabDr_HadF==5 "
    CutFlav["C"]=" && jet_LabDr_HadF==4 "
    CutFlav["L"]=" && (jet_LabDr_HadF!=4 && jet_LabDr_HadF!=5 && jet_LabDr_HadF!=15) && jet_dRminToB>0.8 && jet_dRminToC>0.8 && jet_dRminToT>0.8"

    t = ROOT.TChain("bTag_AntiKt4EMTopoJets")
    t.Add("/eos/user/c/crizzi/HGTD/btagging/output/ITK/user.crizzi.mc15_14TeV.117050.PowhegPythia_P2011C_ttbar.recon.AOD.e2176_s3348_s3347_r10900_r11003.btag_ITKonly_v0_Akt4EMTo/user.crizzi.16988484.Akt4EMTo._0014*")
    var = ("fabs(jet_eta)",10,0.,4)
    name_can="eff_btag"
    title_x_axis = "|#eta|"
    setups = ["ITK"]

    flavours = ["B","C","L"]
    colors = [410, 856, 607, 801, 629, 879, 602, 921, 622]
    for flav in flavours:
        sel = CutBase+CutFlav[flav]
        n_allbkg=[]
        for b in setups:
            n_thisbkg=[]
            i=0
            for s in ["1","jet_ip3d_llr>2.007"]:
                print s
                htmp_name = name_can+"_"+str(i)+"_"+b+flav
                htmp=ROOT.TH1D(htmp_name, htmp_name, var[1], var[2], var[3])
                htmp.GetXaxis().SetTitle(title_x_axis)
                htmp.Sumw2()
                sel_slice=sel+" && ("+s+")"
                string_draw = var[0]+">>"+name_can+"_"+str(i)+"_"+b+flav
                print string_draw, sel_slice
                t.Draw(string_draw,sel_slice,"goff")
                print htmp.Integral()
                n_thisbkg.append(htmp.Clone())
                i+=1
            n_allbkg.append(n_thisbkg)

        heff=[]
        i=0
        for b in setups:
            heff.append(n_allbkg[i][1].Clone())
            heff[i].Divide(n_allbkg[i][0])
            heff[i].SetLineColor(colors[i])
            i+=1


        c = ROOT.TCanvas("can"+name_can+flav,"can"+name_can,600,600)
        pad1 = ROOT.TPad("pad1", "pad1",0.0,0.35,1.0,1.0,21)
        pad2 = ROOT.TPad("pad2", "pad2",0.0,0.0,1.0,0.35,22)
        pad2.SetGridy()
        pad1.SetFillColor(0)
        pad1.Draw()
        pad2.SetFillColor(0)
        pad2.Draw()

        pad1.cd()
        i=0
        for b in setups:
            if i==0:
                heff[i].GetYaxis().SetTitleFont(43)
                heff[i].GetYaxis().SetTitle("B-tagging eff")
                heff[i].GetYaxis().SetTitleSize(19)
                heff[i].GetYaxis().SetLabelFont(43)
                heff[i].GetYaxis().SetLabelSize(15)
                heff[i].GetYaxis().SetTitleOffset(1.3)        
                
                heff[i].GetXaxis().SetTitle(title_x_axis)
                heff[i].GetXaxis().SetTitleFont(43)
                heff[i].GetXaxis().SetTitleSize(19)
                heff[i].GetXaxis().SetLabelFont(43)
                heff[i].GetXaxis().SetLabelSize(15)
                heff[i].GetXaxis().SetTitleOffset(1.5)    
                
                heff[i].Draw("histo")
            else:
                heff[i].Draw("histo same")
            i+=1

        leg = make_leg(heff, setups)
        leg.Draw()
        pad1.Update()

        hratio = []
        i=0
        for b in setups:
            hratio.append(heff[i].Clone())
            hratio[i].Divide(heff[0])
            i+=1

        pad2.cd()
        i=0
        for b in setups:
            if i==0:
                hratio[i].SetFillColor(622)
                hratio[i].SetFillStyle(3244)
                hratio[i].SetMaximum(1.8)
                hratio[i].SetMinimum(0.2)
                hratio[i].GetXaxis().SetTitle("")
                
                hratio[i].GetYaxis().SetTitleFont(43)
                hratio[i].GetYaxis().SetTitleSize(19)
                
                hratio[i].GetYaxis().SetLabelFont(43)
                hratio[i].GetYaxis().SetLabelSize(15)
                hratio[i].GetYaxis().SetTitleOffset(1.3) 
                
                hratio[i].GetYaxis().SetTitle("Ratio to ITK")
                hratio[i].Draw("E2")        
            else:
                hratio[i].Draw("same")

        c.SaveAs(name_can+"_"+flav+".pdf")
