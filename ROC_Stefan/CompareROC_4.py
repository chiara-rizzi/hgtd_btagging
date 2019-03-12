#import argparse
from glob import glob
import os, sys
from array import *
import ROOT 
from ROOT import * #TCanvas,TLegend,gROOT
import math
import time
from math import sqrt

ROOT.gROOT.LoadMacro("AtlasStyle.C")
ROOT.gROOT.LoadMacro("AtlasUtils.C")
ROOT.gROOT.SetBatch(True)
SetAtlasStyle()


tmpVale=sys.argv
def Helper():
    print " "
    print " Usage: python PrintROC.py file1 file2 leg1 leg2 outfolder "
    print " "
    sys.exit(1)
#if len(sys.argv) <6: Helper()

#odir = sys.argv[1]
#if ".root" in odir:
#    print " ..... outputfolder is a file .... please correct"
#    sys.exit(1)

leg =["ITk", "HGTD Initial","HGTD Int. pre-repl.","HGTD Int. post-repl.","HGTD Final"]

name1=""
name2=""
name3=""
name4=""
name5=""

def findStart(bin,target,curve):
    start=bin-10
    Tpx2=Double(0)
    Tpy2=Double(0)
    curve2.GetPoint(start,Tpx2,Tpy2)
    while Tpx2<target and Tpx2!=0 and start>-10:
        #print "      mah: "+str(start)+"  with: "+str(Tpx2)
        start-=10
        curve2.GetPoint(start,Tpx2,Tpy2)
    if start<0: start=0
    #print " ---> start is: "+str(start)+"  with: "+str(Tpx2)
    return start

print sys.argv
print len(sys.argv)
print sys.argv[-1]
files= []
infile1=sys.argv[1]
files.append(TFile(infile1,"R"))
if len(sys.argv) > 3:
    infile2=sys.argv[2]
    files.append(TFile(infile2,"R"))
if len(sys.argv) > 4:
    infile3=sys.argv[3]
    files.append(TFile(infile3,"R"))
if len(sys.argv) > 5:
    infile4=sys.argv[4]
    files.append(TFile(infile4,"R"))
if len(sys.argv) > 6:
    infile5=sys.argv[5]
    files.append(TFile(infile5,"R"))

odir= sys.argv[-1]
gSystem.Exec("mkdir -p "+odir)

taggerList=[] ### getting it from the first file
histoList= files[0].GetListOfKeys()
iterator=TIter(histoList)
obj=iterator()
while obj!=None:
    name=obj.GetName()
    if "---bc" in name:
        pass
    else:
        #if "MVb" not in name and "MV1" not in name:
        if "MV1" not in name: 
            taggerList.append( name.split("---")[0] )
    obj=iterator()
if name1=="MV1" or name2=="MV1":
    taggerList.append("MV1")
if name1=="MV1c" or name2=="MV1c":
    taggerList.append("MV1c")
print taggerList
if name1!="":
    taggerList=[]
    taggerList.append(name1)
if name2!="": taggerList.append(name2)
if name3!="": taggerList.append(name3)
if name4!="": taggerList.append(name4)
    
taggerList=list(set(taggerList))
#taggerList=["MV1","IP3D","IP2D","SV1","JetFitter","MV2c10","MV2c20"]
taggerList=["IP3D"]#,"IP3D","SV1"]

print taggerList

light=TH1F("b VS light","b VS light",100,0.5,1);
light.SetTitle(";b-jet efficiency;light-jet rejection;")
lightCurve=[]

cj=TH1F("b VS c","b VS c",100,0.5,1);
cj.SetTitle(";b-jet efficiency;c-jet rejection;")
cCurve=[]


myC=TCanvas( "bVSl", "bVSl",900,900);
myC.SetLogy()
myC.SetGridy()
myC.SetGridx()
light.SetMinimum(0.5)
light.SetMaximum(5e4)
light.GetXaxis().SetRangeUser(0.5,1.0)
light.Draw()
myLumi= "ATLAS Simulation Internal"
myLumi2= " t#bar{t} simulation, jet p_{T} > 20 GeV, |#eta|>2.4"
myText(0.20,0.24,1,myLumi,0.045)
myText(0.20,0.19,1,myLumi2,0.045)
legend4=TLegend(0.48,0.53,0.92,0.95)
legend4.SetTextFont(42)
legend4.SetTextSize(0.030)
legend4.SetFillColor(0)
legend4.SetLineColor(0)
legend4.SetFillStyle(0)
legend4.SetBorderSize(0)
count=1
colors = [kBlack, kBlue, kGreen, kRed, kYellow]
for tag in taggerList:
    curves=[]
    if name1!="" and name1!=tag: continue
    print tag
    myC.cd()
    for f in files:
        f.cd()
        curves.append(f.Get(tag+"---bl"))
    print curves
    ic=0
    for curve in curves:
        curve.SetLineStyle(1)
        curve.SetLineColor(colors[ic])
        curve.SetLineWidth(3)
        curve.Draw("C")
        ic+=1

    myCx=TCanvas( "bVSl"+tag, "bVSl"+tag,800,800);
    pad_1=TPad("pad_1", "up", 0., 0.35, 1., 1.);
    pad_1.SetBottomMargin(0);
    pad_1.Draw();   
    pad_1.SetGridy()
    pad_1.SetGridx()
    pad_1.SetLogy()
    pad_2=TPad("pad_2", "down", 0.0, 0.00, 1.0, 0.35);
    pad_2.SetTopMargin(0);
    pad_2.SetBottomMargin(0.28);
    pad_2.Draw();
    pad_1.cd();

    light.Draw()
    myLumi3=tag
    myText(0.20,0.85,1,myLumi3,0.045)
    myText(0.20,0.14,1,myLumi,0.045)
    myText(0.20,0.09,1,myLumi2,0.045)
    for curve in curves:
        curve.Draw("C")
    
    legend5=TLegend(0.55,0.7,0.92,0.90)
    legend5.SetTextFont(42)
    legend5.SetTextSize(0.042)
    legend5.SetFillColor(0)
    legend5.SetLineColor(0)
    legend5.SetFillStyle(0)
    legend5.SetBorderSize(0)
    for ic in range(len(curves)):
        legend5.AddEntry(curves[ic] ,leg[ic],"L")
    legend5.Draw("SAME")

    myCx.Update()
    
    pad_2.cd()
    ratio=light.Clone("ratio")
    ratio.GetYaxis().SetTitle("rel. diff.")
    ratio.GetYaxis().SetTitleOffset(1.1)
    ratio.GetXaxis().SetLabelSize(0.10)
    ratio.GetXaxis().SetTitleSize(0.10)
    ratio.GetYaxis().SetTitleOffset(0.7)
    ratio.GetYaxis().SetLabelSize(0.09)
    ratio.GetYaxis().SetTitleSize(0.09)
    ratio.SetMaximum(4.2)
    ratio.SetMinimum(0.3)
    ratio.SetLineColor(kGray)
    ratio.Draw("HIST")


    ratC=[]
    for i in range(len(curves)):
        ratC.append(TGraph())

    print " looping over: "+str(curves[0].GetN())+"    points"

    for ic in range(1,len(curves)):
        countPoint=0
        for bin in xrange(1,curves[0].GetN()+1):
            px1=Double(0)
            py1=Double(0)
            curves[0].GetPoint(bin-1,px1,py1)
            if bin%1500==0: print str(bin)+" ....: "+str(px1)+" , "+str(py1)
            if bin%3==0: continue
            if px1>0.999: continue
            if px1<light.GetXaxis().GetXmin(): break
            clo_py2=100000
            clo_px2=100000
            endV=0
            start=0
            clo_py2=curves[ic].Eval(px1)
            if clo_py2!=100000 :
                ratC[ic].SetPoint(countPoint,px1,clo_py2/py1)
                countPoint+=1

        ratC[ic].SetLineWidth(3)
        ratC[ic].SetLineColor(colors[ic])
        line1=TLine(light.GetXaxis().GetXmin(),1.,1.0,1.)
        line1.SetLineWidth(2)
        line1.SetLineStyle(2)
        line1.SetLineColor(922)
        line1.Draw("SAME")
        ratC[ic].Draw("SAMEL")
        myCx.Update()
    

    myCx.Update()
    myCx.Print(odir+"/bVSlight__"+tag+".pdf")
    myCx.Print(odir+"/bVSlight__"+tag+".eps")
    #myCx.Print(odir+"/bVSlight__"+tag+".png")
    #myCx.Print(odir+"/bVSlight__"+tag+".C")
    
    
myC.cd()
legend4.Draw()
myText(0.20,0.24,1,myLumi,0.045)
myText(0.20,0.19,1,myLumi2,0.045)
myC.Update()
#myC.Print(odir+"/bVSlightALL.eps")

################################################################

myC2=TCanvas( "cVSl", "cVSl",900,900);
myC2.SetLogy()
myC2.SetGridy()
myC2.SetGridx()
cj.SetMinimum(1)
cj.SetMaximum(2e2)
cj.GetXaxis().SetRangeUser(0.5,1.0)
cj.Draw()
#myText(0.20,0.24,1,myLumi,0.045)
#myText(0.20,0.19,1,myLumi2,0.045)
count=1
for tag in taggerList:
    if name1!="" and name1!=tag: continue
    myC2.cd()
    f1.cd()
    curve=f1.Get(tag+"---bc")
    f2.cd()
    if name2!="": tag=name2
    curve2=f2.Get(tag+"---bc")
    print curve
    if curve2==None: continue
    count+=1
    if count==5: count+=1
    if count==9: count=kOrange
    curve2.SetLineStyle(2)
    curve.SetLineColor(count)
    if (name1!="" and name2!=""): curve.SetLineColor(1)
    curve2.SetLineColor(count)
    curve.SetLineWidth(3)
    curve2.SetLineWidth(3)
    curve.Draw("C")
    curve2.Draw("C")

    curve3=None
    if name3!="":
        curve3=f2.Get(name3+"---bc")
        print curve3
        curve3.SetLineStyle(2)
        curve3.SetLineColor(4)
        curve3.SetLineWidth(3)
        curve3.Draw("C")
    
        
    myCx=TCanvas( "bVSc"+tag, "bVSc"+tag,800,800);
    pad_1=TPad("pad_1", "up", 0., 0.35, 1., 1.);
    pad_1.SetBottomMargin(0);
    pad_1.Draw();   
    pad_1.SetGridy()
    pad_1.SetGridx()
    pad_1.SetLogy()
    pad_2=TPad("pad_2", "down", 0.0, 0.00, 1.0, 0.35);
    pad_2.SetTopMargin(0);
    pad_2.SetBottomMargin(0.28);
    pad_2.Draw();
    pad_1.cd();
    
    cj.Draw()
    myText(0.20,0.24,1,myLumi,0.045)
    myText(0.20,0.19,1,myLumi2,0.045)
    curve.Draw("C")
    curve2.Draw("C")
    if name3!="": curve3.Draw("C")
    
    legend5=TLegend(0.55,0.7,0.92,0.90)
    legend5.SetTextFont(42)
    legend5.SetTextSize(0.038)
    legend5.SetFillColor(0)
    legend5.SetLineColor(0)
    legend5.SetFillStyle(0)
    legend5.SetBorderSize(0)
    if not (name1!="" and name2!=""):
        legend5.AddEntry(curve ,tag.replace("_","+")+" "+leg1,"L")
        legend5.AddEntry(curve2,tag.replace("_","+")+" "+leg2,"L")
    else:
        legend5.AddEntry(curve ,leg1+" : "+name1,"L")
        legend5.AddEntry(curve2,leg1+" : "+name2,"L")
        if name3!="": legend5.AddEntry(curve3,leg2+" : "+name3,"L")
    legend5.Draw("SAME")

    pad_2.cd()
    ratio=cj.Clone("ratio")
    ratio.GetYaxis().SetTitle("rel. diff.")
    ratio.GetYaxis().SetTitleOffset(1.1)
    ratio.GetXaxis().SetLabelSize(0.10)
    ratio.GetXaxis().SetTitleSize(0.10)
    ratio.GetYaxis().SetTitleOffset(0.7)
    ratio.GetYaxis().SetLabelSize(0.09)
    ratio.GetYaxis().SetTitleSize(0.09)
    ratio.SetMaximum(4.2)
    ratio.SetMinimum(-0.6)
    ratio.SetLineColor(kGray)
    ratio.Draw("HIST")

    ratC=TGraph()
    countPoint=0
    ratC2=TGraph()

    for bin in xrange(1,curve.GetN()+1):
        px1=Double(0)
        py1=Double(0)
        curve.GetPoint(bin-1,px1,py1)
        if bin%1500==0: print str(bin)+" ....: "+str(px1)
        if bin%2==0: continue
        if px1>0.999: continue
        if px1<cj.GetXaxis().GetXmin(): continue
        clo_py2=100000
        clo_px2=100000
        endV=0
        start=0
        clo_py2=curve2.Eval(px1)
        if clo_py2!=100000 :
            ratC.SetPoint(countPoint,px1,clo_py2/py1)
            #if "Errors" in ratC.ClassName(): ratC.SetPointError(countPoint,0,0)
            #ratC.SetPointError(countPoint,0,0)
            countPoint+=1
    for point in xrange(countPoint+1,ratC.GetN()+1):
        ratC.RemovePoint(point)
    ratC.SetLineWidth(3)
    ratC.SetLineColor(4)
    if (name1!="" and name2!=""): ratC.SetLineColor(2)
    line1=TLine(cj.GetXaxis().GetXmin(),1.,1.0,1.)
    line1.SetLineWidth(2)
    line1.SetLineStyle(2)
    line1.SetLineColor(922)
    line1.Draw("SAME")
    ratC.Draw("SAMEL")
    
    if name3!="":
        countPoint=0
        for bin in xrange(1,curve.GetN()+1):
            px1=Double(0)
            py1=Double(0)
            curve.GetPoint(bin-1,px1,py1)
            if bin%1500==0: print str(bin)+" ....: "+str(px1)
            if bin%3==0: continue
            if px1>0.999: continue
            if px1<cj.GetXaxis().GetXmin(): continue
            clo_py2=100000
            clo_px2=100000
            endV=0
            start=0
            clo_py2=curve2.Eval(px1)
            if clo_py2!=100000 :
            	ratC.SetPoint(countPoint,px1,clo_py2/py1)
            	#if "Errors" in ratC.ClassName(): ratC.SetPointError(countPoint,0,0)
                ratC.SetPointError(countPoint,0,0)
            	countPoint+=1
    
        for point in xrange(countPoint+1,ratC2.GetN()+1): ratC2.RemovePoint(point)
        ratC2.SetLineWidth(2)
        ratC2.SetLineColor(4)
        ratC2.Draw("C")
    
    myCx.Update()
    
    myCx.Print(odir+"/bVSc__"+tag+".pdf")
    #myCx.Print(odir+"/bVSc__"+tag+".eps")
    #myCx.Print(odir+"/bVSc__"+tag+".png")
    myCx.Print(odir+"/bVSc__"+tag+".C")

myC2.cd()
legend4.Draw()
#myText(0.20,0.24,1,myLumi,0.045)
#myText(0.20,0.19,1,myLumi2,0.045)
myC2.Update()
##myC2.Print(odir+"/bVScALL.pdf")

