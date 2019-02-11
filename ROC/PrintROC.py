#import argparse
from glob import glob
import os, sys
from array import *
import ROOT
from ROOT import *#TCanvas,TLegend,gROOT
import math
import time
from math import sqrt

ROOT.gROOT.LoadMacro("AtlasStyle.C")
ROOT.gROOT.LoadMacro("AtlasUtils.C")
SetAtlasStyle()


tmpVale=sys.argv

def Helper():
    print " "
    print " Usage: python PrintROC.py <outFolder>  <fileName> or <folderName>"
    print "  ... outFolder MUST containt LC or EM "
    print " "
    sys.exit(1)

if len(sys.argv) !=3: Helper()

    
########################################################################################
#class MyParser(argparse.ArgumentParser):
#    def error(self, message):
#        sys.stderr.write('error: %s\n' % message)
#        self.print_help()
#        sys.exit(2)
#parser = argparse.ArgumentParser(description='Debugging info in yields and shapes.')
#parser.add_argument("-i","--idir",
#                    help="Directory containing input files",
#                    dest="idir",
#                    default="")
#parser.add_argument("-o","--odir",
#                    help="Directory containing all the outputs",
#                    dest="odir",
#                    default="testOoutput")
#parser.add_argument("-f","--file",
#                    help="file to run on",
#                    dest="file",
#                    default="a.root")

#options= parser.parse_args()
#print parser.parse_args()

odir = sys.argv[1]
if ".root" in odir:
    print " ..... outputfolder is a file .... please correct"
    sys.exit(1)

second=sys.argv[2]
fileList=[]
if ".root" in second:
    test=TFile(second)
    if test==None:
        print " ..... input file: "+test+" ... NOT FOUND"
        sys.exit(1)
    fileList.append(second)
else:
    print "ls -la "+second+"/*/*.root*"
    fileList=glob(second+"/*/*.root*")
    
if len(fileList)==0:
    print " ..... no files selected .... please check"
    sys.exit(1)

#if "EM" not in odir and "LC" not in odir:
#    Helper()
#    sys.exit(1)

gROOT.SetBatch(True) 
########################################################################################

taggers=[]

is8TeV=False
isXAOD=False
if "8TeV" in odir:
    is8TeV=True
if "XAOD" in odir:
    isXAOD=True

if is8TeV and isXAOD:
    taggers=[
        ["IP3D"    , "ip3d_llr", -12. ,   30,  3000,   8 ],
        ["SV1"     , "sv1_llr" ,  -4. ,   13,  3000,   6 ],
        ["IP3D+SV1", "sv1ip3d" , -10. ,   35,  3000, 797 ],
        ]
elif isXAOD:
    taggers=[
        ["MV1"     , "mv1"     ,   0.0 ,  0.9945 ,2000, 1 ],   #20000
##        ["MV1c"    , "mv1c"    ,   0   ,  1      ,  3000, 3 ],
#        ["MV2c00"  , "mv2c00"  ,  -0.5 ,  0.5    ,  3000, 2 ],
        ["MV2c10"  , "mv2c10"  ,  -0.5 ,  0.5    ,  3000, 4 ],
        ["MV2c20"  , "mv2c20"  ,  -0.5 ,  0.5    ,  3000, 7 ],
#        ["MV2c00"  , "mv2c00"  ,  -1.01 ,  1.01    ,  1000, 2 ],
#        ["MV2c10"  , "mv2c10"  ,  -1.01 ,  1.01    ,  1000, 4 ],
#        ["MV2c20"  , "mv2c20"  ,  -1.01 ,  1.01    ,  1000, 7 ],
        ["IP3D"    , "ip3d_llr", -12.  ,   30    ,  1000, 8 ],
        ["IP2D"    , "ip2d_llr", -12.  ,   30    ,  1000, 10 ],
        ["SV1"     , "sv1_llr" ,  -4.  ,   13    ,  1000, 6 ],
        ["IP3D+SV1", "sv1ip3d" , -10.  ,   35    ,  3000, 797 ],
#        ["MVb"     , "mvb"     ,  -1.01,  1.01   ,  1000, 920 ],
        ["JetFitter"     , "jf_llr"     ,  -15,  10    ,  1000, 40 ]
        ]
else:
     taggers=[
              ["MV1"     , "mv1"     ,   0.0  ,  0.9945  , 2000, 1 ],   #20000
   #           ["MV1c"    , "mv1c"    ,   0  ,  1  ,  2000, 3 ],
   #           ["MV2c00"  , "mv2c00"  ,  -0.5,  0.5,  2000, 2 ],
              ["MV2c10"  , "mv2c10"  ,  -1.01,  1.01,  1000, 4 ],
              ["MV2c20"  , "mv2c20"  ,  -1.01,  1.01,  1000, 7 ],
              ["IP3D"    , "ip3d_llr"    , -12.   ,   30,  2000, 8 ],
              ["IP2D"    , "ip2d_llr",     -12.  ,   30    ,  1000, 10 ],
              ["SV1"     , "sv1_llr"     ,  -4.   ,   13,  2000, 6 ],
              ["JetFitter"     , "jf_llr"     ,  -15,  10    ,  1000, 40 ]
              #["MVb"   , "mvb"       ,  -1.05 ,  0.8,  2000, 920 ],
              ]
effThreshold=0.70


def GetHisto(tag, intree, val):
    tmpH=TH1F(tag[1]+str(val),tag[1]+str(tag[5]),tag[4],tag[2],tag[3])
    tmpH.Sumw2()
    var="jet_"+tag[1]+">>"+tmpH.GetName()
    print val

#     // track confirmation
  #  rpt = jet_sumtrkS_pt / jet_pt 
#    print rp

#      if (fabs(jeteta[j])<2.4) {
#	if (jetpt[j]<40e3) Rptcut = 0.12;
#	else if (jetpt[j]<50e3) Rptcut = 0.13;
#	else if (jetpt[j]<70e3) Rptcut = 0.15;
#	else if (jetpt[j]<100e3) Rptcut = 0.11;
 #     } else if (fabs(jeteta[j])<3.2) {
#	if (jetpt[j]<40e3) Rptcut = 0.16;
#	else if (jetpt[j]<50e3) Rptcut = 0.17;
#	else if (jetpt[j]<70e3) Rptcut = 0.14;
#	else if (jetpt[j]<100e3) Rptcut = 0.12;
#      } else {
#	if (jetpt[j]<40e3) Rptcut = 0.17;
#	else if (jetpt[j]<50e3) Rptcut = 0.18;
#	else if (jetpt[j]<70e3) Rptcut = 0.18;
#	else if (jetpt[j]<100e3) Rptcut = 0.15;
 #     }

    if val<4: 
        cut=" jet_dRminToB>0.8 && jet_dRminToC>0.8 && jet_dRminToT>0.8  && jet_LabDr_HadF=="+str(val)+"&& jet_pt>20e3 && jet_isPU==0 && jet_truthMatch==1 && fabs(PVz-truth_PVz)<0.1 && fabs(jet_eta)>2.4"
    else: 
        cut=" jet_LabDr_HadF=="+str(val)+"&& jet_pt>20e3 && jet_truthMatch==1 && jet_isPU==0 && fabs(PVz-truth_PVz)<0.1 && fabs(jet_eta)>2.4" 

    if"jet_pt<40e3": cut+=" && jet_sumtrkS_pt/jet_pt > 0.15"
    
    print cut

    ####cut=" abs(jet_eta)<2.5 && jet_pt>25e3 &&  jet_truthMatch==1 && jet_truthflav=="+str(val)
    ##if not is8T: cut="jet_truthflav=="+str(val)+" && jet_pt>25e3  &&  (jet_JVT>0.2 || jet_pt>50e3 || abs(jet_eta)>2.4) "
    ##             cut="jet_LabDr_HadF=="+str(val)+" && jet_pt>25e3  &&  (jet_JVT>0.2 || jet_pt>50e3 || abs(jet_eta)>2.4) "
    ##             cut="jet_GhostL_HadF=="+str(val)+" && jet_pt>25e3  && jet_truthMatch==1 "
    ##             cut="jet_GhostL_HadF=="+str(val)+" && jet_pt>25e3  && jet_truthMatch==1 && jet_ip3d_ntrk>=2 "
    ##             cut="jet_trueFlav=="+str(val)+" && jet_pt>25e3 && jet_truthmatched==1 "




    intree.Draw( var, cut,"goof") 
    tmpH.SetBinContent(1,tmpH.GetBinContent(1)+tmpH.GetBinContent(0))
    tmpH.SetBinError(1,sqrt(pow(tmpH.GetBinError(1),2)+pow(tmpH.GetBinError(0),2)))
    tmpH.SetBinContent(0,0.0)
    tmpH.SetBinError(0,0.0)
    tmpH.SetBinContent(tmpH.GetNbinsX(),tmpH.GetBinContent(tmpH.GetNbinsX())+tmpH.GetBinContent(tmpH.GetNbinsX()+1))
    tmpH.SetBinError(tmpH.GetNbinsX(),sqrt(pow(tmpH.GetBinError(tmpH.GetNbinsX()),2)+pow(tmpH.GetBinError(tmpH.GetNbinsX()+1),2)))
    tmpH.SetBinContent(tmpH.GetNbinsX()+1,0.0)
    tmpH.SetBinError(tmpH.GetNbinsX()+1,0.0)
    if val!=4:
        c = TCanvas("_Discr","_Discr",800,600)    
        tmpH.GetXaxis().SetTitle(tag[0])
        tmpH.Draw("E")
        c.Update()
        c.Print(odir+"/test.eps")
        time.sleep(0.5)
    tmpInt=tmpH.Integral(-10,tmpH.GetNbinsX()+10)
    if tmpH.Integral():
        tmpH.Scale(1./tmpH.Integral(-10,tmpH.GetNbinsX()+10))
    else: print tmpH.Integral(-10,tmpH.GetNbinsX()+10)
    return tmpH,tmpInt

    
def GetROC(tag,intree, bVSlight):
    intSig=-1
    intBkgd=-1
    hsig,intSig=GetHisto(tag, intree,5)
    hbkgd=None
    if not bVSlight: hbkgd,intBkgd=GetHisto(tag, intree,4)
    else           : hbkgd,intBkgd=GetHisto(tag, intree,0)

    if bVSlight:
        found=False
        for bin in xrange(1,hsig.GetNbinsX()+2):
            partInt=hsig.Integral(bin,hsig.GetNbinsX()+10)
            #print str(bin)+"    "+str(hsig.GetBinCenter(bin))+"    "+str(partInt)
            if partInt<effThreshold and not found:
                print " CUT= "+str(hsig.GetBinCenter(bin))+" has eff: "+str(partInt)
                found=True
                
    myROC=TGraphErrors();##hsig.GetNbinsX()-2 );
    maxRej=1
    count=-1
    for bin in xrange(2,hsig.GetNbinsX()):
        sigEff =hsig.Integral(bin,hsig.GetNbinsX()+10)
        bkgdEff=hbkgd.Integral(bin,hsig.GetNbinsX()+10)
        ##if bVSlight: print str(bkgdEff)+"   "+str(sigEff)+" CUT: "+str(hsig.GetBinCenter(bin))
        if bkgdEff!=0 and sigEff!=0 and sigEff<0.99:
            ex=0#sqrt( sigEff*(1-sigEff)/intSig )
            ey=0#sqrt( bkgdEff*(1-bkgdEff)/intBkgd )
            ##if bVSlight: print str(bkgdEff)+"   "+str(sigEff)
            count+=1
            myROC.SetPoint(count,sigEff,1/bkgdEff);
            myROC.SetPointError(count, ex, ey/(bkgdEff*bkgdEff) )
            if 1/bkgdEff>maxRej: maxRej=1/bkgdEff
    myROC.SetLineWidth(3)
    return myROC,maxRej
  
#########################################################################################
gSystem.Exec("mkdir -p "+odir)

intree=None
intree=TChain("bTag_AntiKt4EMTopoJets")
##if "EM" in odir: intree=TChain("bTag_AntiKt4EMTopoJets")
##else           : intree=TChain("bTag_AntiKt4LCTopoJets")
for file in fileList:
    intree.Add(file)

print "    "
print "TOTAL NUMBER OF ENTRIES IS: "+str(intree.GetEntries())
print "    "

countT=-1
light=TH1F("b VS light","b VS light",100,0.3,1);
light.SetTitle(";b-jet efficiency;light-jet rejection;")
lightCurve=[]

cj=TH1F("b VS c","b VS c",100,0.3,1);
cj.SetTitle(";b-jet efficiency;c-jet rejection;")
cCurve=[]

for tag in taggers:
    countT+=1
    print "<<<<<<<  "+tag[0]+"   >>>>>>>>>>>>>>"
    curve,Rej=GetROC(tag,intree,True)
    curve.SetLineColor(tag[5])
    lightCurve.append(curve)
    if Rej*2>light.GetMaximum(): light.SetMaximum(Rej*2)

    curve2,Rej=GetROC(tag,intree,False)
    curve2.SetLineColor(tag[5])
    cCurve.append(curve2)
    if Rej*2>cj.GetMaximum(): cj.SetMaximum(Rej*2)

ofile=TFile(odir+"/output.root","RECREATE")        

myC=TCanvas( "bVSl", "bVSl",900,900);
myC.SetLogy()
myC.SetGridy()
myC.SetGridx()
light.SetMinimum(1)
light.SetMaximum(1e5)
light.Draw()
myLumi= "t#bar{t} simulation,"
if is8TeV and isXAOD: myLumi+=" 8TeV, rel19"
elif isXAOD: myLumi+=" 13 TeV, rel19"
else       : myLumi+=" 8TeV, rel17"
myLumi2= "jet p_{T}>25 GeV, |#eta|<2.5"
#myText(0.20,0.24,1,myLumi,0.045)
#myText(0.20,0.19,1,myLumi2,0.045)
legend4=TLegend(0.70,0.60,0.92,0.92)
legend4.SetTextFont(42)
legend4.SetTextSize(0.04)
legend4.SetFillColor(0)
legend4.SetLineColor(0)
legend4.SetFillStyle(0)
legend4.SetBorderSize(0)
count=-1
for curve in lightCurve:
    curve.Draw("C")
    count+=1
    legend4.AddEntry(curve,taggers[count][0],"L")
    #ofile.WriteObject(curve,taggers[count][0].replace("+","_")+"---bl")
    curve.SetName(taggers[count][0].replace("+","_")+"---bl")
    curve.Write()
legend4.Draw()
myC.Update()
myC.Print(odir+"/bVSlight.eps")
myC.Print(odir+"/bVSlight.png")
myC.Print(odir+"/bVSlight.C")

myC2=TCanvas( "cVSl", "cVSl",900,900);
myC2.SetLogy()
myC2.SetGridy()
myC2.SetGridx()
cj.SetMinimum(1)
cj.SetMaximum(1e3)
cj.Draw()
count=-1
for curve in cCurve:
    curve.Draw("C")
    count+=1
    ofile.WriteObject(curve,taggers[count][0].replace("+","_")+"---bc")
legend4.Draw()
#myText(0.20,0.24,1,myLumi,0.045)
#myText(0.20,0.19,1,myLumi2,0.045)
myC2.Update()
myC2.Print(odir+"/cVSlight.eps")
myC2.Print(odir+"/cVSlight.png")
myC2.Print(odir+"/cVSlight.C")

ofile.Close()

##time.sleep(111111)
