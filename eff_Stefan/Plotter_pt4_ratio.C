#include <string>
#include <iostream>
#include <sstream>

#include <dirent.h>

#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TSystem.h"

#include "AtlasStyle.C"
#include "AtlasUtils.C"

using namespace std;

TChain* myT_1;
TChain* myT_2;
TChain* myT_3;
TChain* myT_4;
TChain* myT_5;
TFile* outF;
string outputFolder;
string workpoint = "70";
/////////////////////////////////////////////////////
/// Plotting a lot of histograms from a given tree
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string getCut(string tagger, bool a8TeV) {
  if (tagger=="MV1")  return  "0.945487725";
  if (tagger=="MV1c") return  "0.779833333333";
  if (tagger=="MV2c00") return  "0.0308333333333";
  if (tagger=="MV2c10") return  "-0.00416666666667";
  if (tagger=="MV2c20" && workpoint == "60") return  "-0.0215";
  if (tagger=="MV2c20" && workpoint == "70") return  "-0.0215";
  if (tagger=="MV2c20" && workpoint == "77") return  "-0.0215";
  if (tagger=="MV2c20" && workpoint == "85") return  "-0.0215";
  if (tagger=="IP3D")      return  "2.007";
  if (tagger=="IP3D+SV1")  return  "4.3625";
  if (tagger=="MVb")       return  "-0.120991666667";
  if (tagger=="SV1")       return  "-97";
  if (tagger=="JetFitter") return  "-1.6125"; ///"-1.215  , -2.3 on Z
  
  cout << "NOT SUPPORTED!!! " << endl;
  return "0";
}

string getVariable(string tagger, bool a8TeV) {
  if (tagger=="MV1")    return "jet_mv1";
  if (tagger=="MV1c")   return "jet_mv1c";
  if (tagger=="MV2c00") return "jet_mv2c00";
  if (tagger=="MV2c10") return "jet_mv2c10";
  if (tagger=="MV2c20") return "jet_mv2c20";
  if (tagger=="MVb")    return "jet_mvb";
  if (tagger=="IP3D")   return  "jet_ip3d_llr";
  if (tagger=="IP3D+SV1")  return  "jet_sv1ip3d";
  if (tagger=="SV1")       return "jet_sv1_m";  //"jet_sv1_llr";
  if (tagger=="JetFitter") return  "jet_jf_m";
  cout << "NOT SUPPORTED!!! " << endl;
  return "0";
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D* GetHisto(string varName, string cutBase, 
	       string varLabel,string yLabel, 
	       int nBin, float Max, float Min,
	       bool normalize=false) {
  
  /// this is over ultra stupid but I am in a rush and I don't manage to get it to work otherwise
  TString tmpName=varName+cutBase;
  tmpName=tmpName.ReplaceAll(" ","").ReplaceAll("&","").ReplaceAll("(","").ReplaceAll(")","").ReplaceAll("=","").ReplaceAll(">","").ReplaceAll("<","").ReplaceAll("/1e3","").ReplaceAll(".","").ReplaceAll("+","");
  string theName=string(tmpName);

  TH1D* den  =new TH1D( theName.c_str(), theName.c_str(), nBin, Max, Min); den->Sumw2();
  string fullVar=varName+">>"+den->GetName();
  //  cout << "fullVar: " << fullVar << endl;
  //  cout << "cutBase: " << cutBase << endl;
  myT_1->Draw( fullVar.c_str(), cutBase.c_str(),"goff");//,1000000);
  //  cout << nBin << " , " << Max << " , " << Min << endl;
  //  cout << "Int: " << den->Integral() << endl; 

  den->SetBinContent(1, den->GetBinContent(0)+den->GetBinContent(1));
  den->SetBinError(1, sqrt(pow(den->GetBinError(0),2)+pow(den->GetBinError(1),2)));
  den->SetBinContent(0, 0.0);
  den->SetBinError(0, 0.0);
  den->SetBinContent(den->GetNbinsX(), den->GetBinContent(den->GetNbinsX()));
  den->SetBinError(den->GetNbinsX(), sqrt(pow(den->GetBinError(den->GetNbinsX()),2)+pow(den->GetBinError(den->GetNbinsX()+1),2)));
  den->SetBinContent(den->GetNbinsX()+1, 0.0);
  den->SetBinError(den->GetNbinsX()+1, 0.0);
  den->SetLineWidth(3);
  den->SetLineColor(2);
  den->SetMarkerStyle(20);
  den->SetMarkerSize(0.6);
  den->SetMarkerColor(2);
 
  if (normalize) {
    float maxV=den->GetMaximum()*1.1/den->Integral();
    den->Scale(1./den->Integral());
    den->SetMaximum(maxV);
  }
  return den;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D* GetHisto2(string varName, string cutBase, 
	       string varLabel,string yLabel, 
	       int nBin, float Max, float Min,
	       bool normalize=false) {
  
  /// this is over ultra stupid but I am in a rush and I don't manage to get it to work otherwise
  TString tmpName=varName+cutBase;
  tmpName=tmpName.ReplaceAll(" ","").ReplaceAll("&","").ReplaceAll("(","").ReplaceAll(")","").ReplaceAll("=","").ReplaceAll(">","").ReplaceAll("<","").ReplaceAll("/1e3","").ReplaceAll(".","").ReplaceAll("+","");
  string theName=string(tmpName);

  TH1D* den  =new TH1D( theName.c_str(), theName.c_str(), nBin, Max, Min); den->Sumw2();
  string fullVar=varName+">>"+den->GetName();
  //  cout << "fullVar: " << fullVar << endl;
  //  cout << "cutBase: " << cutBase << endl;
  myT_2->Draw( fullVar.c_str(), cutBase.c_str(),"goff");//,1000000);
  //  cout << nBin << " , " << Max << " , " << Min << endl;
  //  cout << "Int: " << den->Integral() << endl; 

  den->SetBinContent(1, den->GetBinContent(0)+den->GetBinContent(1));
  den->SetBinError(1, sqrt(pow(den->GetBinError(0),2)+pow(den->GetBinError(1),2)));
  den->SetBinContent(0, 0.0);
  den->SetBinError(0, 0.0);
  den->SetBinContent(den->GetNbinsX(), den->GetBinContent(den->GetNbinsX()));
  den->SetBinError(den->GetNbinsX(), sqrt(pow(den->GetBinError(den->GetNbinsX()),2)+pow(den->GetBinError(den->GetNbinsX()+1),2)));
  den->SetBinContent(den->GetNbinsX()+1, 0.0);
  den->SetBinError(den->GetNbinsX()+1, 0.0);
  den->SetLineWidth(3);
  den->SetLineColor(2);
  den->SetMarkerStyle(20);
  den->SetMarkerSize(0.6);
  den->SetMarkerColor(2);
 
  if (normalize) {
    float maxV=den->GetMaximum()*1.1/den->Integral();
    den->Scale(1./den->Integral());
    den->SetMaximum(maxV);
  }
  return den;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D* GetHisto3(string varName, string cutBase, 
	       string varLabel,string yLabel, 
	       int nBin, float Max, float Min,
	       bool normalize=false) {
  
  /// this is over ultra stupid but I am in a rush and I don't manage to get it to work otherwise
  TString tmpName=varName+cutBase;
  tmpName=tmpName.ReplaceAll(" ","").ReplaceAll("&","").ReplaceAll("(","").ReplaceAll(")","").ReplaceAll("=","").ReplaceAll(">","").ReplaceAll("<","").ReplaceAll("/1e3","").ReplaceAll(".","").ReplaceAll("+","");
  string theName=string(tmpName);

  TH1D* den  =new TH1D( theName.c_str(), theName.c_str(), nBin, Max, Min); den->Sumw2();
  string fullVar=varName+">>"+den->GetName();
  //  cout << "fullVar: " << fullVar << endl;
  //  cout << "cutBase: " << cutBase << endl;
  myT_3->Draw( fullVar.c_str(), cutBase.c_str(),"goff");//,1000000);
  //  cout << nBin << " , " << Max << " , " << Min << endl;
  //  cout << "Int: " << den->Integral() << endl; 

  den->SetBinContent(1, den->GetBinContent(0)+den->GetBinContent(1));
  den->SetBinError(1, sqrt(pow(den->GetBinError(0),2)+pow(den->GetBinError(1),2)));
  den->SetBinContent(0, 0.0);
  den->SetBinError(0, 0.0);
  den->SetBinContent(den->GetNbinsX(), den->GetBinContent(den->GetNbinsX()));
  den->SetBinError(den->GetNbinsX(), sqrt(pow(den->GetBinError(den->GetNbinsX()),2)+pow(den->GetBinError(den->GetNbinsX()+1),2)));
  den->SetBinContent(den->GetNbinsX()+1, 0.0);
  den->SetBinError(den->GetNbinsX()+1, 0.0);
  den->SetLineWidth(3);
  den->SetLineColor(2);
  den->SetMarkerStyle(20);
  den->SetMarkerSize(0.6);
  den->SetMarkerColor(2);
 
  if (normalize) {
    float maxV=den->GetMaximum()*1.1/den->Integral();
    den->Scale(1./den->Integral());
    den->SetMaximum(maxV);
  }
  return den;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D* GetHisto4(string varName, string cutBase, 
	       string varLabel,string yLabel, 
	       int nBin, float Max, float Min,
	       bool normalize=false) {
  
  /// this is over ultra stupid but I am in a rush and I don't manage to get it to work otherwise
  TString tmpName=varName+cutBase;
  tmpName=tmpName.ReplaceAll(" ","").ReplaceAll("&","").ReplaceAll("(","").ReplaceAll(")","").ReplaceAll("=","").ReplaceAll(">","").ReplaceAll("<","").ReplaceAll("/1e3","").ReplaceAll(".","").ReplaceAll("+","");
  string theName=string(tmpName);

  TH1D* den  =new TH1D( theName.c_str(), theName.c_str(), nBin, Max, Min); den->Sumw2();
  string fullVar=varName+">>"+den->GetName();
  //  cout << "fullVar: " << fullVar << endl;
  //  cout << "cutBase: " << cutBase << endl;
  myT_4->Draw( fullVar.c_str(), cutBase.c_str(),"goff");//,1000000);
  //  cout << nBin << " , " << Max << " , " << Min << endl;
  //  cout << "Int: " << den->Integral() << endl; 

  den->SetBinContent(1, den->GetBinContent(0)+den->GetBinContent(1));
  den->SetBinError(1, sqrt(pow(den->GetBinError(0),2)+pow(den->GetBinError(1),2)));
  den->SetBinContent(0, 0.0);
  den->SetBinError(0, 0.0);
  den->SetBinContent(den->GetNbinsX(), den->GetBinContent(den->GetNbinsX()));
  den->SetBinError(den->GetNbinsX(), sqrt(pow(den->GetBinError(den->GetNbinsX()),2)+pow(den->GetBinError(den->GetNbinsX()+1),2)));
  den->SetBinContent(den->GetNbinsX()+1, 0.0);
  den->SetBinError(den->GetNbinsX()+1, 0.0);
  den->SetLineWidth(3);
  den->SetLineColor(2);
  den->SetMarkerStyle(20);
  den->SetMarkerSize(0.6);
  den->SetMarkerColor(2);
 
  if (normalize) {
    float maxV=den->GetMaximum()*1.1/den->Integral();
    den->Scale(1./den->Integral());
    den->SetMaximum(maxV);
  }
  return den;
}


TH1D* GetHisto5(string varName, string cutBase, 
	       string varLabel,string yLabel, 
	       int nBin, float Max, float Min,
	       bool normalize=false) {
  
  /// this is over ultra stupid but I am in a rush and I don't manage to get it to work otherwise
  TString tmpName=varName+cutBase;
  tmpName=tmpName.ReplaceAll(" ","").ReplaceAll("&","").ReplaceAll("(","").ReplaceAll(")","").ReplaceAll("=","").ReplaceAll(">","").ReplaceAll("<","").ReplaceAll("/1e3","").ReplaceAll(".","").ReplaceAll("+","");
  string theName=string(tmpName);

  TH1D* den  =new TH1D( theName.c_str(), theName.c_str(), nBin, Max, Min); den->Sumw2();
  string fullVar=varName+">>"+den->GetName();
  //  cout << "fullVar: " << fullVar << endl;
  //  cout << "cutBase: " << cutBase << endl;
  myT_5->Draw( fullVar.c_str(), cutBase.c_str(),"goff");//,1000000);
  //  cout << nBin << " , " << Max << " , " << Min << endl;
  //  cout << "Int: " << den->Integral() << endl; 

  den->SetBinContent(1, den->GetBinContent(0)+den->GetBinContent(1));
  den->SetBinError(1, sqrt(pow(den->GetBinError(0),2)+pow(den->GetBinError(1),2)));
  den->SetBinContent(0, 0.0);
  den->SetBinError(0, 0.0);
  den->SetBinContent(den->GetNbinsX(), den->GetBinContent(den->GetNbinsX()));
  den->SetBinError(den->GetNbinsX(), sqrt(pow(den->GetBinError(den->GetNbinsX()),2)+pow(den->GetBinError(den->GetNbinsX()+1),2)));
  den->SetBinContent(den->GetNbinsX()+1, 0.0);
  den->SetBinError(den->GetNbinsX()+1, 0.0);
  den->SetLineWidth(3);
  den->SetLineColor(2);
  den->SetMarkerStyle(20);
  den->SetMarkerSize(0.6);
  den->SetMarkerColor(2);
 
  if (normalize) {
    float maxV=den->GetMaximum()*1.1/den->Integral();
    den->Scale(1./den->Integral());
    den->SetMaximum(maxV);
  }
  return den;
}



///////////////////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors* GetEfficiency2ITk(string varName, 
				 string cutBase, string effCut, 
				 string varLabel,string yLabel, 
				 int nBin, float Max, float Min, float &numC) {
  
  TH1D* num  =GetHisto(varName, cutBase+effCut, varLabel, "num", nBin, Max, Min);
  TH1D* den  =GetHisto(varName, cutBase, varLabel, "num", nBin, Max, Min);
  // cout << " num: " << num->Integral() << " , den: " << den->Integral() << endl;
  // numC=num->Integral();

  TH1D* num2  =GetHisto2(varName, cutBase+effCut, varLabel, "num", nBin, Max, Min);
  TH1D* den2  =GetHisto2(varName, cutBase, varLabel, "num", nBin, Max, Min);

  num->Divide(den);
  num2->Divide(den2);

  TH1D *newhist = (TH1D*) num->Clone();
  newhist -> Divide(num2);

  TGraphAsymmErrors* graphHisto= new TGraphAsymmErrors(newhist);
  graphHisto->SetLineWidth(3);
  graphHisto->SetLineColor(2);
  graphHisto->SetMarkerStyle(20);
  graphHisto->SetMarkerSize(0.6);
  graphHisto->SetMarkerColor(2);

  return graphHisto;
}



///////////////////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors* GetEfficiency3ITk(string varName, 
				 string cutBase, string effCut, 
				 string varLabel,string yLabel, 
				 int nBin, float Max, float Min, float &numC) {
  
  TH1D* num  =GetHisto(varName, cutBase+effCut, varLabel, "num", nBin, Max, Min);
  TH1D* den  =GetHisto(varName, cutBase, varLabel, "num", nBin, Max, Min);
  // cout << " num: " << num->Integral() << " , den: " << den->Integral() << endl;
  // numC=num->Integral();

  TH1D* num2  =GetHisto3(varName, cutBase+effCut, varLabel, "num", nBin, Max, Min);
  TH1D* den2  =GetHisto3(varName, cutBase, varLabel, "num", nBin, Max, Min);

  num->Divide(den);
  num2->Divide(den2);

  TH1D *newhist = (TH1D*) num->Clone();
  newhist -> Divide(num2);

  TGraphAsymmErrors* graphHisto= new TGraphAsymmErrors(newhist);
  graphHisto->SetLineWidth(3);
  graphHisto->SetLineColor(2);
  graphHisto->SetMarkerStyle(20);
  graphHisto->SetMarkerSize(0.6);
  graphHisto->SetMarkerColor(2);

  return graphHisto;
}




///////////////////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors* GetEfficiency4ITk(string varName, 
				 string cutBase, string effCut, 
				 string varLabel,string yLabel, 
				 int nBin, float Max, float Min, float &numC) {
  
  TH1D* num  =GetHisto(varName, cutBase+effCut, varLabel, "num", nBin, Max, Min);
  TH1D* den  =GetHisto(varName, cutBase, varLabel, "num", nBin, Max, Min);
  // cout << " num: " << num->Integral() << " , den: " << den->Integral() << endl;
  // numC=num->Integral();

  TH1D* num2  =GetHisto4(varName, cutBase+effCut, varLabel, "num", nBin, Max, Min);
  TH1D* den2  =GetHisto4(varName, cutBase, varLabel, "num", nBin, Max, Min);

  num->Divide(den);
  num2->Divide(den2);

  TH1D *newhist = (TH1D*) num->Clone();
  newhist -> Divide(num2);

  TGraphAsymmErrors* graphHisto= new TGraphAsymmErrors(newhist);
  graphHisto->SetLineWidth(3);
  graphHisto->SetLineColor(2);
  graphHisto->SetMarkerStyle(20);
  graphHisto->SetMarkerSize(0.6);
  graphHisto->SetMarkerColor(2);

  return graphHisto;
}

///////////////////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors* GetEfficiency(string varName, 
				 string cutBase, string effCut, 
				 string varLabel,string yLabel, 
				 int nBin, float Max, float Min, float &numC) {
  
  TH1D* num  =GetHisto(varName, cutBase+effCut, varLabel, "num", nBin, Max, Min);
  TH1D* den  =GetHisto(varName, cutBase, varLabel, "num", nBin, Max, Min);
  cout << " num: " << num->Integral() << " , den: " << den->Integral() << endl;
  numC=num->Integral();

  TGraphAsymmErrors* graphHisto= new TGraphAsymmErrors(num,den);
  graphHisto->SetLineWidth(3);
  graphHisto->SetLineColor(2);
  graphHisto->SetMarkerStyle(20);
  graphHisto->SetMarkerSize(0.6);
  graphHisto->SetMarkerColor(2);

  return graphHisto;
}

///////////////////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors* GetEfficiency2(string varName, 
				 string cutBase, string effCut, 
				 string varLabel,string yLabel, 
				 int nBin, float Max, float Min, float &numC) {
  
  TH1D* num  =GetHisto2(varName, cutBase+effCut, varLabel, "num", nBin, Max, Min);
  TH1D* den  =GetHisto2(varName, cutBase, varLabel, "num", nBin, Max, Min);
  cout << " num: " << num->Integral() << " , den: " << den->Integral() << endl;
  numC=num->Integral();

  TGraphAsymmErrors* graphHisto= new TGraphAsymmErrors(num,den);
  graphHisto->SetLineWidth(3);
  graphHisto->SetLineColor(2);
  graphHisto->SetMarkerStyle(20);
  graphHisto->SetMarkerSize(0.6);
  graphHisto->SetMarkerColor(2);

  return graphHisto;
}
///////////////////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors* GetEfficiency3(string varName, 
				 string cutBase, string effCut, 
				 string varLabel,string yLabel, 
				 int nBin, float Max, float Min, float &numC) {
  
  TH1D* num  =GetHisto3(varName, cutBase+effCut, varLabel, "num", nBin, Max, Min);
  TH1D* den  =GetHisto3(varName, cutBase, varLabel, "num", nBin, Max, Min);
  cout << " num: " << num->Integral() << " , den: " << den->Integral() << endl;
  numC=num->Integral();

  TGraphAsymmErrors* graphHisto= new TGraphAsymmErrors(num,den);
  graphHisto->SetLineWidth(3);
  graphHisto->SetLineColor(2);
  graphHisto->SetMarkerStyle(20);
  graphHisto->SetMarkerSize(0.6);
  graphHisto->SetMarkerColor(2);

  return graphHisto;
}
///////////////////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors* GetEfficiency4(string varName, 
				 string cutBase, string effCut, 
				 string varLabel,string yLabel, 
				 int nBin, float Max, float Min, float &numC) {
  
  TH1D* num  =GetHisto4(varName, cutBase+effCut, varLabel, "num", nBin, Max, Min);
  TH1D* den  =GetHisto4(varName, cutBase, varLabel, "num", nBin, Max, Min);
  cout << " num: " << num->Integral() << " , den: " << den->Integral() << endl;
  numC=num->Integral();

  TGraphAsymmErrors* graphHisto= new TGraphAsymmErrors(num,den);
  graphHisto->SetLineWidth(3);
  graphHisto->SetLineColor(2);
  graphHisto->SetMarkerStyle(20);
  graphHisto->SetMarkerSize(0.6);
  graphHisto->SetMarkerColor(2);

  return graphHisto;
}

///////////////////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors* GetEfficiency5(string varName, 
				 string cutBase, string effCut, 
				 string varLabel,string yLabel, 
				 int nBin, float Max, float Min, float &numC) {
  
  TH1D* num  =GetHisto5(varName, cutBase+effCut, varLabel, "num", nBin, Max, Min);
  TH1D* den  =GetHisto5(varName, cutBase, varLabel, "num", nBin, Max, Min);
  cout << " num: " << num->Integral() << " , den: " << den->Integral() << endl;
  numC=num->Integral();

  TGraphAsymmErrors* graphHisto= new TGraphAsymmErrors(num,den);
  graphHisto->SetLineWidth(3);
  graphHisto->SetLineColor(2);
  graphHisto->SetMarkerStyle(20);
  graphHisto->SetMarkerSize(0.6);
  graphHisto->SetMarkerColor(2);

  return graphHisto;
}



void GetComparison(string file1, string cutBase, string effCut,
		   string cut_flav,
		   string varName, string varLabel,string yLabel, 
		   int nBin, float Max, float Min, string tagger) {
  TCanvas* can=new TCanvas( (varLabel+"  "+yLabel).c_str(), (varLabel+"  "+yLabel).c_str(), 900, 900);
  
  TPad *pad_1 = new TPad("pad_1", "up", 0., 0.35, 1., 1.);
  pad_1->SetBottomMargin(0.02);
  pad_1->Draw();   

  TPad *pad_2 = new TPad("pad_2", "down", 0.0, 0.00, 1.0, 0.35);
  pad_2->SetTopMargin(0.0);
  pad_2->SetBottomMargin(0.28);
  pad_2->Draw();
  
  pad_1->cd();
  float numC;
  
  string tmpVarName=varName;
 
  TH1D* mainH=new TH1D("histo_den", "2", nBin, Max, Min);
  mainH->Reset();
  mainH->SetTitle( (";"+varLabel+";"+yLabel+";").c_str() );
  mainH->GetXaxis()->SetLabelOffset(0.1);
  //  if(varName=="jet_pt/1e3") mainH->SetMaximum(0.02);
  //  else mainH->SetMaximum(0.2); // chiara: maximum set here
  //else mainH->SetMaximum(0.014);
  
  // mainH->SetAxisRange(0.0,4.5,"X");
  mainH->SetMinimum(0.0);
  mainH->SetDirectory(0);
  //////// mainH->GetYaxis()->SetNdivisions(505);

    if (myT_1==0) {
    myT_1=new TChain("bTag_AntiKt4EMTopoJets");
    cout << " OPENING FILE: " << file1 << endl;
    if ( file1.find("root")!=string::npos )  {
      myT_1->Add( file1.c_str() );
    } else {
      cout << "Input is a directory: ging fancy: " << endl;
      DIR*     dir;
      dirent*  pdir;
      dir = opendir( file1.c_str() );     // open current directory
      while (pdir = readdir(dir))  {
	string foldName=pdir->d_name;
	if (foldName.find("mc")==string::npos && foldName.find("valid")==string::npos) continue;
	cout << pdir->d_name << endl;
	DIR*     dir2;
	dirent*  pdir2;
	dir2 = opendir( (file1+"/"+foldName).c_str() );     // open current directory
	while (pdir2 = readdir(dir2))  {
	  string fName=pdir2->d_name;
	  //cout << " fName: " << fName << endl;
	  if (fName.find("root")==string::npos) continue;
	  //cout << fName << endl;
	  myT_1->Add( (file1+"/"+foldName+"/"+fName).c_str() );
	}
      }
    }
  }

  cout << "TOTAL number of events is: " << myT_1->GetEntries() << endl;

  string file2="/eos/user/c/crizzi/HGTD/btagging/input_eff_plot/Initial_trkEffselfTag/file.root";
  myT_2=new TChain("bTag_AntiKt4EMTopoJets");
  cout << " OPENING FILE: " << file2 << endl;
    if ( file2.find("root")!=string::npos )  {
      myT_2->Add( file2.c_str() );
    } else {
      cout << "Input is a directory: ging fancy: " << endl;
      DIR*     dir;
      dirent*  pdir;
      dir = opendir( file2.c_str() );     // open current directory
      while (pdir = readdir(dir))  {
	string foldName=pdir->d_name;
	if (foldName.find("mc")==string::npos && foldName.find("valid")==string::npos) continue;
	cout << pdir->d_name << endl;
	DIR*     dir2;
	dirent*  pdir2;
	dir2 = opendir( (file2+"/"+foldName).c_str() );     // open current directory
	while (pdir2 = readdir(dir2))  {
	  string fName=pdir2->d_name;
	  //cout << " fName: " << fName << endl;
	  if (fName.find("root")==string::npos) continue;
	  //cout << fName << endl;
	  myT_2->Add( (file2+"/"+foldName+"/"+fName).c_str() );
	}
      }
    }

  cout << "TOTAL number of events 2 is : " << myT_2->GetEntries() << endl;

  string file3="/eos/user/c/crizzi/HGTD/btagging/input_eff_plot/Intermediate_pre_Replacement_trkEffselfTag/file.root";
  myT_3=new TChain("bTag_AntiKt4EMTopoJets");
  cout << " OPENING FILE: " << file3 << endl;
  myT_3->Add( file3.c_str() );
 

  cout << "TOTAL number of events 3 is : " << myT_3->GetEntries() << endl;

  string file4="/eos/user/c/crizzi/HGTD/btagging/input_eff_plot/Intermediate_post_Replacement_trkEffselfTag/file.root";
  myT_4=new TChain("bTag_AntiKt4EMTopoJets");
  cout << " OPENING FILE: " << file4 << endl;
  myT_4->Add( file4.c_str() );

  string file5="/eos/user/c/crizzi/HGTD/btagging/input_eff_plot/Final_trkEffselfTag/file.root";
  myT_5=new TChain("bTag_AntiKt4EMTopoJets");
  cout << " OPENING FILE: " << file5 << endl;
  myT_5->Add( file5.c_str() );


  cout << "TOTAL number of events 4 is : " << myT_4->GetEntries() << endl;


 
  TGraphAsymmErrors* gra1=GetEfficiency( varName, 
					 (cutBase+cut_flav),  effCut, 
					 varLabel, yLabel,
					 nBin, Max, Min, numC);

 
  TGraphAsymmErrors* gra2=GetEfficiency2( varName, 
					 (cutBase+cut_flav), effCut, 
					 varLabel, yLabel,
					 nBin, Max, Min, numC);
 
  TGraphAsymmErrors* gra3=GetEfficiency3( varName, 
					 cutBase+cut_flav, effCut, 
					 varLabel, yLabel,
					 nBin, Max, Min, numC);

  TGraphAsymmErrors* gra4=GetEfficiency4( varName, 
					 cutBase+cut_flav, effCut, 
					 varLabel, yLabel,
					 nBin, Max, Min, numC);

  TGraphAsymmErrors* gra5=GetEfficiency5( varName, 
					 cutBase+cut_flav, effCut, 
					 varLabel, yLabel,
					 nBin, Max, Min, numC);


  double hmax = 0;
  int nPoints = gra1->GetN(); // number of points in your TGraph
  for(int i=0; i < nPoints; ++i) {
    double x,y;
    gra1->GetPoint(i, x, y);
    if(y>hmax)
      hmax=y;
  }
  nPoints = gra2->GetN(); // number of points in your TGraph
  for(int i=0; i < nPoints; ++i) {
    double x,y;
    gra2->GetPoint(i, x, y);
    if(y>hmax)
      hmax=y;
  }
  nPoints = gra3->GetN(); // number of points in your TGraph
  for(int i=0; i < nPoints; ++i) {
    double x,y;
    gra3->GetPoint(i, x, y);
    if(y>hmax)
      hmax=y;
  }
  nPoints = gra4->GetN(); // number of points in your TGraph
  for(int i=0; i < nPoints; ++i) {
    double x,y;
    gra4->GetPoint(i, x, y);
    if(y>hmax)
      hmax=y;
  }

  mainH->SetMaximum(1.5*hmax);
  mainH->Draw("HIST");  

  gra1->SetLineColor(1);
  gra1->SetMarkerColor(1);
  gra1->Draw("LX");
  gra2->SetLineColor(kGreen+2);
  gra2->SetMarkerColor(kGreen+2);
  gra2->Draw("LX");
  gra3->SetLineColor(kRed+1);
  gra3->SetMarkerColor(kRed+1);
  gra3->Draw("LX");
  gra4->SetLineColor(kGray+1);
  gra4->SetMarkerColor(kGray+1);
  gra4->Draw("LX");
 



  TLine* myL=new TLine(Min,1.0,Max,1.0);
  myL->SetLineStyle(2);
  myL->SetLineColor(kGray);
  myL->SetLineWidth(2);
  myL->Draw("SAME");
  if ( varName.find("Lxy")!=string::npos) {
    TLine* myL5=new TLine(33,0.0,33,1.05);
    myL5->SetLineStyle(2);
    myL5->SetLineColor(1);
    myL5->SetLineWidth(2);
    myL5->Draw("SAME");
    TLine* myL3=new TLine(50.5,0.0,50.5,1.05);
    myL3->SetLineStyle(2);
    myL3->SetLineColor(1);
    myL3->SetLineWidth(2);
    myL3->Draw("SAME");
    TLine* myL4=new TLine(88.5,0.0,88.5,1.05);
    myL4->SetLineStyle(2);
    myL4->SetLineColor(1);
    myL4->SetLineWidth(2);
    myL4->Draw("SAME");
  }

  TLegend* legend4=new TLegend(0.23,0.63,0.45,0.80);
  legend4->SetTextFont(42);
  legend4->SetTextSize(0.04);
  legend4->SetFillColor(0);
  legend4->SetLineColor(0);
  legend4->SetFillStyle(0);
  legend4->SetBorderSize(0);

  //legend4->AddEntry(gra1 ,"ITk"   ,"l");
  //legend4->AddEntry(gra2 ,"ITk+HGTD (Full #eta)","l");
  //legend4->AddEntry(gra3 ,"ITk+HGTD (2.0 < |#eta| < 4.0)","l");
  //legend4->AddEntry(gra4 ,"ITk+HGTD (2.4 < |#eta| < 4.0)","l");

  legend4->AddEntry(gra1 ,"ITk-only"   ,"l");
  legend4->AddEntry(gra2 ,"Initial","l");
  legend4->AddEntry(gra3 ,"Int. pre-repl.","l");
  legend4->AddEntry(gra4 ,"Int. post-repl.","l");
  legend4->AddEntry(gra5 ,"Final","l");
  
  legend4->Draw("SAME");
 

  //ATLAS_LABEL(0.22,0.88,1,0.035); 

  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize(0.035);
  l.SetTextColor(1);
  l.DrawLatex(0.22,0.88,"#bf{#it{ATLAS}} Simulation Internal");

  TLatex l2; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l2.SetNDC();
  l2.SetTextFont(42);
  l2.SetTextSize(0.035);
  l2.SetTextColor(1);
  l2.DrawLatex(0.22,0.83,"t#bar{t} simulation, jet p_{T} > 20 GeV");

  TLatex l3; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l3.SetNDC();
  l3.SetTextFont(42);
  l3.SetTextSize(0.035);
  l3.SetTextColor(1);
  l3.DrawLatex(0.65,0.83,TString(tagger)+TString(" 70% fixed cut WP"));

  TLatex l4; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l4.SetNDC();
  l4.SetTextFont(42);
  l4.SetTextSize(0.035);
  l4.SetTextColor(1);
  l4.DrawLatex(0.65,0.88,"HGTD (2.4 < |#eta| < 4.0)");


  // RATIO
  pad_2->cd();
 
  TH1D* mainH2=new TH1D("histo_den2", "2", nBin, Max, Min);
  mainH2->Reset();
  mainH2->SetMaximum(4.2);
  mainH2->SetTitle( (";"+varLabel+";Ratio to ITk;").c_str() );
  mainH2->SetTitleSize(0.09);
  mainH2->GetXaxis()->SetLabelSize(0.1);
  mainH2->GetYaxis()->SetLabelSize(0.09);
  mainH2->GetYaxis()->SetTitleOffset(0.85);
  mainH2->GetYaxis()->SetTitleSize(0.08);
  
  // mainH->SetAxisRange(0.0,4.5,"X");
  mainH2->SetMinimum(0.0);
  mainH2->Draw("HIST");
  mainH2->SetDirectory(0);
  mainH2->GetYaxis()->SetNdivisions(505);

  //// Make ratio

  TGraphAsymmErrors* gradiv1=GetEfficiency2ITk( varName, 
					 cutBase+cut_flav, effCut, 
					 varLabel, yLabel,
					 nBin, Max, Min, numC);
					
  gradiv1->SetLineColor(kGreen+2);
  gradiv1->SetMarkerColor(kGreen+2);
  gradiv1->Draw("LX");

  TGraphAsymmErrors* gradiv2=GetEfficiency3ITk( varName, 
					 cutBase+cut_flav, effCut, 
					 varLabel, yLabel,
					 nBin, Max, Min, numC);

					
  gradiv2->SetLineColor(kRed+1);
  gradiv2->SetMarkerColor(kRed+1);
  gradiv2->Draw("LX");

  TGraphAsymmErrors* gradiv3=GetEfficiency4ITk( varName, 
					 cutBase+cut_flav, effCut, 
					 varLabel, yLabel,
				       	nBin, Max, Min, numC);


  gradiv3->SetLineColor(kGray+1);
  gradiv3->SetMarkerColor(kGray+1);
  gradiv3->Draw("LX  ");

  TLine* myL2=new TLine(Min,1.0,Max,1.0);
  myL2->SetLineStyle(2);
  myL2->SetLineColor(kGray);
  myL2->SetLineWidth(2);
  myL2->Draw("SAME");
  TLine* myL3=new TLine(Min,5.0,Max,5.0);
  myL3->SetLineStyle(2);
  myL3->SetLineColor(kGray);
  myL3->SetLineWidth(2);
  // myL3->Draw("SAME");

  varName=tmpVarName+"_"+tagger;
  TString varToPrint="Eff__"+varName+"__"+yLabel+"_"+".eps";
  varToPrint=varToPrint.ReplaceAll("/1e3","").ReplaceAll("abs(","").ReplaceAll(")","").ReplaceAll("=","").ReplaceAll(">","_").ReplaceAll("eff.","").ReplaceAll("%","").ReplaceAll("@","_").ReplaceAll(" ","").ReplaceAll("+","");
  
  outF->cd();
  TString baseName="Base__"+varName;
  baseName= baseName.ReplaceAll("/1e3","").ReplaceAll("abs(","").ReplaceAll(")","").ReplaceAll("=","").ReplaceAll(">","_").ReplaceAll("eff.","").ReplaceAll("%","").ReplaceAll("@","_").ReplaceAll(" ","").ReplaceAll("+","");
  
  TString plotName= varToPrint.ReplaceAll(".eps","");
  if (numC!=0) outF->WriteObject(mainH,baseName);

  TString bEff="Eff_b__"+plotName;
  TString cEff="Eff_c__"+plotName;
  TString lEff="Eff_l__"+plotName;
  TString lEff2="Eff2_l__"+plotName;
  if (numC!=0) {
  gra1->SetName(bEff);
  gra1->Write();
  gra2->SetName(cEff);
  gra2->Write();
  gra3->SetName(lEff);
  gra3->Write();
  gra4->SetName(lEff2);
  gra4->Write();

  gradiv1->Write();
  gradiv2->Write();
  gradiv3->Write();

  }
  //outF->WriteObject(gra1,bEff);
  //outF->WriteObject(gra2,cEff);
  //outF->WriteObject(gra3,lEff);
  
  TString varToPrint1=outputFolder+"/"+varToPrint+"_3.eps";
  TString varToPrint2=outputFolder+"/"+varToPrint+"_3.png";
  TString varToPrint3=outputFolder+"/"+varToPrint+"_3.C";
  TString varToPrint4=outputFolder+"/"+varToPrint+"_3.pdf";
  can->Print( varToPrint1 );
  can->Print( varToPrint2 );
  can->Print( varToPrint3 );
  can->Print( varToPrint4 );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GetComparisonSimple(string file1, string cutBase,  
			 string Cut1 , string Cut2 , string Cut3 , 
			 string varName, string varLabel,string yLabel, 
			 int nBin, float Max, float Min, bool log=false) {

  TCanvas* can=new TCanvas( (varLabel+"  "+yLabel).c_str(), (varLabel+"  "+yLabel).c_str(), 900, 900);
  if (log) can->SetLogy();
  
  if (myT_1==0) {
    myT_1=new TChain("bTag");
    cout << " OPENING FILE: " << file1 << endl;
    if ( file1.find("root")!=string::npos )  {
      myT_1->Add( file1.c_str() );
    } else {
      cout << "Input is a directory: ging fancy: " << endl;
      DIR*     dir;
      dirent*  pdir;
      dir = opendir( file1.c_str() );     // open current directory
      int number_of_words=0;
      int text_length = 30;
      char filename[300];
      while (pdir = readdir(dir))  {
	cout << pdir->d_name << endl;
      }
    }
  }
    
  cout << "TOTAL number of events is: " << myT_1->GetEntries() << endl;
  cout << cutBase+Cut1 << endl;

  TH1D* gra1=GetHisto( varName,
		       (cutBase+Cut1+" ").c_str(), 
		       varLabel, yLabel,
		       nBin, Max, Min,true);
  gra1->SetLineColor(2);
  gra1->SetLineWidth(3);
  gra1->SetMarkerColor(2);
  gra1->SetTitle( (";"+varLabel+";"+yLabel+";").c_str() );
  gra1->SetMinimum(0.001);
  TH1D* gra1b=GetHisto( varName, 
			(cutBase+Cut1+" && bH_Lxy>33").c_str(), 
			varLabel, yLabel,
			nBin, Max, Min,true);
  gra1b->SetLineColor(6);
  gra1b->SetMarkerColor(6);


  TH1D* gra2=GetHisto( varName, 
		       cutBase+Cut2, 
		       varLabel, yLabel,
		       nBin, Max, Min,true);
  gra2->SetLineColor(8);
  gra2->SetMarkerColor(8);
 
  
  TH1D* gra3=GetHisto( varName, 
		       cutBase+Cut3, 
		       varLabel, yLabel,
		       nBin, Max, Min,true);
  gra3->SetLineColor(4);
  gra3->SetMarkerColor(4);

  gra1->Draw("HIST");
  gra2->Draw("SAMEHIST");
  gra3->Draw("SAMEHIST");
  gra1b->Draw("SAMEHIST");
  gra1->Draw("SAMEHIST");

  TLegend* legend4=new TLegend(0.67,0.68,0.920,0.93);
  legend4->SetTextFont(42);
  legend4->SetTextSize(0.04);
  legend4->SetFillColor(0);
  legend4->SetLineColor(0);
  legend4->SetFillStyle(0);
  legend4->SetBorderSize(0);
  legend4->AddEntry(gra1 ,"b-jets","l");
  legend4->AddEntry(gra1b,"b-jets (after IBL)","l");
  legend4->AddEntry(gra2 ,"c-jets"   ,"l");
  legend4->AddEntry(gra3 ,"light-jet","l");
  legend4->Draw("SAME");

  TString varToPrint="Var__"+varName+".eps";
  varToPrint=varToPrint.ReplaceAll("/1e3","");
  can->Print( varToPrint );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PrintTagger(string tagger, string file1, string cut_flav, string yLabel)  {
  /// even more ugly
  
  string CutBase="";
  CutBase=" jet_pt>20e3 && jet_truthMatch==1 && jet_isPU==0 && abs(PVz-truth_PVz)<0.1 && abs(jet_eta)<4";  
  
  // MV1: quite detailed info
  string effCut=" && "+getVariable(tagger, false)+">"+getCut(tagger, false)+" ";
 
  GetComparison(file1,CutBase, effCut,
		cut_flav,
	        "abs(jet_eta)", "jet |#eta|", yLabel,  
		19, -0.1, 3.7, tagger);

}
  


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Plotter_pt4(const char* infile,
		 const char* outfolder, string wkpoint) {
  gStyle->SetOptStat(0);
  SetAtlasStyle();
  
  outputFolder=outfolder;
  workpoint=wkpoint;
  
  string file1=infile;
  gSystem->Exec( ("mkdir -p "+outputFolder).c_str());
  outF=new TFile( (outputFolder+"/effPlots.root").c_str(),"RECREATE");
  cout << "Created file: " << outF->GetName() << endl;

  string Cut1=" "; 
  string Cut2=" "; 
  string Cut3=" ";  
  Cut1=" && jet_LabDr_HadF==5 "; 
  Cut2=" && jet_LabDr_HadF==4 "; 
  Cut3=" && (jet_LabDr_HadF!=4 && jet_LabDr_HadF!=5 && jet_LabDr_HadF!=15) && jet_dRminToB>0.8 && jet_dRminToC>0.8 && jet_dRminToT>0.8"; 

  string yLabel="Light-jet mis-tagging efficiency";
  PrintTagger("MV1",file1, Cut3, yLabel);
  PrintTagger("IP3D+SV1",file1, Cut3, yLabel);
  PrintTagger("IP3D",file1, Cut3, yLabel);
  PrintTagger("SV1",file1, Cut3, yLabel);


  yLabel="B-jet tagging efficiency";
  PrintTagger("MV1",file1, Cut1, yLabel);
  PrintTagger("IP3D+SV1",file1, Cut1, yLabel);
  PrintTagger("IP3D",file1, Cut1, yLabel);
  PrintTagger("SV1",file1, Cut1, yLabel);

  //PrintTagger("MV1c",file1);
  //PrintTagger("MV2c00",file1);
  //PrintTagger("MV2c10",file1);
  //PrintTagger("MV2c20",file1);
  //PrintTagger("MVb",file1);
  //  PrintTagger("IP3D+SV1",file1);
  //  PrintTagger("IP3D",file1);
  //  PrintTagger("SV1",file1);

  // PrintTagger("JetFitter",file1);

  outF->Close();
}

