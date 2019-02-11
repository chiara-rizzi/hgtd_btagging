{
  //gROOT->ProcessLine(".L Plotter_pt3.C++g");
  cout << "##############################################################################################" << endl;

  //gROOT->ProcessLine("Plotter_pt2(\"/tmp/guindon/user.guindon.noHGTD.root\",\"NO_HGTD\",\"70\")");
  //gROOT->ProcessLine("Plotter_pt3(\"/eos/user/g/guindon/user.guindon.user.guindon.noHGTD_fixedJetFitter.root\",\"HGTD_Forward_Checks\",\"70\")");

  //gROOT->ProcessLine("Plotter_pt3(\"/tmp/guindon/user.guindon.user.guindon.noHGTD_fixedJetFitter_2.root\",\"HGTD_Forward_Checks\",\"70\")");

  // gROOT->ProcessLine(".L Plotter_pt3_vale_ratio.C++g");
  //gROOT->ProcessLine("Plotter_pt3(\"/eos/user/g/guindon/user.guindon.noHGTD_trk.root\",\"HGTD_Forward_Checks_New_ITk\",0.7)");

  gROOT->SetBatch(kTRUE);
  gROOT->ProcessLine(".L Plotter_pt4_ratio.C++g");
  gROOT->ProcessLine("Plotter_pt4(\"input/HGTD/file.root\",\"HGTD_Forward_Checks_Step2p2\",\"70\")");

  //gROOT->ProcessLine(".L Plotter_pt4_ratio_eff.C++g");
  //gROOT->ProcessLine("Plotter_pt4(\"/eos/user/g/guindon/user.guindon.user.guindon.noHGTD_fixedJetFitter.root\",\"HGTD_Forward_Checks_ITk\",\"70\")");

  //gROOT->ProcessLine(".L Plotter_pt4_ratio.C++g");
  //gROOT->ProcessLine("Plotter_pt4(\"/eos/user/g/guindon/user.guindon.user.guindon.noHGTD_fixedJetFitter.root\",\"HGTD_Forward_Checks_ITk\",\"70\")");

  gROOT->ProcessLine(".q");
}
