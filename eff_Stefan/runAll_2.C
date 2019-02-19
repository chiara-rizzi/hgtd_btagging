{
  gROOT->ProcessLine(".L Plotter_pt2.C++g");
  cout << "##############################################################################################" << endl;

  //gROOT->ProcessLine("Plotter_pt2(\"/tmp/guindon/user.guindon.noHGTD.root\",\"NO_HGTD\",\"70\")");
   gROOT->ProcessLine("Plotter_pt2(\"input/HGTD/file.root\",\"HGTD_Forward_Checks\",\"70\")");


  gROOT->ProcessLine(".q");
}
