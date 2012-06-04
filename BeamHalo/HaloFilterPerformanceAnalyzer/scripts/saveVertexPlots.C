/**
 * This macro pretties up the nVtx plots produced by
 * BeamHalo/HaloFilterPerformanceAnalyzer
 *
 * Marissa Rodenburg
 * May 2012
 */

{
  //
  // MET
  //
  /**
  TH1F* hall = _file0.Get("nVtx_METAll");
  TH1F* htight = _file0.Get("nVtx_METHaloTight");
  TH1F* hloose = _file0.Get("nVtx_METHaloLoose");
  
  hall->Sumw2();
  htight->Sumw2();
  hloose->Sumw2();
  
  htight->Divide(hall);
  hloose->Divide(hall);
  
  htight->Draw("e");
  htight->SetStats(0);
  htight->SetTitle(";nVertex;CSCTightHalo fraction");
  htight->SetMarkerStyle(20);
  htight->SetMarkerColor(kBlue+2);
  htight->SetLineColor(kBlue+2);
  htight->Draw("e");
  */
  hloose->SetAxisRange(0,27);
  hloose->Draw("e");
  hloose->SetStats(0);
  hloose->SetTitle(";nVertex;CSCLooseHalo fraction");
  hloose->SetMarkerStyle(20);
  hloose->SetMarkerColor(kBlue+2);
  hloose->SetLineColor(kBlue+2);
  hloose->Draw("e");
  

  //
  // MinimumBias(total rate)
  //
  
  /**
  TH1F* hall = _file0.Get("nVtx_All");
  TH1F* htight = _file0.Get("nVtx_HaloTight");
  TH1F* hloose = _file0.Get("nVtx_HaloLoose");
  
  hall->Sumw2();
  htight->Sumw2();
  hloose->Sumw2();
  
  htight.Divide(hall);
  hloose->Divide(hall);
  
  htight->Draw("e");
  htight->SetStats(0);
  htight->SetTitle(";nVertex;CSCTightHalo fraction");
  htight->SetMarkerStyle(20);
  htight->SetMarkerColor(kGreen+3); //(kRed+2);
  htight->SetLineColor(kGreen+3);   //(kRed+2);
  htight->Draw("e");
  
  hloose->SetAxisRange(0,27);
  hloose->Draw("e");
  hloose->SetStats(0);
  hloose->SetTitle(";nVertex;CSCLooseHalo fraction");
  hloose->SetMarkerStyle(20);
  hloose->SetMarkerColor(kGreen+3);  //(kRed+2);
  hloose->SetLineColor(kGreen+3);    //(kRed+2);
  hloose->Draw("e");
  */

}
