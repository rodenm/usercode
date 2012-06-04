// -*- C++ -*-
//
// Package:    HaloFilterPerformanceAnalyzer
// Class:      HaloFilterPerformanceAnalyzer
// 
/**\class HaloFilterPerformanceAnalyzer HaloFilterPerformanceAnalyzer.cc BeamHalo/HaloFilterPerformanceAnalyzer/src/HaloFilterPerformanceAnalyzer.cc

 Description: This analyzer handles the basic functions needed to test the beam halo filter.
       1. Performs a basic analysis of the efficiency of loose & tight
       2. Produces histograms showing efficiency of each flag
       3. Repeats #1 with a requirement of 50GeV MET opposite in phi to a CSCSegment

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marissa Rodenburg
//         Created:  Wed May 16 15:22:57 CDT 2012
// $Id: HaloFilterPerformanceAnalyzer.cc,v 1.2 2012/05/22 03:12:44 rodenm Exp $
//
//

//
// TODO:
// - Add code to calculate halo rate by fill
// - Break out some of this to an EDFilter
// - Have event list write out to a text file. Make this crab-compatible.



// system include files
#include <memory>
#include <sstream>

// root files
#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLegend.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/METReco/interface/CSCHaloData.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"

#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

//
// class declaration
//

class HaloFilterPerformanceAnalyzer : public edm::EDAnalyzer {
public:
  explicit HaloFilterPerformanceAnalyzer(const edm::ParameterSet&);
  ~HaloFilterPerformanceAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // helper functions
  virtual void analyzeHaloWithMET(const edm::Event&, const edm::EventSetup&, 
				  bool, const bool, const bool, const int);

      // ----------member data ---------------------------

  TH1F hHaloWithMET_;
  TH1F hNotHaloFlags_;
  TH1F hLooseHaloFlags_;
  TH1F hTightHaloFlags_;
  TH1F hHaloSummary_;
  TH1F hHaloMET_;
  TH1F hMETPhi_;
  TH1F hCSCSegmentPhi_;
  TH1F hMET_;
  TH1F hL1HaloTrigger_;
  TH1F hLooseHaloMET_;
  TH1F hTightHaloMET_;
  TH1F hCrudeDphiCSCSegmentMET_;
  TH1F hDphiCSCSegmentMET_;
  TH1F hDphiProbeCSCSegmentMET_;

  TH1F hnVtxAll_;
  TH1F hnVtxHaloLoose_;
  TH1F hnVtxHaloTight_;
  TH1F hnVtxMETAll_;
  TH1F hnVtxMETHaloLoose_;
  TH1F hnVtxMETHaloTight_;


  // TODO: this should be a cfg file parameter
  bool printEventInfo_;
  bool printMETEventList_; // If false, only print regular event lists
  float minimumMET_;    // in GeV
  std::stringstream cscEvents_;
  std::stringstream cscLooseEvents_;
  std::stringstream cscTightEvents_;
  std::stringstream cscNeitherEvents_;

  std::stringstream metCSCEvents_;
  std::stringstream metCSCLooseEvents_;
  std::stringstream metCSCTightEvents_;
  std::stringstream metCSCNeitherEvents_;
  std::string outputFileName_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HaloFilterPerformanceAnalyzer::HaloFilterPerformanceAnalyzer(const edm::ParameterSet& iConfig) :
  printEventInfo_(iConfig.getParameter<bool>("printEventInfo")),
  minimumMET_(iConfig.getParameter<double>("minimumMET")),
  outputFileName_(iConfig.getParameter<std::string>("outputFile"))
{
  //now do what ever initialization is needed
  //printEventInfo_ = false;
  //minimumMET_ = 50.0; 
  printMETEventList_ = false;
  
  // TODO: only initialize this when analyzeMET will be called
  // metCSCEvents_;
}


HaloFilterPerformanceAnalyzer::~HaloFilterPerformanceAnalyzer() {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void HaloFilterPerformanceAnalyzer::beginJob() {
  hHaloWithMET_ = TH1F("HaloWithMET", "", 5, 0.5, 5.5);
  hHaloWithMET_.GetXaxis()->SetBinLabel(1,"MET>50GeV");
  hHaloWithMET_.GetXaxis()->SetBinLabel(2,"MET>50GeV && Opposite CSCSegment");
  hHaloWithMET_.GetXaxis()->SetBinLabel(2,"MET>50GeV && Opposite NC-CSCSegment");
  hHaloWithMET_.GetXaxis()->SetBinLabel(3,"MET>50GeV && Opposite NC-CSCSegment && CSCLoose");
  hHaloWithMET_.GetXaxis()->SetBinLabel(4,"MET>50GeV && Opposite NC-CSCSegment && CSCTight");
  hHaloWithMET_.SetTitle("MET requirements & halo filter results");

  hNotHaloFlags_   = TH1F("NotHaloFlags", "", 8, 0.5, 8.5);
  hNotHaloFlags_.GetXaxis()->SetBinLabel(1,"All");
  hNotHaloFlags_.GetXaxis()->SetBinLabel(2,"CSCLooseHaloId || CSCTightHaloId");
  hNotHaloFlags_.GetXaxis()->SetBinLabel(3,"NumberOfHaloTriggers > 0");
  hNotHaloFlags_.GetXaxis()->SetBinLabel(4,"NumberOfHaloTracks > 0");
  hNotHaloFlags_.GetXaxis()->SetBinLabel(5,"NOutOfTimeHits() > 10");
  hNotHaloFlags_.GetXaxis()->SetBinLabel(6,"NFlatHaloSegments() > 2");
  hNotHaloFlags_.GetXaxis()->SetBinLabel(7,"SegmentsInBothEndcaps()==1");
  hNotHaloFlags_.GetXaxis()->SetBinLabel(8,"NTracksSmalldT() > 0");

  hLooseHaloFlags_   = TH1F("LooseHaloFlags", "", 8, 0.5, 8.5);
  hLooseHaloFlags_.GetXaxis()->SetBinLabel(1,"All");
  hLooseHaloFlags_.GetXaxis()->SetBinLabel(2,"CSCLooseHaloId==1");
  hLooseHaloFlags_.GetXaxis()->SetBinLabel(3,"NumberOfHaloTriggers > 0");
  hLooseHaloFlags_.GetXaxis()->SetBinLabel(4,"NumberOfHaloTracks > 0");
  hLooseHaloFlags_.GetXaxis()->SetBinLabel(5,"NOutOfTimeHits() > 10");
  hLooseHaloFlags_.GetXaxis()->SetBinLabel(6,"NFlatHaloSegments() > 2");
  hLooseHaloFlags_.GetXaxis()->SetBinLabel(7,"SegmentsInBothEndcaps()==1");
  hLooseHaloFlags_.GetXaxis()->SetBinLabel(8,"NTracksSmalldT() > 0");

  hTightHaloFlags_   = TH1F("TightHaloFlags", "", 8, 0.5, 8.5);
  hTightHaloFlags_.GetXaxis()->SetBinLabel(1,"All");
  hTightHaloFlags_.GetXaxis()->SetBinLabel(2,"CSCTightHaloId==1");
  hTightHaloFlags_.GetXaxis()->SetBinLabel(3,"NumberOfHaloTriggers > 0");
  hTightHaloFlags_.GetXaxis()->SetBinLabel(4,"NumberOfHaloTracks > 0");
  hTightHaloFlags_.GetXaxis()->SetBinLabel(5,"NOutOfTimeHits() > 10");
  hTightHaloFlags_.GetXaxis()->SetBinLabel(6,"NFlatHaloSegments() > 3");
  hTightHaloFlags_.GetXaxis()->SetBinLabel(7,"SegmentsInBothEndcaps()==1");
  hTightHaloFlags_.GetXaxis()->SetBinLabel(8,"NTracksSmalldT() > 0");

  hHaloSummary_ = TH1F("HaloSummary", "", 3, 0.5, 3.5);
  hHaloSummary_.GetXaxis()->SetBinLabel(1,"All");
  hHaloSummary_.GetXaxis()->SetBinLabel(2,"CSCLooseHaloId==1");
  hHaloSummary_.GetXaxis()->SetBinLabel(3,"CSCTightHaloId==1");

  hLooseHaloMET_ = TH1F("LooseHaloMET", "", 300, 0, 300);
  hTightHaloMET_ = TH1F("TightHaloMET", "", 300, 0, 300);

  hMETPhi_ = TH1F("METPhi", "", 72 ,-TMath::Pi() , TMath::Pi());
  hCSCSegmentPhi_ = TH1F("CSCSegmentPhi", "", 72 ,-TMath::Pi() , TMath::Pi());
  hMET_ = TH1F("AllMET", "", 200, 0.0, 200.0);

  const double dphiBins[12] = {0., .3 ,.6 ,.9 ,1.2 ,1.5, 1.8, 2.1, 2.4, 2.7, 3.0, TMath::Pi()};
  hCrudeDphiCSCSegmentMET_ = TH1F("crudeDphi_CSCSegmentMET", "", 11, dphiBins);
  hDphiCSCSegmentMET_ = TH1F("Dphi_CSCSegmentMET", "" , 36, 0, TMath::Pi());
  hDphiProbeCSCSegmentMET_ = TH1F("DphiProbe_CSCSegmentMET", "" , 36, 0, TMath::Pi());

  hL1HaloTrigger_ = TH1F("L1HaloTrigger", "", 2, 0.5, 2.5);
  hL1HaloTrigger_.GetXaxis()->SetBinLabel(1,"No Accept");
  hL1HaloTrigger_.GetXaxis()->SetBinLabel(2,"Accept");

  hnVtxAll_ = TH1F("nVtx_All", "", 30, 0.5, 30.5);
  hnVtxHaloLoose_ = TH1F("nVtx_HaloLoose", "", 30, 0.5, 30.5);
  hnVtxHaloTight_ = TH1F("nVtx_HaloTight", "", 30, 0.5, 30.5);
  hnVtxMETAll_ = TH1F("nVtx_METAll", "", 30, 0.5, 30.5);
  hnVtxMETHaloLoose_ = TH1F("nVtx_METHaloLoose", "", 30, 0.5, 30.5);
  hnVtxMETHaloTight_ = TH1F("nVtx_METHaloTight", "", 30, 0.5, 30.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void HaloFilterPerformanceAnalyzer::endJob() {

  // TODO: this should be a cfg file parameter
  //TFile *outputFile = new TFile("haloFilterPerformance.root","RECREATE");
  TFile *outputFile = new TFile(outputFileName_.c_str(),"RECREATE");

  outputFile->cd();
  hHaloWithMET_.Write();
  hNotHaloFlags_.Write();
  hLooseHaloFlags_.Write();
  hTightHaloFlags_.Write();
  hHaloSummary_.Write();
  
  hLooseHaloMET_.Write();
  hTightHaloMET_.Write();

  hnVtxAll_.Write();
  hnVtxHaloLoose_.Write();
  hnVtxHaloTight_.Write();
  hnVtxMETAll_.Write();
  hnVtxMETHaloLoose_.Write();
  hnVtxMETHaloTight_.Write();

  hMETPhi_.Write();
  hCSCSegmentPhi_.Write();
  hCrudeDphiCSCSegmentMET_.Write();
  hDphiCSCSegmentMET_.Write();
  hDphiProbeCSCSegmentMET_.Write();
  hMET_.Write();
  hL1HaloTrigger_.Write();


  // TODO: make this contingent on a bool that is passed in from the config file
  if (printMETEventList_) {
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Events with MET > " << minimumMET_ << " and opposite NC-CSCSegment" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << metCSCEvents_.str() << std::endl << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Events with MET > " << minimumMET_ << " and opposite NC-CSCSegment && CSCLoose" 
	      << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << metCSCLooseEvents_.str() << std::endl << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Events with MET > " << minimumMET_ << " and opposite NC-CSCSegment && CSCTight" 
	      << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << metCSCTightEvents_.str() << std::endl << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Events with MET > " << minimumMET_ << " and opposite NC-CSCSegment && NotHalo" 
	      << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << metCSCNeitherEvents_.str() << std::endl << std::endl;
  } else {
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "All events with CSCSegment" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << cscEvents_.str() << std::endl << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Events with CSCLoose" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << cscLooseEvents_.str() << std::endl << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Events with CSCTight" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << cscTightEvents_.str() << std::endl << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Events with NotHalo" << std::endl;
    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << cscNeitherEvents_.str() << std::endl << std::endl;
  }    

  std::cout << "---------------Beam Halo Summary-----------------" << std::endl;
  std::cout << "All events:     " << hHaloSummary_.GetBinContent(1) << std::endl;
  std::cout << "Loose events:   " << hHaloSummary_.GetBinContent(2) << std::endl;
  std::cout << "Tight events:   " << hHaloSummary_.GetBinContent(3) << std::endl;
  std::cout << "L1Halo Trigger: " << hL1HaloTrigger_.GetBinContent(2) << std::endl <<std::endl;

  std::cout << "---------------MET Filter Summary-----------------" << std::endl;
  std::cout << "MET>50GeV:                                       " << hHaloWithMET_.GetBinContent(1) << std::endl;
  std::cout << "MET>50GeV && Opposite CSCSegment:                " << hHaloWithMET_.GetBinContent(2) << std::endl;
  std::cout << "MET>50GeV && Opposite NC-CSCSegment:             " << hHaloWithMET_.GetBinContent(3) << std::endl;
  std::cout << "MET>50GeV && Opposite NC-CSCSegment && CSCLoose: " << hHaloWithMET_.GetBinContent(4) << std::endl;
  std::cout << "MET>50GeV && Opposite NC-CSCSegment && CSCTight: " << hHaloWithMET_.GetBinContent(5) << std::endl;
  std::cout << std::endl;
  // Use a crude sideband subtraction to estimate the non-halo contamination 
  // of the "MET>50GeV && Opposite NC-CSCSegment" signal region. This crude because
  // Dphi isn't exactly flat (I'm assuming it is here); there's a bit of a peak
  // near Dphi = 0. However, it's small for pfmet so I'm rocking the easy method.
  //
  // Documentation: http://lphe.epfl.ch/publications/diplomas/fd.master.pdf page 17
  //
  double nTotal = hCrudeDphiCSCSegmentMET_.Integral(1,11);
  double nBackground = hCrudeDphiCSCSegmentMET_.Integral(2,10);
  double nBackgroundInSignal = nBackground*(TMath::Pi()-3.0)/2.7;
  std::cout << "---------------MET Filter Summary-----------------" << std::endl;
  std::cout << "N_total    = " << nTotal << std::endl;
  std::cout << "N_bkgd     = " << nBackground << std::endl;
  std::cout << "Sigma_sig  = " << TMath::Pi()-3.0 << std::endl;
  std::cout << "Sigma_bkgd = " << 2.7 << std::endl;
  std::cout << "N_bkgd,sig = " << nBackgroundInSignal << std::endl;
  std::cout << "Events in signal region: " << hCrudeDphiCSCSegmentMET_.GetBinContent(11) << std::endl;
  std::cout << "Non-halo contamination rate: " << 100.0*nBackgroundInSignal/hCrudeDphiCSCSegmentMET_.GetBinContent(11)
	    << "%" << std::endl << std::endl;
}

// ------------ method called for each event  ------------
void HaloFilterPerformanceAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // TODO: make this a class variable
  // TODO: make these labels parameters in the cfg file
  edm::Handle<reco::BeamHaloSummary> beamHaloSummaryHandle;
  iEvent.getByLabel("BeamHaloSummary",beamHaloSummaryHandle);
  
  // CSC Specific Halo Data
  edm::Handle<reco::CSCHaloData> cscHaloDataHandle;
  iEvent.getByLabel("CSCHaloData", cscHaloDataHandle);
  
  edm::Handle< reco::CaloMETCollection > TheCaloMET;
  iEvent.getByLabel("met", TheCaloMET); 
  const reco::CaloMETCollection *calometcol = TheCaloMET.product();
  const reco::CaloMET *calomet = &(calometcol->front());

  edm::Handle< reco::PFMETCollection > ThePFMET;
  iEvent.getByLabel("pfMet", ThePFMET); 
  const reco::PFMETCollection *pfmetcol = ThePFMET.product();
  const reco::PFMET *pfmet = &(pfmetcol->front());

  edm::Handle<reco::VertexCollection> recoVertices;
  iEvent.getByLabel("offlinePrimaryVertices", recoVertices);
  
  bool eventPrinted = false;
  bool cscLooseId = false;
  bool cscTightId = false;
  int nVtx = 0;
  
  if( beamHaloSummaryHandle.isValid() && cscHaloDataHandle.isValid() && recoVertices.isValid() ) {
    const reco::BeamHaloSummary haloSummary = (*beamHaloSummaryHandle.product() );
    const reco::CSCHaloData cscHaloData = (*cscHaloDataHandle.product() );
    
    hNotHaloFlags_.SetBinContent(1, hNotHaloFlags_.GetBinContent(1)+1);
    hHaloSummary_.SetBinContent(1, hHaloSummary_.GetBinContent(1)+1);
    //hLooseHaloFlags_.SetBinContent(1, hLooseHaloFlags_.GetBinContent(1)+1); // bad for scale
    //hTightHaloFlags_.SetBinContent(1, hTightHaloFlags_.GetBinContent(1)+1); // bad for scale

    cscEvents_ << iEvent.id().run() << ":" 
	       << iEvent.luminosityBlock() << ":" 
               << iEvent.id().event()
               << std::endl;
 
    // Count halo triggers - curiosity 
    if ( cscHaloData.NumberOfHaloTriggers() )
      hL1HaloTrigger_.SetBinContent(2, hL1HaloTrigger_.GetBinContent(2)+1);
    else
      hL1HaloTrigger_.SetBinContent(1, hL1HaloTrigger_.GetBinContent(1)+1);

    // Count # of vertices
    reco::VertexCollection::const_iterator it;
    for(it = recoVertices->begin(); it!=recoVertices->end(); ++it) {
      if (!it->isFake())
	nVtx++;
    }
    hnVtxAll_.Fill(nVtx);
  
    // Loose Halo
    if (haloSummary.CSCLooseHaloId()){
      cscLooseId = true;
      hHaloSummary_.SetBinContent(2, hHaloSummary_.GetBinContent(2)+1);
      hLooseHaloFlags_.SetBinContent(2, hLooseHaloFlags_.GetBinContent(2)+1);
      //hLooseHaloMET_.Fill(calomet->pt());
      hLooseHaloMET_.Fill(pfmet->pt());
      hnVtxHaloLoose_.Fill(nVtx);
      cscLooseEvents_ << iEvent.id().run() << ":" 
		 << iEvent.luminosityBlock() << ":" 
		 << iEvent.id().event()
		 << std::endl;
  
     
      if ( cscHaloData.NumberOfHaloTriggers() )
	hLooseHaloFlags_.SetBinContent(3, hLooseHaloFlags_.GetBinContent(3)+1);
      
      if ( cscHaloData.NumberOfHaloTracks() ) 
	hLooseHaloFlags_.SetBinContent(4, hLooseHaloFlags_.GetBinContent(4)+1);

      if ( cscHaloData.NOutOfTimeHits() > 10 )
	hLooseHaloFlags_.SetBinContent(5, hLooseHaloFlags_.GetBinContent(5)+1);

      if ( cscHaloData.NFlatHaloSegments() > 2 )
	hLooseHaloFlags_.SetBinContent(6, hLooseHaloFlags_.GetBinContent(6)+1);
      
      if ( cscHaloData.GetSegmentsInBothEndcaps() )
	hLooseHaloFlags_.SetBinContent(7, hLooseHaloFlags_.GetBinContent(7)+1);

      if ( cscHaloData.NTracksSmalldT() )
	hLooseHaloFlags_.SetBinContent(8, hLooseHaloFlags_.GetBinContent(8)+1);
    }

    // Tight Halo
    if (haloSummary.CSCTightHaloId()) {
      cscTightId = true;
      hHaloSummary_.SetBinContent(3, hHaloSummary_.GetBinContent(3)+1);
      hTightHaloFlags_.SetBinContent(2, hTightHaloFlags_.GetBinContent(2)+1);
      //hTightHaloMET_.Fill(calomet->pt());
      hTightHaloMET_.Fill(pfmet->pt());
      hnVtxHaloTight_.Fill(nVtx);
      cscTightEvents_ << iEvent.id().run() << ":" 
		      << iEvent.luminosityBlock() << ":" 
		      << iEvent.id().event()
		      << std::endl;
  
      if ( cscHaloData.NumberOfHaloTriggers() )
	hTightHaloFlags_.SetBinContent(3, hTightHaloFlags_.GetBinContent(3)+1);
      
      if ( cscHaloData.NumberOfHaloTracks() ) 
	hTightHaloFlags_.SetBinContent(4, hTightHaloFlags_.GetBinContent(4)+1);

      if ( cscHaloData.NOutOfTimeHits() > 10 )
	hTightHaloFlags_.SetBinContent(5, hTightHaloFlags_.GetBinContent(5)+1);

      if ( cscHaloData.NFlatHaloSegments() > 3 )
	hTightHaloFlags_.SetBinContent(6, hTightHaloFlags_.GetBinContent(6)+1);
      
      if ( cscHaloData.GetSegmentsInBothEndcaps() )
	hTightHaloFlags_.SetBinContent(7, hTightHaloFlags_.GetBinContent(7)+1);

      if ( cscHaloData.NTracksSmalldT() )
	hTightHaloFlags_.SetBinContent(8, hTightHaloFlags_.GetBinContent(8)+1);
    }

    // Not Halo
    if ( !haloSummary.CSCTightHaloId() && !haloSummary.CSCLooseHaloId() ) {
        cscNeitherEvents_ << iEvent.id().run() << ":" 
			  << iEvent.luminosityBlock() << ":" 
			  << iEvent.id().event()
			  << std::endl;

      if ( haloSummary.CSCTightHaloId() || haloSummary.CSCLooseHaloId() )  // this should be redundant
	hNotHaloFlags_.SetBinContent(2, hLooseHaloFlags_.GetBinContent(2)+1);
 
      if ( cscHaloData.NumberOfHaloTriggers() )
	hNotHaloFlags_.SetBinContent(3, hNotHaloFlags_.GetBinContent(3)+1);
	
      if ( cscHaloData.NumberOfHaloTracks() ) 
	hNotHaloFlags_.SetBinContent(4, hNotHaloFlags_.GetBinContent(4)+1);

      if ( cscHaloData.NOutOfTimeHits() > 10 )
	hNotHaloFlags_.SetBinContent(5, hNotHaloFlags_.GetBinContent(5)+1);

      if ( cscHaloData.NFlatHaloSegments() > 2 )
	hNotHaloFlags_.SetBinContent(6, hNotHaloFlags_.GetBinContent(6)+1);
      
      if ( cscHaloData.GetSegmentsInBothEndcaps() )
	hNotHaloFlags_.SetBinContent(7, hNotHaloFlags_.GetBinContent(7)+1);

      if ( cscHaloData.NTracksSmalldT() )
	hNotHaloFlags_.SetBinContent(8, hNotHaloFlags_.GetBinContent(8)+1);
    }

    // If printEventInfo_ == true, print out event number and results of loose and tight haloId
    if (printEventInfo_ && (haloSummary.CSCTightHaloId() || haloSummary.CSCLooseHaloId()) ) {
      eventPrinted = true;
      std::cout << "================= " 
		<<  iEvent.id().run() << ":" <<iEvent.luminosityBlock() << ":" << iEvent.id().event()
		<< " =================" << std::endl;
      std::cout << "CSCLooseHalo: " << haloSummary.CSCLooseHaloId() << std::endl;
      std::cout << "CSCTightHalo: " << haloSummary.CSCTightHaloId() << std::endl;
    }
  } else { 
    if (!beamHaloSummaryHandle.isValid())
      edm::LogError("HaloFilterPerformanceAnalyzer") << "BeamHaloSummary is invalid.";
    
    if (!cscHaloDataHandle.isValid())
      edm::LogError("HaloFilterPerformanceAnalyzer") << "CSCHaloData is invalid.";

    if ( !recoVertices.isValid() )
      edm::LogError("MissingProduct") << "Vertices not found. Branch will not be filled" << std::endl;
    
    return;
  }
  
  // TODO: make this contingent on a bool that is passed in from the config file
  analyzeHaloWithMET(iEvent, iSetup, eventPrinted, cscLooseId, cscTightId, nVtx);
}

void HaloFilterPerformanceAnalyzer::analyzeHaloWithMET(const edm::Event& iEvent, 
						       const edm::EventSetup& iSetup,
						       bool eventPrinted, 
						       const bool cscLooseId, 
						       const bool cscTightId, 
						       const int nVtx) {

  edm::Handle< reco::CaloMETCollection > TheCaloMET;
  iEvent.getByLabel("met", TheCaloMET); 
  const reco::CaloMETCollection *calometcol = TheCaloMET.product();
  const reco::CaloMET *calomet = &(calometcol->front());

  edm::Handle< reco::PFMETCollection > ThePFMET;
  iEvent.getByLabel("pfMet", ThePFMET); 
  const reco::PFMETCollection *pfmetcol = ThePFMET.product();
  const reco::PFMET *pfmet = &(pfmetcol->front());

  //Get CSC Geometry
  edm::ESHandle<CSCGeometry> cscGeometry;
  iSetup.get<MuonGeometryRecord>().get(cscGeometry);
  
  //Get CSC Segments
  edm::Handle<CSCSegmentCollection> cscSegments;
  iEvent.getByLabel("cscSegments", cscSegments);
  if (cscSegments.isValid() && cscSegments->size() < 1) // No csc activity, continue
    return;

  // TODO: Move this logic out into an EDFilter
  //////////  Unbiased Halo Tag and Probe
  if( pfmet->pt() > minimumMET_ ) {
    //if( calomet->pt() > minimumMET_ ) {

    hHaloWithMET_.SetBinContent(1, hHaloWithMET_.GetBinContent(1)+1);

    //float met_phi = calomet->phi();
    float met_phi = pfmet->phi();
    bool tagPlus=false;
    bool tagMinus = false;
    bool oppositeMuon = false;
    bool dphi_bins[11]; for(int i = 0 ; i < 11 ;i++ ) dphi_bins[i]=false; 
      
    for(CSCSegmentCollection::const_iterator iSegment = cscSegments->begin();
	iSegment != cscSegments->end();
	iSegment++) {
  
      CSCDetId iId  = (CSCDetId)iSegment->cscDetId();
      LocalPoint iLocalPos = iSegment->localPosition();
      GlobalPoint iGlobalPos = cscGeometry->chamber(iId)->toGlobal(iLocalPos);
      float iSegment_phi = iGlobalPos.phi();
      
      iSegment_phi = iSegment_phi > TMath::Pi() ? iSegment_phi - 2.* TMath::Pi() : iSegment_phi;
      float dphi = TMath::ACos( TMath::Cos( iSegment_phi - met_phi ));
      hDphiCSCSegmentMET_.Fill(dphi);
      hCrudeDphiCSCSegmentMET_.Fill(dphi);
      // This is the condition for CSCSegment opposite the MET vector, but before
      // tagging it, make sure this isn't a CSCSegment that is matched to a 
      // collision muon.
      if( dphi >= 3. && dphi <= TMath::Pi() ) {	
 
        bool hasCollisionMatch = false;

	/******************* BEGIN MUON SEARCH *******************/
	// TODO: tags should be moved to cfg file
	edm::Handle<reco::MuonCollection> collisionMuons;
	iEvent.getByLabel("muons", collisionMuons);

	if( collisionMuons.isValid()) {
	  reco::MuonCollection::const_iterator mu;
	  for(mu = collisionMuons->begin() ; mu != collisionMuons->end() ; mu++ ) {
	    if( !mu->isTrackerMuon() ) continue;	  
	    
	    //float dphi = -0.05;
	    if( fabs(mu->eta()) > 1.0  && fabs(mu->eta()) <= 3.0 ) {
	      const std::vector<reco::MuonChamberMatch> chambers = mu->matches();

	      std::vector<reco::MuonChamberMatch>::const_iterator iChamber;
	      for(iChamber = chambers.begin(); iChamber != chambers.end(); iChamber ++ ) {
		if( iChamber->detector() != MuonSubdetId::CSC ) continue;

		std::vector<reco::MuonSegmentMatch>::const_iterator jSegmentMatch;
		for(jSegmentMatch = iChamber->segmentMatches.begin(); 
		    jSegmentMatch != iChamber->segmentMatches.end(); ++jSegmentMatch ) {
		  edm::Ref<CSCSegmentCollection> jSegment = jSegmentMatch->cscSegmentRef;

		  // Check if iSegment == jSegment by looking at global coordinates
		  CSCDetId jId  = (CSCDetId)jSegment->cscDetId();
		  LocalPoint jLocalPos = jSegment->localPosition();
		  GlobalPoint jGlobalPos = cscGeometry->chamber(jId)->toGlobal(jLocalPos);
		  
		  if (iGlobalPos == jGlobalPos) hasCollisionMatch = true;
		}
	      }
	    }
	  }
	}
	/** END MUON SEARCH **/
	oppositeMuon = true;
	if (hasCollisionMatch) {
	  // Do we want to record anything about the segments that are associated 
	  // with collision muons?
	} else {
	  if( iGlobalPos.z() < 0. ) 
	    tagMinus = true;
	  else                     
	    tagPlus = true;
	  dphi_bins[10]=true;
	  hCSCSegmentPhi_.Fill(iSegment_phi); // histogram phi for all tagged CSCSegments
	  hDphiProbeCSCSegmentMET_.Fill(dphi);

	}
      } 

      // CSCSegment NOT opposite the MET vector, so don't tag this.
      else if ( iGlobalPos.z() > 0. ){ 
	int bin = dphi/0.3;
	if( bin < 10 ) 
	  dphi_bins[bin] = true;
	else
	  std::cerr << "ERROR dphi bin is out of range,  dphi = " << dphi << "    bin = " << bin <<std::endl;
      }

      if (printEventInfo_) {
	if (!eventPrinted) {
	  std::cout << "================= " 
		    <<  iEvent.id().run() << ":" <<iEvent.luminosityBlock() << ":" << iEvent.id().event()
		    << " =================" << std::endl;
	  eventPrinted = true;
	}
	std::cout << pfmet->pt() << "\t"
		  << pfmet->phi() << "\t"
		  << calomet->pt() << "\t" 
		  << calomet->phi() << "\t"
		  << iSegment_phi << "\t"
		  << dphi << "\t"
		  << tagMinus << "\t"
		  << tagPlus << "\t"
		  << std::endl;
	
      }
    }
    
    // Record MET and MET phi for all events tagged.
    if (oppositeMuon) {
      hHaloWithMET_.SetBinContent(2, hHaloWithMET_.GetBinContent(2)+1); // MET>50GeV && Opposite CSCSegment
      hnVtxMETAll_.Fill(nVtx);
    }

    if (tagPlus || tagMinus) {
      hHaloWithMET_.SetBinContent(3, hHaloWithMET_.GetBinContent(3)+1); // MET>50GeV && Opposite NC-CSCSegment
      //hMET_.Fill(calomet->pt());
      hMET_.Fill(pfmet->pt());
      hMETPhi_.Fill(met_phi);
      metCSCEvents_ <<  iEvent.id().run() << ":" 
		    <<iEvent.luminosityBlock() << ":" 
		    << iEvent.id().event()
		    << std::endl;
      
      // Probe (ie, if event has CSCSegment opposite MET vector, is it halo?)
      if (cscLooseId) {
	hHaloWithMET_.SetBinContent(4, hHaloWithMET_.GetBinContent(4)+1);
	hnVtxMETHaloLoose_.Fill(nVtx);
	metCSCLooseEvents_ << iEvent.id().run() << ":" 
			   << iEvent.luminosityBlock() << ":" 
			   << iEvent.id().event()
			   << std::endl;
      }
      if (cscTightId) { 
	hHaloWithMET_.SetBinContent(5, hHaloWithMET_.GetBinContent(5)+1);
	hnVtxMETHaloTight_.Fill(nVtx);
	metCSCTightEvents_ << iEvent.id().run() << ":" 
			   << iEvent.luminosityBlock() << ":" 
			   << iEvent.id().event()
			   << std::endl;
      }
      if (!cscLooseId && !cscTightId) { 
	metCSCNeitherEvents_ << iEvent.id().run() << ":" 
			     << iEvent.luminosityBlock() << ":" 
			     << iEvent.id().event()
			     << std::endl;
      }
    }
  }
}

// ------------ method called when starting to processes a run  ------------
void HaloFilterPerformanceAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run  ------------
void HaloFilterPerformanceAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void HaloFilterPerformanceAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void HaloFilterPerformanceAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HaloFilterPerformanceAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HaloFilterPerformanceAnalyzer);
