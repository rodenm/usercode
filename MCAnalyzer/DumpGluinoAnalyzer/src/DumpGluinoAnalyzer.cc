// -*- C++ -*-
//
// Package:    DumpGluinoAnalyzer
// Class:      DumpGluinoAnalyzer
// 
/**\class DumpGluinoAnalyzer DumpGluinoAnalyzer.cc MCAnalyzer/DumpGluinoAnalyzer/src/DumpGluinoAnalyzer.cc

 Description: Dump info from MC decay chains into a text file. Similar to SimG4Core/CustomPhysics/plugins/RHStopDump.cc

 Implementation:
     Based on https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookWriteFrameworkModule
*/
//
// Original Author:  Marissa Rodenburg
//         Created:  Mon Aug 15 10:18:57 CDT 2011
// $Id$
//
//


// system include files
#include <memory>
#include <fstream>
#include <cmath>
#include <TMath.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "GeneratorInterface/ThePEGInterface/interface/HepMCConverter.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/Math/interface/LorentzVector.h"


//
// class declaration
//
class DumpGluinoAnalyzer : public edm::EDAnalyzer {
public:
  explicit DumpGluinoAnalyzer(const edm::ParameterSet&);
  ~DumpGluinoAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  std::ofstream mStream_;
  std::string mProducer_;
  std::string hepProducer_;
      // ----------member data ---------------------------
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
DumpGluinoAnalyzer::DumpGluinoAnalyzer(const edm::ParameterSet& iConfig)
  : mStream_ (iConfig.getParameter<std::string>("stoppedFile").c_str()),
    mProducer_ (iConfig.getUntrackedParameter<std::string>("producer", "g4SimHits")),
    hepProducer_ (iConfig.getUntrackedParameter<std::string>("producer", "generator"))
{
  //now do what ever initialization is needed 
}


DumpGluinoAnalyzer::~DumpGluinoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void DumpGluinoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //using namespace edm;

  // Collect value of phi for the stopped points. Save for both gluino and r-hadron
  // comparison to monitor how energy changes over the full decay/fragmation chain.
  std::vector<float> stoppedPhis_g;
  std::vector<float> stoppedPhis_r;

  // Save details of gluino and r-hadrons in strings to be written to mStream_ at the end
  std::ostringstream rStream;
  std::ostringstream stop1, stop2;
  bool filled1 = false;
  bool filled2 = false;
  bool filled3 = false;
  bool filled4 = false;
  bool filled5 = false;
  bool filled6 = false;

  edm::Handle<std::vector<std::string> > names;
  iEvent.getByLabel (mProducer_, "StoppedParticlesName", names);
  edm::Handle<std::vector<float> > xs;
  iEvent.getByLabel (mProducer_, "StoppedParticlesX", xs);
  edm::Handle<std::vector<float> > ys;
  iEvent.getByLabel (mProducer_, "StoppedParticlesY", ys);
  edm::Handle<std::vector<float> > zs;
  iEvent.getByLabel (mProducer_, "StoppedParticlesZ", zs);
  
  if (names->size() != xs->size() || xs->size() != ys->size() || ys->size() != zs->size()) {
    edm::LogError ("DumpGluinoAnalyzer") << "mismatch array sizes name/x/y/z:"
					 << names->size() << '/' << xs->size() << '/' << ys->size() << '/' << zs->size()
					 << std::endl;
    return;
  } else {
    if (names->size() > 0) {
      // mStream_ << "### " << names->size() << " stopped particle(s)." << std::endl;
      for (size_t i = 0; i < names->size(); ++i) {
	float phi = ((*ys)[i]==0 && (*xs)[i]==0) ? 0 : atan2((*ys)[i],(*xs)[i]);
	//mStream_ << std::endl << (*names)[i] << ' ' << (*xs)[i] << ' ' << (*ys)[i] << ' ' << (*zs)[i] << ' ' << phi ;
	if (!filled1) {
	  stop1 << std::endl << (*names)[i] << ' ' << (*xs)[i] << ' ' << (*ys)[i] << ' ' << (*zs)[i] << ' ' << phi ;
	  filled1 = true;
	} else {
	  stop2 << std::endl << (*names)[i] << ' ' << (*xs)[i] << ' ' << (*ys)[i] << ' ' << (*zs)[i] << ' ' << phi ;
	  filled2 = true;
	}
	stoppedPhis_g.push_back(phi);
	stoppedPhis_r.push_back(phi);
      }
      
      edm::Handle<edm::HepMCProduct> hepMC;
      iEvent.getByLabel(hepProducer_, "", hepMC);
      const HepMC::GenEvent* mc = hepMC->GetEvent();
      if( mc == 0 ) 
	throw edm::Exception( edm::errors::InvalidReference ) << "HepMC has null pointer to GenEvent" << std::endl;
      //const size_t size = mc->particles_size();
      //std::cout << "particles: " << size << std::endl;
      //mc->print( std::cout );
      
      // Iterate over the HepMC vertices and particles, look for the gluinos that produce r-hadrons,
      // print kinematic info for the earliest gluino that has the same phi as the stopped r-hadron
      for ( HepMC::GenEvent::particle_const_iterator piter  = mc->particles_begin();
	    piter != mc->particles_end(); 
	    ++piter ) {
	HepMC::GenParticle* p = *piter;
	int partId = p->pdg_id();

	// I'm making the (possibly crap) assumption here that R-hadrons all have numbers between 
	// The SUSY particles and 2000000. 
	if (partId == 1000021) {
	  //int barcode = p->barcode();
	  math::XYZTLorentzVector momentum1(p->momentum().px(),
					    p->momentum().py(),
					    p->momentum().pz(),
					    p->momentum().e());
	  double phi = momentum1.phi();
	  double pt = momentum1.pt();
	  double e = momentum1.e();
	  double eta = momentum1.eta();

	  unsigned iPhi = 0;
	  std::vector<float>::iterator iter;
	  for (iter = stoppedPhis_g.begin(); iter < stoppedPhis_g.end(); iter++, iPhi++) {
	    if (TMath::ACos(TMath::Cos(phi - stoppedPhis_g[iPhi])) <= M_PI/1.5) {	      
	      if (iPhi == 0 && !filled3) {
		stop1 << " " << partId << " " << phi << " " << eta << " " << pt << " " << e;
		filled3 = true;
	      } else if (iPhi == 1 && !filled4) {
		stop2 << " " << partId << " " << phi << " " << eta << " " << pt << " " << e;
		filled4 = true;
	      } else if (stoppedPhis_r.size() > 2) {
		edm::LogError ("DumpGluinoAnalyzer") << "Whaaaat? Too many entries in stoppedPhis_g: " 
						     << iPhi << " " << stoppedPhis_g.size() << std::endl;
	      }	      
	      //  barcode  partId  phi  eta  pT  E(GeV)
	      //mStream_ << "\n\t" << barcode << "\t" << abs(partId) << "\t" << phi << "\t" << eta << "\t" << pt << "\t" << e;
	      
	    }
	  }
	}
	// I'm making the (possibly crap) assumption here that R-hadrons all have numbers between 
	// The SUSY particles and 2000000. 
	else if (fabs(partId) > 1000100 && fabs(partId) < 2000000) {
	  int barcode = p->barcode();
	  
	  math::XYZTLorentzVector momentum1(p->momentum().px(),
					    p->momentum().py(),
					    p->momentum().pz(),
					    p->momentum().e());
	  double eta = momentum1.eta();
	  double phi = momentum1.phi();
	  double pt = momentum1.pt();
	  double e = momentum1.e();

	  //mStream_ << "\n\tbarcode = " << barcode << "\tpdgid = " << partId << "\tphi = " << phi 
	  //	   << "\teta = " << eta << "\tpt = " << pt << "\te(GeV) = " << e;
	  
	  rStream << "\n\t" << barcode << "\t" << abs(partId) << "  \t" << phi << "  \t" << eta << "  \t" << pt << "  \t" << e;

	  unsigned iPhi = 0;
	  std::vector<float>::iterator iter;
	  for (iter = stoppedPhis_r.begin(); iter < stoppedPhis_r.end(); iter++, iPhi++) {
            if (TMath::ACos(TMath::Cos(phi - stoppedPhis_r[iPhi])) <= M_PI/1.5) {
	      if (iPhi==0 && !filled5) {
		stop1 << " " << partId << " " << phi << " " << eta << " " << pt << " " << e;
		filled5 = true;
	      } else if (iPhi==1 && !filled6) {
		stop2 << " " << partId << " " << phi << " " << eta << " " << pt << " " << e;
		filled6 = true;
	      } else if (stoppedPhis_r.size() > 2) {
		edm::LogError ("DumpGluinoAnalyzer") << "Whaaaat? Too many entries in stoppedPhis_r: " 
						     << iPhi << " " << stoppedPhis_r.size() << std::endl;
	      }
	    }
	  }
	  
	  // UNCOMMENT the next two lines for production runs
	  //mStream_ << " " << partId << " " << eta << " " << phi << " " << pt << " " << e; //<< std::endl;
	}
      }
    } 
    if (filled1) mStream_ << stop1.str();
    if (filled2) mStream_ << stop2.str();
    mStream_ << rStream.str();
  }
  
  
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}





// ------------ method called once each job just before starting event loop  ------------
void 
DumpGluinoAnalyzer::beginJob(){
  mStream_ << "r-hadron x  y  z  phi  pdgid  gluino_phi  gluino_eta  gluino_pT  gluino_E\n"
	   << "\tpdgid  gluino_phi  gluino_eta  gluino_pT  gluino_E  pdgid  rhadron_phi  rhadron_eta  rhadron_pT  rhadron_E";
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DumpGluinoAnalyzer::endJob() {
  mStream_ << std::endl;
}

// ------------ method called when starting to processes a run  ------------
void 
DumpGluinoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&){
}

// ------------ method called when ending the processing of a run  ------------
void 
DumpGluinoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&){
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DumpGluinoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DumpGluinoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DumpGluinoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DumpGluinoAnalyzer);
