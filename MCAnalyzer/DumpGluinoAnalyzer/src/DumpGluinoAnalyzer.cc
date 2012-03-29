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
// $Id: DumpGluinoAnalyzer.cc,v 1.3 2012/02/29 20:09:23 rodenm Exp $
//
//


// system include files
#include <memory>
#include <fstream>
#include <cmath>
#include <vector>
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

#include "TLorentzVector.h"


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
  
  // Helper functions
  double eta(double x, double y, double z, double time);
  void getMasses(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  std::ofstream mStream_;
  std::string mProducer_;
  std::string hepProducer_;
  unsigned stopped_count_;
  unsigned hb_count_;
  unsigned he_count_;
  unsigned eb_count_;
  unsigned ee_count_;
  unsigned mb_count_;
  unsigned me_count_;
  unsigned tracker_count_;
  unsigned detector_count_;
  unsigned cavern_count_;
  unsigned other_count_;
  unsigned total_count_;
  // ----------member data ---------------------------
}
;

//
// constants, enums and typedefs
//
int sparticleID_ = 1000006;//1000021;
int particleID_ = 6;
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
  stopped_count_ = 0;
  total_count_ = 0;
  hb_count_ = 0;
  he_count_ = 0;
  eb_count_ = 0;
  ee_count_ = 0;
  mb_count_ = 0;
  me_count_ = 0;
  tracker_count_ = 0;
  detector_count_ = 0;
  cavern_count_ = 0;
  other_count_ = 0;
}


DumpGluinoAnalyzer::~DumpGluinoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//
void DumpGluinoAnalyzer::getMasses(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::HepMCProduct> hepMC;
  iEvent.getByLabel(hepProducer_, "", hepMC);
  const HepMC::GenEvent* mc = hepMC->GetEvent();
  if( mc == 0 ) 
    throw edm::Exception( edm::errors::InvalidReference ) << "HepMC has null pointer to GenEvent" << std::endl;
  
  //const size_t size = mc->particles_size();
  //std::cout << "particles: " << size << std::endl;
  mc->print( std::cout );      
  
  unsigned sparticle_count = 0;
  unsigned neutralino_count = 0;
  unsigned particle_count = 0;
  std::vector<double> sparticleM;
  std::vector<double> neutralinoM;
  std::vector<double> particleM;
  std::vector<double> neutralinoE;
  std::vector<double> particleE;
  // Iterate over the HepMC vertices and particles, look for the gluinos that produce r-hadrons,
  // print kinematic info for the earliest gluino that has the same phi as the stopped r-hadron
  for ( HepMC::GenEvent::particle_const_iterator piter  = mc->particles_begin();
	piter != mc->particles_end(); 
	++piter ) {
    HepMC::GenParticle* p = *piter;
    int partId = p->pdg_id();
    
    if (partId == sparticleID_ && sparticle_count < 1) {
      math::XYZTLorentzVector momentum1(p->momentum().px(),
					p->momentum().py(),
					p->momentum().pz(),
					p->momentum().e());
      //double phi = momentum1.phi();
      //double pt = momentum1.pt();
      //double e = momentum1.e();
      //double eta = momentum1.eta();
      double mass = p->momentum().m();  

      sparticleM.push_back(mass);
      sparticle_count++;
    }
    if (partId > 1000021 && partId < 1000038 && neutralino_count < 1) {
      math::XYZTLorentzVector momentum1(p->momentum().px(),
					p->momentum().py(),
					p->momentum().pz(),
					p->momentum().e());
      //double phi = momentum1.phi();
      //double pt = momentum1.pt();
      double e = momentum1.e();
      //double eta = momentum1.eta();
      double mass = p->momentum().m();  

      neutralinoM.push_back(mass);
      neutralinoE.push_back(e);
      neutralino_count++;
    }
    if ( (partId == particleID_) && particle_count < 1) {
      math::XYZTLorentzVector momentum1(p->momentum().px(),
					p->momentum().py(),
					p->momentum().pz(),
					p->momentum().e());
      //double phi = momentum1.phi();
      //double pt = momentum1.pt();
      double e = momentum1.e();
      //double eta = momentum1.eta();
      double mass = p->momentum().m();  

      particleM.push_back(mass);
      particleE.push_back(e);
      particle_count++;
    }
    if ( (partId >= 1000600) ) {
      math::XYZTLorentzVector momentum1(p->momentum().px(),
					p->momentum().py(),
					p->momentum().pz(),
					p->momentum().e());
      //double phi = momentum1.phi();
      //double pt = momentum1.pt();
      double e = momentum1.e();
      //double eta = momentum1.eta();
      double mass = p->momentum().m();  
      
      //mStream_ << partId << "\t" << mass << "\t" << e << std::endl;
    }
  }
  // If the sparticle never decays, the following print statement will segfault
  if (sparticle_count > particle_count || sparticle_count > neutralino_count) {
    for (unsigned i = 0; i < sparticle_count; i++) {
      mStream_ << sparticleM[i] << std::endl;
    }
  } else {
    for (unsigned i = 0; i < sparticle_count; i++) {
      mStream_ << sparticleM[i] << "\t" << neutralinoM[i] << "\t" << particleM[i] 
	       << "\t" << neutralinoE[i] << "\t" << particleE[i] << std::endl;
    }
  }
  
}

// ------------ method called for each event  ------------
void DumpGluinoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //using namespace edm;
  total_count_++;

  // Use this function to print the sparticle, neutralino, and daughter masses
  getMasses(iEvent, iSetup);
  return;

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
  edm::Handle<std::vector<float> > ts;
  iEvent.getByLabel (mProducer_, "StoppedParticlesTime", ts);
  
  if (names->size() != xs->size() || xs->size() != ys->size() || ys->size() != zs->size()) {
    edm::LogError ("DumpGluinoAnalyzer") << "mismatch array sizes name/x/y/z:"
					 << names->size() << '/' << xs->size() << '/' << ys->size() << '/' << zs->size()
					 << std::endl;
    return;
  } else {
    if (names->size() > 0) {
      stopped_count_++;
      // mStream_ << "### " << names->size() << " stopped particle(s)." << std::endl;
      for (size_t i = 0; i < names->size(); ++i) {
	float phi = ((*ys)[i]==0 && (*xs)[i]==0) ? 0 : atan2((*ys)[i],(*xs)[i]);
	//mStream_ << std::endl << (*names)[i] << ' ' << (*xs)[i] << ' ' << (*ys)[i] << ' ' << (*zs)[i] << ' ' << phi ;
	if (!filled1) {
	  stop1 << std::endl << (*names)[i] << ' ' << (*xs)[i]/10.0 << ' ' << (*ys)[i]/10.0 
		<< ' ' << (*zs)[i]/10.0 << ' ' << phi ;
	  filled1 = true;
	} 
	//else {
	//stop2 << std::endl << (*names)[i] << ' ' << (*xs)[i]/10.0 << ' ' << (*ys)[i]/10.0 
	//<< ' ' << (*zs)[i]/10.0 << ' ' << phi ;
	//  filled2 = true;
	//}
	stoppedPhis_g.push_back(phi);
	stoppedPhis_r.push_back(phi);
      }
      
      double r = sqrt((*xs)[0]*(*xs)[0] + (*ys)[0]*(*ys)[0])/10.0;
      double z = (*zs)[0]/10.0;
      double particle_eta = eta((*xs)[0], (*ys)[0], (*zs)[0], (*ts)[0]);
      
      // Identify which detector region the particles stopped in. For ME and MB, this
      // definition includes the entire muon system, not just the yokes.
      if (r < 131.0 && fabs(particle_eta) <= 2.5 && fabs(z) < 300.0) { // TRACKER
        tracker_count_++;
      } else if (r>=131.0 && r<184.0 && fabs(z)<376.0 && fabs(particle_eta)<1.479) { // EB
        eb_count_++;
      } else if (fabs(z)<376.0 && fabs(z) >= 300.0 && fabs(particle_eta)>=1.479
                 && fabs(particle_eta)<3.0){ // EE
        ee_count_++;
      } else if (r>=184.0 && r<295.0 && fabs(particle_eta)<1.3 && fabs(z)<500.0) { // HB
        hb_count_++;
      } else if (fabs(z)<560.0 && fabs(z)>=376.0 && fabs(particle_eta)>=1.3
                 && fabs(particle_eta)<3.0) { // HE
	he_count_++;
      } else if (r>=295.0 && r<728.5 && fabs(z)<675.0) { // MB
        mb_count_++;
      } else if (r>=267.3 && r<728.5 && fabs(z)>=675.0 && fabs(z)<1080.0) { // ME-top
	me_count_++;
      } else if (r<267.3 && fabs(particle_eta)<3.0 && fabs(z)>=560.0
                 && fabs(z)<1080.0) { // ME-bottom
        me_count_++;
      } else if (r<728.5 && fabs(z)<1080.0) { // other regions?
        other_count_++;
      }

      if (r >= 728.5 || fabs(z) > 1080)
	cavern_count_++;
      else
	detector_count_++;

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
      /**
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
      */
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
  //mStream_ << "r-hadron x  y  z  phi  pdgid  gluino_phi  gluino_eta  gluino_pT  gluino_E\n"
  //	   << "\tpdgid  gluino_phi  gluino_eta  gluino_pT  gluino_E  pdgid  rhadron_phi  rhadron_eta  rhadron_pT  rhadron_E";
  
  // for dumpMasses
  mStream_ << "m_stop\tm_neutralino\tm_top\tE_neutralino\tE_top\n";
  //mStream_ << "m_gluino\tm_neutralino\tm_gluon\tE_neutralino\tE_gluon\n";
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DumpGluinoAnalyzer::endJob() {
  mStream_ << std::endl << std::endl;
  mStream_ << "-------------------------------" << std::endl;
  mStream_ << "Total count = " << total_count_ << std::endl;
  mStream_ << "Stopped count = " << stopped_count_ << std::endl;

  mStream_ << "Tracker count = " << tracker_count_ << std::endl;
  mStream_ << "EB count = " << eb_count_ << std::endl;
  mStream_ << "EE count = " << ee_count_ << std::endl;
  mStream_ << "HB count = " << hb_count_ << std::endl;
  mStream_ << "HE count = " << he_count_ << std::endl;
  mStream_ << "MB count = " << mb_count_ << std::endl;
  mStream_ << "ME count = " << me_count_ << std::endl << std::endl;
  
  mStream_ << "detector_count = " << detector_count_ << std::endl;
  mStream_ << "cavern count = " << cavern_count_ << std::endl;
  mStream_ << "other count = " << other_count_ << std::endl;
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

/**
 * eta()
 * 
 * Calculates eta (pseudorapidity) given the cartesian coordinates
 */
// TODO: there's some weirdness where if pt = 0, PsuedoRapitidy() returns
//       a bad value. I don't know how it's calculating pt considering
//       the input is just the position...sort this out.
double DumpGluinoAnalyzer::eta(double x, double y, double z, double time) {
  TLorentzVector v = TLorentzVector(x, y, z, time);
  return v.PseudoRapidity();
}



//define this as a plug-in
DEFINE_FWK_MODULE(DumpGluinoAnalyzer);
