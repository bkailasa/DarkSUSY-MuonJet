// -*- C++ -*-
//
// Package:    FastjetEx/FastJetSimple1
// Class:      FastJetSimple1
//
/**\class FastJetSimple1 FastJetSimple1.cc FastjetEx/FastJetSimple1/plugins/FastJetSimple1.cc
 Description: [one line class summary]
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Balashangar Kailasapathy.
//         Created:  Sat, 07 Nov 2020 01:13:32 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//======Additional header files==================
#include "FWCore/Framework/interface/EventSetup.h" 
#include "FWCore/Framework/interface/ESHandle.h" 
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" //to work with reco::GenParticle
#include "DataFormats/PatCandidates/interface/Jet.h" //to work with pat::Jets
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" //to work with particles inside jets
#include <vector> 
#include <string> 
#include <map>
#include <iostream>
#include <fstream>
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/ServiceRegistry/interface/Service.h" // to use TFileService
#include "CommonTools/UtilAlgos/interface/TFileService.h" // to use TFileService



class FastJetSimple1 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit FastJetSimple1(const edm::ParameterSet&);
      ~FastJetSimple1();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
	edm::EDGetTokenT<std::vector<pat::Jet>> patjetToken;
	edm::Service<TFileService> fs;

// For Jets------------------------------------
    TH1F *hist_njets; 
	TH1F *hist_jetspt;
	TH1F *hist_dausPID;
};


FastJetSimple1::FastJetSimple1(const edm::ParameterSet& iConfig)
 :
patjetToken(consumes<std::vector<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("jetpat")))
{
    // Histogram for number of jets----------------------------------------
	hist_njets = fs->make<TH1F>("NJets", "Number of Jets", 12, -1.5, 10.5);
	hist_njets->SetTitle("Number of Jets in events");
	hist_njets->GetXaxis()->SetTitle("Number of Jets");
	hist_njets->GetYaxis()->SetTitle("Number of events");
	hist_njets->SetFillStyle( 3001);
    hist_njets->SetFillColor( kRed);
	//---------------------------------------------------------------------
	
	hist_jetspt = fs->make<TH1F>("Jetspt", "Jets", 200, 0.0, 200.0);	
	hist_jetspt->SetTitle("Transverse momentum of jets");
	hist_jetspt->GetXaxis()->SetTitle("Number of Jets");
	hist_jetspt->GetYaxis()->SetTitle("Transverse Momentum");
	hist_jetspt->SetFillStyle( 3011);
    hist_jetspt->SetFillColor( kGreen);
	//---------------------------------------------------------------------
	
	hist_dausPID = fs->make<TH1F>("DausPID", "Daus PdgID",250,0,250);
}

FastJetSimple1::~FastJetSimple1()
{
}

void FastJetSimple1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

	edm::Handle<std::vector<pat::Jet>> patjet;
	iEvent.getByToken(patjetToken, patjet);
	
    int jets = 0;  // this counts the number of jets
    	
    for (std::vector<pat::Jet>::const_iterator itJets=patjet->begin(); itJets!=patjet->end(); ++itJets) {
       jets=jets+1; //counters for jets
	   
	   //Transverse Momentum of Jets------------------------------------
	   float jetspt = itJets->pt();
	   hist_jetspt -> Fill(jetspt);
	   
	   // ---------------
	   
	   
	   std::vector daus(itJets->daughterPtrVector());
	   for (unsigned int k =0; k < daus.size(); k++){
           const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[k]);
           int pid =fabs(cand.pdgId());
           //std::cout<<"Particle IDs in the jet"<<pid<<std::endl;
		   hist_dausPID -> Fill(pid);
             } 
	   
    }
    hist_njets->Fill(jets); //filling histogram with the number of jets
}

//End of Jets----------------------------------------------------------------------------------------------------------------------


// ------------ method called once each job just before starting event loop  ------------
void
FastJetSimple1::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
FastJetSimple1::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FastJetSimple1::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FastJetSimple1);
