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
// Original Author:  Balashangar Kailasapathy. New
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

//Missing Energy
//==============
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/CorrMETData.h"



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
	edm::EDGetTokenT<std::vector<pat::Jet>		> patjetToken;
	edm::EDGetTokenT<std::vector<pat::MET>		> patMetToken;
	
	edm::Service<TFileService> fs;

// For Jets------------------------------------
    TH1F *hist_njets; 
	TH1F *hist_jetspt;
	TH1F *hist_jetseta;
	TH1F *hist_jetsphi;
	TH1F *hist_jetsmass;
	TH1F *hist_dausPID;
	
// For Mets-----------------------------------
	TH1F *hist_metsumEt;
	TH1F *hist_metet;
	TH1F *hist_meteta;
	TH1F *hist_metphi;
	//----------------------------------------	
};


FastJetSimple1::FastJetSimple1(const edm::ParameterSet& iConfig)
 :
patjetToken(consumes<std::vector<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("jetTag"))),
patMetToken(consumes<std::vector<pat::MET> >(iConfig.getUntrackedParameter<edm::InputTag>("metTag")))
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
	hist_jetspt->GetXaxis()->SetTitle("Transverse Momentum");
	hist_jetspt->GetYaxis()->SetTitle("Number of Jets");
	hist_jetspt->SetFillStyle( 3011);
    	hist_jetspt->SetFillColor( kGreen);
	//---------------------------------------------------------------------
	hist_jetseta = fs->make<TH1F>("Jetseta", "Jetseta",20,-10,10);
	hist_jetseta->SetTitle("Eta of jets");
	hist_jetseta->GetXaxis()->SetTitle("eta");
	hist_jetseta->GetYaxis()->SetTitle("Number of Jets");
	hist_jetseta->SetFillStyle( 3021);
    	hist_jetseta->SetFillColor( kBlue);
	//---------------------------------------------------------------------
	hist_jetsphi = fs->make<TH1F>("Jetsphi", "Jetsphi",10,-5,5);
	hist_jetsphi->SetTitle("Phi of jets");
	hist_jetsphi->GetXaxis()->SetTitle("Phi");
	hist_jetsphi->GetYaxis()->SetTitle("Number of Jets");
	hist_jetsphi->SetFillStyle( 3012);
    	hist_jetsphi->SetFillColor( kPink);
	//---------------------------------------------------------------------
	hist_jetsmass = fs->make<TH1F>("Jetsmass", "Jetsmass",15,0,15);
	hist_jetsmass->SetTitle("Mass of jets");
	hist_jetsmass->GetXaxis()->SetTitle("Mass");
	hist_jetsmass->GetYaxis()->SetTitle("Number of Jets");
	hist_jetsmass->SetFillStyle( 3022);
    	hist_jetsmass->SetFillColor( kOrange);
	//---------------------------------------------------------------------

	
	hist_dausPID = fs->make<TH1F>("DausPID", "Daus PdgID",250,0,250);
	
	hist_metsumEt= fs->make<TH1F>("metsumEt", "metsumEt",250,0,250);
	hist_metet= fs->make<TH1F>("metet", "metet",250,0,250);
	hist_meteta= fs->make<TH1F>("meteta", "meteta",250,0,250);
	hist_metphi= fs->make<TH1F>("metphi", "metphi",250,0,250);
	
	
	
	
}

FastJetSimple1::~FastJetSimple1()
{
}

void FastJetSimple1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

	edm::Handle<std::vector<pat::Jet>> patjet;
	iEvent.getByToken(patjetToken, patjet);
	
	edm::Handle<std::vector<pat::MET>> patmet;
	iEvent.getByToken(patMetToken, patmet); 
	
//Jets
	
    int jets = 0;  // this counts the number of jets
    	
    for (std::vector<pat::Jet>::const_iterator itJets=patjet->begin(); itJets!=patjet->end(); ++itJets) {
       jets=jets+1; //counters for jets
	   
	   float jetspt = itJets->pt();
	   float jetseta  = itJets->eta();
	   float jetsphi  = itJets->phi();
	   float jetsmass = itJets->mass();
	   
	   if(abs(jetseta)<2){
	   hist_jetspt -> Fill(jetspt);
	   hist_jetseta -> Fill(jetseta);
	   hist_jetsphi -> Fill(jetsphi);
	   hist_jetsmass -> Fill(jetsmass);
	   }
	   // ---------------
	   
	   
	   std::vector daus(itJets->daughterPtrVector());
	   for (unsigned int k =0; k < daus.size(); k++){
           const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[k]);
           int pid =fabs(cand.pdgId());
           //std::cout<<"Particle IDs in the jet"<<pid<<std::endl;
		   hist_dausPID -> Fill(pid);
             } 
//Missing Energy
	  
	    
	    
	    
	float metsumEtMax = 0
	for(std::vector<pat::MET>::const_iterator itMets = patmet->begin(); itMets != patmet->end(); ++itMets) {
		    
	   float metsumEt = itMets->sumEt();
	   float metet  = itMets->et();
	   float meteta  = itMets->eta();
	   float metphi = itMets->phi();

	   hist_metsumEt -> Fill(metsumEt);
	   hist_metet  -> Fill(metet);
	   hist_meteta  -> Fill(meteta);
	   hist_metphi -> Fill(metphi);
		
		if(metsumEt>metsumEtMax){
			metsumEtMax = metsumEt
			}
	}
	std::cout<<metsumEtMax<<std::endl;
  
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
