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

//Muons
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"      //not sure  is it needed?


//Missing Energy
//==============
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/CorrMETData.h"



//Header file for Fastjet Analysis
//====================================

#include "fastjet/ClusterSequence.hh"
#include "fastjet/config.h"
#include "fastjet/SISConePlugin.hh"

#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "FastjetEx/FastJetSimple1/interface/Myheaderfile1.h"       // My header file
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"



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
	edm::EDGetTokenT<std::vector<pat::Muon>			> patmuonToken;
	edm::EDGetTokenT<std::vector<pat::Jet>			> patjetToken;
	edm::EDGetTokenT<std::vector<pat::MET>			> patMetToken;
	edm::EDGetTokenT<std::vector<pat::IsolatedTrack>	> patIsolatedTrackToken;
	
	
	edm::Service<TFileService> fs;
	

// Defining Histogtams-------------------------------------------------------------------------
	
// For muons----------------------------------
	
	

// For Jets------------------------------------
    TH1F *hist_njets; 
	TH1F *hist_jetspt;
	TH1F *hist_jetseta;
	TH1F *hist_jetsphi;
	TH1F *hist_jetsmass;
	TH1F *hist_dausPID;
	
	
	
	
	
	TH1F *hist_n_inc_jets;
	TH1F *hist_n_exc_jets;
// For Mets-----------------------------------
	
	
	
	
/*	TH1F *hist_metsumEt;
	TH1F *hist_metet;
	TH1F *hist_meteta;
	TH1F *hist_metphi;    
*/
	
	TH1F *hist_metpt;
	//----------------------------------------
	
	
	
	

};


FastJetSimple1::FastJetSimple1(const edm::ParameterSet& iConfig)
 :
patmuonToken(consumes<std::vector<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>		("muonTag"))),
patjetToken(consumes<std::vector<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>		("jetTag"))),
patMetToken(consumes<std::vector<pat::MET> >(iConfig.getUntrackedParameter<edm::InputTag>		("metTag"))),
patIsolatedTrackToken(consumes<std::vector<pat::IsolatedTrack> >(iConfig.getUntrackedParameter<edm::InputTag>	("trackTag")))	
{
    
	// For muons
	
	
	
	// For jets----------------------------------------
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

	
	hist_dausPID = fs->make<TH1F>("DausPID", "Daus PdgID",500,-250,250);
	
	
	hist_n_inc_jets = fs->make<TH1F>("NInc_Jets", "Number of Inclusive Jets", 10, 0, 10);
	hist_n_inc_jets->SetTitle("Number of Inclusive Jets in events");
	hist_n_inc_jets->GetXaxis()->SetTitle("Number of Inclusive Jets");
	hist_n_inc_jets->GetYaxis()->SetTitle("Number of events");
	hist_n_inc_jets->SetFillStyle( 3001);
    hist_n_inc_jets->SetFillColor( kRed);
		
		
		
	hist_n_exc_jets = fs->make<TH1F>("NExc_Jets", "Number of Exclusive Jets", 10, 0, 10);
	
	
	
/*	// For Mets----------------------------------------

	hist_metsumEt= fs->make<TH1F>("metsumEt", "metsumEt",600,0,300);
	hist_metet= fs->make<TH1F>("metet", "metet",250,0,250);
	hist_meteta= fs->make<TH1F>("meteta", "meteta",250,0,250);
	hist_metphi= fs->make<TH1F>("metphi", "metphi",20,0,5);                    
*/	
	hist_metpt=fs->make<TH1F>("metpt", "metpt",200,0,200);
	
	
}

FastJetSimple1::~FastJetSimple1()
{
}

void FastJetSimple1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	edm::Handle<std::vector<pat::Muon>> patMuon;
	iEvent.getByToken(patmuonToken, patMuon);
	
	edm::Handle<std::vector<pat::Jet>> patJet;
	iEvent.getByToken(patjetToken, patJet);
	
	edm::Handle<std::vector<pat::MET>> patMet;
	iEvent.getByToken(patMetToken, patMet); 
	
	edm::Handle<std::vector<pat::IsolatedTrack>> patIsolatedTrack;
	iEvent.getByToken(patIsolatedTrackToken, patIsolatedTrack); 
	
	std::vector<fastjet::PseudoJet> input_particles;
	
	
	
//===========================Muons==========================Muons=============================Muons====================
/*
	int m=0;
	std::cout << "Number of RECO muons: " << patMuon->size() << std::endl;
	
	for (std::vector<pat::Muon>::const_iterator itMuon=patMuon->begin(); itMuon!=patMuon->end(); ++itMuon) 
	{
		m=m+1; 
	}
	
	std::cout<<m<<std::endl;
	
*/
	
	

//===========================Fastjet==============================Fastjet==============================Fastjet======================	
	
	for(std::vector<pat::IsolatedTrack>::const_iterator itTrack = patIsolatedTrack->begin(); itTrack != patIsolatedTrack->end(); ++itTrack)
		{
			//int charge = itTrack->pt();
			//std::cout<<charge<<std::endl;
		
			input_particles.push_back(fastjet::PseudoJet(itTrack->px(),itTrack->py(),itTrack->pz(),itTrack->energy()));
   		}
	 
		//input_particles.rap(), input_particles.pt(), input_particles.mt2());
	 
	
		std::cout <<  " Number of particles before applying cuts (ie, before using selector) : " <<input_particles.size() << std::endl;
	
	
		//defining basic set of jet cuts using fastjet::Selector
	
		//fastjet::Selector particle_selector = fastjet::SelectorAbsRapRange(1.0,2.5) || (fastjet::SelectorAbsRapMax(1.0) && fastjet::SelectorPtMin(1.0));
		fastjet::Selector particle_selector = fastjet::SelectorAbsRapRange(1.0,2.5);
		
			
		
		input_particles = particle_selector(input_particles);   // Particles with cut and these are the particles to be used for jet clustering
	
		std::cout <<  " Number of particles after applying cuts : " <<input_particles.size() << std::endl;
		
			
		// Jet definition: specifying how to carry out the clustering
	
		double R = 0.7;
		fastjet::JetDefinition jet_def(fastjet::kt_algorithm,R);
		std::cout<<"Jet definition used here: "<<jet_def.description()<<std::endl;
	
	
		//Jet clustering for particles with above jet definition
		fastjet::ClusterSequence clust_seq(input_particles, jet_def);

	
		
		//Retriving required information from the clustered jet
		
		
		//Inclusive Jets 
		std::vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets();
		std::cout<< "Number  of Inclusive jets = "<<inclusive_jets.size()<<std::endl;
		int incJetSize = inclusive_jets.size();
		hist_n_inc_jets->Fill(incJetSize);	
	
		int pdg_id = 13;   	//  - pdg_id        the PDG id of the particle
   		int vertex_no = 1;	//  - vertex_number the id of the vertex it originates from
		
		printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");   // label the columns
		
		for (unsigned int i = 0; i < inclusive_jets.size(); i++)
		{
			Myheaderfile1* infojet = new Myheaderfile1(pdg_id,vertex_no);	
			inclusive_jets[i].set_user_info(infojet);
			const int & pdgid = inclusive_jets[i].user_info<Myheaderfile1>().pdg_id();
		
			std::cout<<"This is the pgdID"<<pdgid<<std::endl;
			
			printf("%5u %15.8f %15.8f %15.8f\n",i, inclusive_jets[i].rap(), inclusive_jets[i].phi(), inclusive_jets[i].mt2());
		}

		//Exclusive Jets
		if (incJetSize >2)
		{
			std::vector<fastjet::PseudoJet> exclusive_jets = clust_seq.exclusive_jets(1);
			std::cout<<"\n";
			std::cout<< "Number of Exclusive jets = "<<exclusive_jets.size()<<std::endl;
			int ExcJetSize = exclusive_jets.size();
			hist_n_exc_jets->Fill(ExcJetSize);
		}
			
	
//===========================Jet==========================Jet=============================Jet============================	

	int jets = 0;  // this counts the number of jets
    	
	for (std::vector<pat::Jet>::const_iterator itJets=patJet->begin(); itJets!=patJet->end(); ++itJets) 
	{
		jets=jets+1; //counters for jets
	   
		float jetspt = itJets->pt();
		float jetseta  = itJets->eta();
		float jetsphi  = itJets->phi();
		float jetsmass = itJets->mass();
	   
		if(abs(jetseta)<2)
		{
			hist_jetspt -> Fill(jetspt);
			hist_jetseta -> Fill(jetseta);
			hist_jetsphi -> Fill(jetsphi);
			hist_jetsmass -> Fill(jetsmass);
		}
		
		std::vector daus(itJets->daughterPtrVector());
	   	for (unsigned int k =0; k < daus.size(); k++)
		{
           		const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[k]);

			//  int pid =fabs(cand.pdgId());
	   		int pid =(cand.pdgId());
           		//std::cout<<"Particle IDs in the jet"<<pid<<std::endl;
		   	hist_dausPID -> Fill(pid);
             	} 
    	}
	
	hist_njets->Fill(jets); //filling histogram with the number of jets
	    

	
	//===========================MET==========================MET=============================MET============================
	std::cout << "\n";
	const pat::MET &met = patMet->front();
	float metpt = met.pt();
	std::cout << " met.pt " <<  met.pt() << " met.px " <<  met.px() << " met.py " <<  met.py()  << " phi " <<  met.phi() << std::endl;
	
	hist_metpt -> Fill(metpt);
	
	
/*  	
	float metsumEtMax = 0;   
	int n=0;
	for(std::vector<pat::MET>::const_iterator itMets = patMet->begin(); itMets != patMet->end(); ++itMets) 
	{
		n=n+1;	    
		float metsumEt = itMets->sumEt();
		float metet  = itMets->et();
		float meteta  = itMets->eta();
		float metphi = itMets->phi();

		hist_metsumEt -> Fill(metsumEt);
		hist_metet  -> Fill(metet);
		hist_meteta  -> Fill(meteta);
		hist_metphi -> Fill(metphi);
		
		if(metsumEt>metsumEtMax)
		{
			metsumEtMax = metsumEt;
		}
	}
*/
	//std::cout<<n<<std::endl; 
	//std::cout<<metsumEtMax<<std::endl;
	std::cout<<"================================"<<std::endl; 

}

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
