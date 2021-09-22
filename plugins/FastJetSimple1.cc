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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"      


//Missing Energy
//==============
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/CorrMETData.h"


#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TTree.h"


//Header file for Fastjet Analysis
//====================================

#include "fastjet/ClusterSequence.hh"
#include "fastjet/config.h"
#include "fastjet/SISConePlugin.hh"

//#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "FastjetEx/FastJetSimple1/interface/Myheaderfile1.h"       // My header file
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"

TFile hfile("htree.root","RECREATE","ROOT file with histograms & trees");
TTree tree("Tree","A ROOT tree ");
	
class FastJetSimple1 : public edm::one::EDAnalyzer<edm::one::SharedResources>  
{
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
	edm::EDGetTokenT<std::vector<pat::MET>			> patMetToken;
		
	edm::Service<TFileService> fs;
	
//For Fastjet---------------------------------	
	TH1F *hist_muonpt;
	TH1F *hist_jetrap;
	TH1F *hist_jetphi;
	TH1F *hist_jetpt;
	TH1F *hist_invmass;
	
	TH1F *hist_n_inc_jets;
	TH1F *hist_n_exc_jets;
	
	
// For Mets-----------------------------------
	TH1F *hist_metpt;

//----------------------------------------
};


FastJetSimple1::FastJetSimple1(const edm::ParameterSet& iConfig):
	patmuonToken(consumes<std::vector<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>	("muonTag"))),
	patMetToken(consumes<std::vector<pat::MET> >(iConfig.getUntrackedParameter<edm::InputTag>	("metTag")))
{
    
	//---------------------hists for Fastjet::pseudojet----------------------------------
	hist_muonpt = fs->make<TH1F>("muonpt", "muon pt (without cuts",10,0,10);
	hist_jetrap =fs->make<TH1F>("jetrap", "Jet Rapidity", 10, 0, 10);
	hist_jetphi =fs->make<TH1F>("jetphi", "Jet phi", 10, 0, 10);
	hist_jetpt = fs->make<TH1F>("jetpt", "Jet pt", 10, 0, 10);
	hist_invmass  = fs->make<TH1F>("invmass", "Jet Invarient mass", 10, 0, 10);
	
	
	hist_n_inc_jets = fs->make<TH1F>("NInc_Jets", "Number of Inclusive Jets", 10, 0, 10);
	hist_n_inc_jets->SetTitle("Number of Inclusive Jets in events");
	hist_n_inc_jets->GetXaxis()->SetTitle("Number of Inclusive Jets");
	hist_n_inc_jets->GetYaxis()->SetTitle("Number of events");
	hist_n_inc_jets->SetFillStyle( 3001);
    	hist_n_inc_jets->SetFillColor( kRed);
		
		
		
	hist_n_exc_jets = fs->make<TH1F>("NExc_Jets", "Number of Exclusive Jets", 10, 0, 10);
	
	
	
	// ---------------------------hists for pat::MET----------------------------------------
	hist_metpt=fs->make<TH1F>("metpt", "metpt",200,0,200);
	
}

FastJetSimple1::~FastJetSimple1()
{
}

void FastJetSimple1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	edm::Handle<std::vector<pat::Muon>> patMuon;
	iEvent.getByToken(patmuonToken, patMuon);
	
	edm::Handle<std::vector<pat::MET>> patMet;
	iEvent.getByToken(patMetToken, patMet); 
	
	std::vector<fastjet::PseudoJet> input_particles;
	
	
	

//===========================Fastjet==============================Fastjet==============================Fastjet======================	
	
//		for(std::vector<pat::IsolatedTrack>::const_iterator itTrack = patIsolatedTrack->begin(); itTrack != patIsolatedTrack->end(); ++itTrack)
		for (std::vector<pat::Muon>::const_iterator itMuon=patMuon->begin(); itMuon!=patMuon->end(); ++itMuon) 
		{
			int muonpt = itMuon->pt();
			hist_muonpt->Fill(muonpt);
					
			input_particles.push_back(fastjet::PseudoJet(itMuon->px(),itMuon->py(),itMuon->pz(),itMuon->energy()));
		}
	 
		std::cout <<  " Number of particles before applying cuts (ie, before using selector) : " <<input_particles.size() << std::endl;
	
	
		//defining basic set of jet cuts using fastjet::Selector
		//fastjet::Selector particle_selector = fastjet::SelectorAbsRapRange(1.0,2.5) || (fastjet::SelectorAbsRapMax(1.0) && fastjet::SelectorPtMin(1.0));
		fastjet::Selector particle_selector = fastjet::SelectorAbsRapRange(1.0,2.5);
		selected_input_particles = particle_selector(input_particles);   // Particles with cut and these are the particles to be used for jet clustering
	
		std::cout <<  " Number of particles after applying cuts : " <<selected_input_particles.size() << std::endl;
		
			
		
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
		
		printf("%5s %15s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt" , "invMass");   // label the columns
		
		for (unsigned int i = 0; i < inclusive_jets.size(); i++)
		{
			Myheaderfile1* infojet = new Myheaderfile1(pdg_id,vertex_no);	
			inclusive_jets[i].set_user_info(infojet);
			const int & pdgid = inclusive_jets[i].user_info<Myheaderfile1>().pdg_id();
		
			std::cout<<"This is the pgdID: "<<pdgid<<std::endl;
			
			printf("%5u %15.8f %15.8f %15.8f %15.8f\n",i, inclusive_jets[i].rap(), inclusive_jets[i].phi(), inclusive_jets[i].pt(), inclusive_jets[i].m());
			
			
			float jetrap  = inclusive_jets[i].rap();
			float jetphi = inclusive_jets[i].phi();
			float jetpt = inclusive_jets[i].pt();
			float invmass = inclusive_jets[i].m();
			
			
			hist_jetrap -> Fill(jetrap);
			
			tree.Branch("JetRapidity",&hist_jetrap);
			
			
			
			
			hist_jetphi -> Fill(jetphi);
			hist_jetpt -> Fill(jetpt);
			hist_invmass -> Fill(invmass);
			
			
			
			tree.Fill();
			
			
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
		
		 tree.Print();
 
    hfile.Write();
 
    hfile.Close();
		
		
		std::cout<<"******************************************************************************************"<<std::endl; 
	
	
//===========================END Fastjet================================================================================================

//===========================pat::MET==========================pat::MET=============================pat::MET============================
	std::cout << "\n";
	const pat::MET &met = patMet->front();
	float metpt = met.pt();
	std::cout << " met.pt " <<  met.pt() << " met.px " <<  met.px() << " met.py " <<  met.py()  << " phi " <<  met.phi() << std::endl;
	hist_metpt -> Fill(metpt);
	//==========================END pat::MET============================================================================================================ 

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
