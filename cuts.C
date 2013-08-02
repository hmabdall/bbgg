#include <iostream>
#include <fstream>
#include "Reader.h"
#include "TNtuple.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMathBase.h"
#include <cmath>
#include "TRefArray.h"
#include "TRef.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "chainMaker.hh"
#include "TMathBase.h"
#include "TLatex.h"
#include "TString.h"
#include "TChain.h"
#include "/home/habdalla/Delphes-3.0.9/classes/DelphesClasses.h"
#include "/home/habdalla/Delphes-3.0.9/modules/Delphes.h"
#include "/home/habdalla/Delphes-3.0.9/external/ExRootAnalysis/ExRootTreeReader.h"
//#include <boost/assign/list_of.hpp>
//#include <boost/tuple/tuple.hpp>
//#include <boost/tuple/tuple_comparison.hpp>

#include <string>
#include <string.h>
#include <iostream>
#include <map>
#include <utility>
using namespace std;

//Default values: cuts(TString c, 30.0, 70.0, 90.0, 0.4, 110.0, 150.0, 30.0, 30.0, 120.0, 130.0, 0.4, 0.4, 100.0, 2.0, 350.0)

/*This function acts upon root samples (specifically, in TChain format). The parameter 'c' which it takes is indiciative of the sample to be used. c = h -> HH signal; c = t -> ttH signal; c = z -> ZH signal. The function writes files containing the attributes of the particles of interest, photons and b-quarks, per event (one event per line). 

bjetPT/F : bjetEta/F : delR_BB/F : w_mass/F : M_bb/F : photonPT/F : photonEta/F : delR_GmGm/F : M_gmgm/F : delR_BBGmGm/F : bbggPT/F : bbggEta : 

*/

void writeEventData(TString c){

  // create a chain of root files
	TChain *chain = new TChain("Delphes");

	// create string for the filename that will be written to
	TString filename;

	// if char c = s this means running cuts on signal
	if (c == "h"){
	  cout << "Now writing data for HH signal.... " << endl;
	  chainMakerHH(chain);
	  filename = "sigHH.txt";
	}
	
	// if char c = b this means runnign cuts on the background m
	else if (c == "t" ){
	  cout << "Now writing data for ttH background.... "<< endl;
	  chainMakerttH(chain);
	  filename = "backttH.txt";
	}
		
	else if (c == "z"){
	  cout << "Now writing data for ZH background..... "<< endl;
	  chainMakerZH(chain);
	  filename = "backZH.txt";
	}
			
	//create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	Long64_t numberOfEntries = treeReader->GetEntries();
	
	// get pointers to branches that will be used in this macro
	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");

	// create file and open it
	ofstream file;
	file.open(filename);

	// begin event loop
	for (Int_t jentry = 0; jentry < numberOfEntries; jentry++){
	
	if (jentry % 5000 == 0){cout << "event no: " << jentry << endl;}
	
	// read entry
	treeReader->ReadEntry(jentry);

	// DECLARE VARIABLES (particle attributes)
	Double_t bjet1PT = 0, bjet2PT, bjet1Eta = 0, bjet2Eta, delR_bb = 0, w_mass = 0, M_bb = 0;
	Double_t photon1PT = 0, photon2PT = 0, photon1Eta = 0, photon2Eta = 0, delR_gg = 0, M_gg = 0;
	Double_t delR_bbgg = 0, bbggPT = 0, bbggEta = 0;

	// OBTAIN VARIABLE VALUES
		
	// get bjet1P,bjet2PT,bjet1Eta,bjet2Eta
	// match a pair of reconstructed jets, store their attributes

	//Instatiate vector of non-btagged jets
	vector<Jet*> otrosjets;
		
	//declare and initilaze the event sizes for every branch, 
	// i.e. set of particles for this event
	Int_t photon_size = branchPhoton->GetEntriesFast();
	Int_t muon_size = branchMuon->GetEntriesFast();
	Int_t electron_size = branchElectron->GetEntriesFast();
	Int_t jet_size = branchJet->GetEntriesFast();
	  
	//Initialize numBJets
	Int_t numBJets = 0;
	vector <Jet*>demJets;
	
	//Loop over jets to match the reconstructed and generated jets
		
		for (Int_t jet_entry=0; jet_entry<jet_size; jet_entry++) {
	        Jet *jet = (Jet*) branchJet->At(jet_entry);
	        // Cut jets with eta > 2.4
	        //if (fabs(jet->P4().Eta()) > 2.4) {continue;} // DOUBLE CHECK WITH MAX IF SHOULD CUT FOR ACCEPTANCE
	            //Cut jets with pt < 30
	            //if (jet->P4().Pt() < jet_pt) {continue;}
	            //Find BTag using Loose tagger
		    Bool_t Btag0k = (jet->BTag & (1 << 1));
	            //Cut jets with Mass within 20 eta < 2.4
	            //if (Btag0k == 0) {
	                // INSIDE HERE IS THE CODE FOR THE W-BOSON CUTS
	                //otrosjets.push_back(jet);
	          		//}
	            //Look at Btagged jets
	        if (Btag0k == 1) {
	                //Add one to numBJets
	                numBJets += 1;
	                //Keep track of where the BJets are
	                demJets.push_back(jet);
	        } //End of Btag check
            
		} //End of jet loop

		//blaha
		// only continue if there are two b-jets in this event
		if (numBJets < 2) {continue;}
		  
		//initialize bjet_index1 and 2
		Int_t BJet_index_1 = 0; //Should this be 1 or 0????
		Int_t BJet_index_2 = 1;

		//keep only the 2 most energetic bjets
		if (numBJets > 2) {
	        Double_t ptmax1 = 0.0;
	        Double_t ptmax2 = 0.0;
	        for (Int_t index = 0; index < numBJets; index++) {
		  		Jet *bjet = demJets[index];
		   		Double_t bjetpt = bjet->P4().Pt();
		   		if (bjetpt > ptmax1) {
					ptmax1 = bjetpt;
					BJet_index_1 = index;
		   		}
	        }

	        for (Int_t index = 0; index < numBJets; index++) {
		   		Jet *Bjet = demJets[index];
		   		Double_t Bjetpt = Bjet->P4().Pt(); 
		   		if (Bjetpt > ptmax2) {
					if (index == BJet_index_1) {continue;}
					ptmax2 = Bjetpt;
					BJet_index_2 = index;
		      	}
	        }
		}

	    //instantiate the two bjets
        Jet *bjet1 = demJets[BJet_index_1];
        Jet *bjet2 = demJets[BJet_index_2];

        // get bjet1PT, bjet2PT, bjet1Eta, bjet2Eta
		bjet1PT = bjet1->P4().Pt(), bjet2PT = bjet2->P4().Pt(), bjet1Eta = bjet1->P4().Eta(), bjet2Eta = bjet2->P4().Eta();

		// get delR_bb, M_bb
		delR_bb = sqrt(pow(bjet1PT-bjet2PT,2)+pow(bjet1Eta-bjet2Eta,2));
		
		//cacluate parent mass
        Double_t nrg1, nrg2, px1, px2, py1, py2, pz1, pz2, nrg, px, py, pz;
        nrg1 = bjet1->P4().Energy();
        nrg2 = bjet2->P4().Energy();
        px1 = bjet1->P4().Px();
        px2 = bjet2->P4().Px();
        py1 = bjet1->P4().Py();
        py2 = bjet2->P4().Py();
        pz1 = bjet1->P4().Pz();
        pz2 = bjet2->P4().Pz();
        nrg = nrg1 + nrg2;
        px = px1 + px2;
        py = py1 + py2;
        pz = pz1 + pz2;
        M_bb = sqrt(nrg*nrg - px*px - py*py - pz*pz);

        /*NOW BEGINS THE PHOTON STUFF*/

        // get photon1PT, photon2PT, photon1Eta, photon2Eta
		Int_t index_of_first = 1;
        Int_t index_of_second = 0; //should this be 1 or 0????
		       
        // create vector
        vector <Photon*>luxons;
        
        Int_t numPhotons = 0;
        // loop through photons
        for (int r = 0; r < photon_size; r++){
            Photon *photon = (Photon*)branchPhoton->At(r);
            //if (photon->PT < photon_pt || abs(photon->Eta > 2.4)){continue;}
            luxons.push_back(photon);
            numPhotons++;
        }
        
        if (numPhotons < 2){continue;}
        // cut for phtons < 2
                
        if (numPhotons == 2){
            index_of_first = 0;
            index_of_second = 1;
        }
        
        int photon_PTs[numPhotons];
        //keep only the 2 most energetic photons
        if (numPhotons > 2) {
            
            Double_t ptmax1 = 0.0;
            for (Int_t index = 0; index < numPhotons; index++) {
                Photon *phot = luxons[index];
                Double_t phot_pt = phot->P4().Pt();
                photon_PTs[index] = phot_pt;
                if (phot_pt > ptmax1) {
                    ptmax1 = phot_pt;
                    index_of_first = index;
                }
            }
			
            Double_t ptmax2 = 0.0;
            for (Int_t index2 = 0; index2 < numPhotons; index2++) {
                Photon *Phot = luxons[index2];
                Double_t Phot_pt = Phot->P4().Pt();
                if ((Phot_pt > ptmax2) && (index2 != index_of_first)) {
                    ptmax2 = Phot_pt;
                    index_of_second = index2;
                }
            }
        }
        
        // create the two photons to be used next
        Photon *photon1 = luxons[index_of_first];
        Photon *photon2 = luxons[index_of_second];
        
        // create vector for higgs_on_shell_from_BB
        TLorentzVector Higgs_gmgm_vec = photon1->P4() + photon2->P4();
		photon1PT = photon1->P4().Pt(), photon2PT = photon2->P4().Pt();
		photon1Eta = photon1->P4().Eta(), photon2Eta = photon2->P4().Eta();

       	// get delR_gg
		delR_gg = sqrt(pow(photon1PT-photon2PT,2)+pow(photon1Eta-photon2Eta,2));
	
		// get M_gg
		// find the mass of these photons's mothers
        double nrg_sum = 0, px_sum = 0, py_sum = 0, pz_sum = 0;
        
        // create a lorentz vector for photon 1
        Double_t phot_nrg1,phot_px1,phot_py1,phot_pz1;
        phot_nrg1=photon1->E;
        phot_px1=photon1->P4().Px();
        phot_py1=photon1->P4().Py();
        phot_pz1=photon1->P4().Pz();
        
        // create a lorentz vector fotr photon 2
        Double_t phot_nrg2,phot_px2,phot_py2,phot_pz2;
        phot_nrg2=photon2->E;
        phot_px2=photon2->P4().Px();
        phot_py2=photon2->P4().Py();
        phot_pz2=photon2->P4().Pz();
        nrg_sum = phot_nrg1+phot_nrg2, px_sum = phot_px1+phot_px2, py_sum = phot_py1+phot_py2, pz_sum = phot_pz1+phot_pz2;
        
        M_gg = sqrt((nrg_sum)*(nrg_sum)-(px_sum)*(px_sum)-(py_sum)*(py_sum)-(pz_sum)*(pz_sum));
        // get delR_bbgg, bbggPT, bbggEta (will do later, these are the advanced cuts)

	  	// write variable names (particle attirbutes)
		file << bjet1PT << " " << bjet2PT << " " << bjet1Eta << " " << bjet2Eta << " " << delR_bb << " ";
		file << M_bb << " " << photon1PT << " " << photon2PT << " " << photon1Eta<< " " <<  photon2Eta << " " << delR_gg << " ";
		file << M_gg << "\n";
		// << delR_bbgg << " " << bbggPT << " " << bbggEta << "\n";

	} // end event loop

	file.close(); // close the file once all events are read in


}

double SampleYield(TString c, map<string,int> MapOfCuts){
  Double_t jet_pt = MapOfCuts["jet_pt"];
  Double_t wmass_min = 80 - MapOfCuts["wmass_width"];
  Double_t wmass_max = 80 + MapOfCuts["wmass_width"];
  Double_t bbdeltaR = MapOfCuts["bbdeltaR"];
  Double_t bbmass_min = 130 - MapOfCuts["bbmass_width"];
  Double_t bbmass_max = 130 + MapOfCuts["bbmass_width"];
  Double_t lepton_pt = MapOfCuts["lepton_pt"];
  Double_t photon_pt = MapOfCuts["photon_pt"];
  Double_t ggmass_min = 125 + MapOfCuts["ggmass_width"];
  Double_t ggmass_max = 125 - MapOfCuts["ggmass_width"];
  Double_t ggdeltaR = MapOfCuts["ggdeltaR"];
  Double_t bgdeltaR = MapOfCuts["bbdeltaR"];
  Double_t bbgg_pt = MapOfCuts["bbgg_pet"];
  Double_t bbgg_eta = MapOfCuts["bbgg_eta"];
  Double_t bbgg_mass = MapOfCuts["bbgg_mass"];

	// create a chain of root files
	TChain *chain = new TChain("Delphes");

	// if char c = s this means running cuts on signal
	if (c == "h"){
	  cout << "Now running cuts on the HH sample.... " << endl;
	chainMakerHH(chain);
		
		}
	
	// if char c = b this means runnign cuts on the background m
	else if (c == "t" ){
	  cout << "Now running cuts on the ttH sample.... "<< endl;
		chainMakerttH(chain);
		
		}
		
	else if (c == "z"){
	  cout << "Now running cuts on the ZH sample..... "<< endl;
		chainMakerZH(chain);
		}
			
	//create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	Long64_t numberOfEntries = treeReader->GetEntries();
	
	// get pointers to branches that will be used in this macro
	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	
	// create two histograms to eventually plot the mass of bb and gmgm in the events that make it past all the cuts	
	TH1F *parent_mass_hist_bb = new TH1F("parent_mass_hist_bb", "Apparent Parent Mass of bb Pair", 200, 0, 400);
	TH1F *parent_mass_hist_gmgm = new TH1F("parent_mass_hist_gmgm", "Apparent Parent Mass of GmGm Pair", 200, 0, 200);

	parent_mass_hist_bb->GetXaxis()->SetTitle("M_{bb}");
	parent_mass_hist_bb->GetYaxis()->SetTitle("# events");	
	
	parent_mass_hist_gmgm->GetXaxis()->SetTitle("M_{GmGm}");
	parent_mass_hist_gmgm->GetYaxis()->SetTitle("# events");	
	
	// declare and initialize variables that will be printed out in the macro summary
	
	// general/other cuts
	Int_t less_than_two_Bs_cuts = 0;
	Int_t less_than_two_Gms_cuts = 0;
	Int_t bgm_isol_cuts = 0;

	// cuts against ttH background/advanced cuts
	Int_t lepton_cuts = 0;
	Int_t w_cuts = 0;
	
	// cuts against ZH background
	Int_t diHiggs_mass_cuts = 0; // also an "advanced cut"
		
	// cuts for HH signal background
	Int_t gmgm_mass_cuts = 0;
	Int_t gmgm_isol_cuts = 0;
	Int_t bb_isol_cuts = 0;
	Int_t bb_mass_cuts = 0;

	// advanced cuts
	Int_t Higgs_PT_cuts = 0;
	Int_t Higgs_Eta_cuts = 0;
	
	// begin the event loop
	for (Int_t jentry = 0; jentry < numberOfEntries; jentry++){
	if (jentry % 5000 == 0){cout << "processing event no: " << jentry << endl;}
	
		treeReader->ReadEntry(jentry);
		
		//Instatiate vector of non-btagged jets
		vector<Jet*> otrosjets;
		
		//declare and initilaze the event sizes for every branch, 
		// i.e. set of particles for this event
		Int_t photon_size = branchPhoton->GetEntriesFast();
		Int_t muon_size = branchMuon->GetEntriesFast();
		Int_t electron_size = branchElectron->GetEntriesFast();
		Int_t jet_size = branchJet->GetEntriesFast();
		
		//THIS IS THE CODE FOR MAKING THE B-CUTS
		//Initialize numBJets
        Int_t numBJets = 0;
		vector <Jet*>demJets;
        
        //Loop over jets to match the reconstructed and generated jets
        for (Int_t jet_entry=0; jet_entry<jet_size; jet_entry++) {
            Jet *jet = (Jet*) branchJet->At(jet_entry);
            // Cut jets with eta > 2.4
            if (fabs(jet->P4().Eta()) > 2.4) { continue;
            }
            //Cut jets with pt < 30
            if (jet->P4().Pt() < jet_pt) { continue;
            }
            //Find BTag using Loose tagger
            Bool_t Btag0k = (jet->BTag & (1 << 1));
            //Cut jets with Mass within 20 eta < 2.4
            if (Btag0k == 0) {
                // INSIDE HERE IS THE CODE FOR THE W-BOSON CUTS
                otrosjets.push_back(jet);
                    }
            //Look at Btagged jets
            if (Btag0k == 1) {
                //Add one to numBJets
                numBJets += 1;
                //Keep track of where the BJets are
                demJets.push_back(jet);
            } //End of Btag check
        } //End of jet loop
        
        // only continue if there are two b-jets in this event
        if (numBJets < 2) { continue;}

		less_than_two_Bs_cuts++;
		
		//calculate number of non btagged jets
        Int_t numOtros = otrosjets.size();
        //Create a bool to indicate whether event passes W cut
        Bool_t wcutfail = false;

		//Loop over pairs of non btagged jets
        for (Int_t jet1_no = 0; jet1_no < numOtros; jet1_no++) {

            for (Int_t jet2_no = jet1_no + 1; jet2_no < numOtros; jet2_no++) {

                //instantiate jet1 and jet 2
                Jet *jet1 = otrosjets[jet1_no];
                Jet *jet2 = otrosjets[jet2_no];
                
                //find parent mass
                Double_t nrg1, nrg2, px1, px2, py1, py2, pz1, pz2, nrg, px, py, pz;
                nrg1 = jet1->P4().Energy();
                nrg2 = jet2->P4().Energy();
                px1 = jet1->P4().Px();
                px2 = jet2->P4().Px();
                py1 = jet1->P4().Py();
                py2 = jet2->P4().Py();
                pz1 = jet1->P4().Pz();
                pz2 = jet2->P4().Pz();
                nrg = nrg1 + nrg2;
                px = px1 + px2;
                py = py1 + py2;
                pz = pz1 + pz2;
                Double_t parent_mass = sqrt(nrg*nrg - px*px - py*py - pz*pz);
                //check whether parent mass falls within 20GeV window of W
                if (parent_mass < wmass_min && parent_mass > wmass_max) {wcutfail = true; } //end of parent mass check
            } //end of jet2 loop
        } //end of jet1 loop
        
        if (wcutfail) { continue;} //cut the event if it fails the W mass check
        w_cuts++;

	//initialize bjet_index1 and 2
        Int_t BJet_index_1 = 0; //Should this be 1 or 0????
        Int_t BJet_index_2 = 1;

        //keep only the 2 most energetic bjets
        if (numBJets > 2) {
            Double_t ptmax1 = 0.0;
            Double_t ptmax2 = 0.0;
            for (Int_t index = 0; index < numBJets; index++) {
                Jet *bjet = demJets[index];
                Double_t bjetpt = bjet->P4().Pt();
                if (bjetpt > ptmax1) {
                    ptmax1 = bjetpt;
                    BJet_index_1 = index;
                }
            }

            for (Int_t index = 0; index < numBJets; index++) {
                Jet *Bjet = demJets[index];
                Double_t Bjetpt = Bjet->P4().Pt(); 
                if (Bjetpt > ptmax2) {
		  if (index == BJet_index_1) {continue; }
                  ptmax2 = Bjetpt;
                  BJet_index_2 = index;
                }
            }
        }

        else {
            BJet_index_1 = 0;
            BJet_index_2 = 1;
        }

        //instantiate the two bjets
        Jet *bjet1 = demJets[BJet_index_1];
        Jet *bjet2 = demJets[BJet_index_2];
        
        // create vector for higgs_on_shell_from_BB
        TLorentzVector Higgs_bb_vec = bjet1->P4() + bjet2->P4();
        
        //restrict deltaR to be greater than 0.4
        Double_t d1 = pow(bjet1->P4().Eta()-bjet2->P4().Eta(), 2);
        Double_t d2 = pow(bjet1->P4().Phi()-bjet2->P4().Phi(), 2);
        Double_t deltaR = sqrt(d1+d2);
        if (deltaR < bbdeltaR) { continue;}
		bb_isol_cuts++;
		
        //cacluate parent mass
        Double_t nrg1, nrg2, px1, px2, py1, py2, pz1, pz2, nrg, px, py, pz;
        nrg1 = bjet1->P4().Energy();
        nrg2 = bjet2->P4().Energy();
        px1 = bjet1->P4().Px();
        px2 = bjet2->P4().Px();
        py1 = bjet1->P4().Py();
        py2 = bjet2->P4().Py();
        pz1 = bjet1->P4().Pz();
        pz2 = bjet2->P4().Pz();
        nrg = nrg1 + nrg2;
        px = px1 + px2;
        py = py1 + py2;
        pz = pz1 + pz2;
        Double_t parent_mass = sqrt(nrg*nrg - px*px - py*py - pz*pz);
        
        //restrict parent mass to be within 25 GeV window of higgs
        parent_mass_hist_bb->Fill(parent_mass);
        
        if (parent_mass > bbmass_max || parent_mass < bbmass_min) { continue;}
		
		bb_mass_cuts++;

		// if it makes it this far, then now it just has to go through the gm gm cuts	
			 
			 Int_t leptons_bool = 1;
			 
			 // loop through muons
			 for (int p=0; p < muon_size; p++){
				 
				 // create isntance of particle
				 Muon *muon = (Muon*)branchMuon->At(p);
				
				if ((muon->PT > lepton_pt) && (abs(muon->Eta) < 2.4)){
						leptons_bool = 0;
					}
				 }
			 
			 // loop through the electrons
			 for (int w = 0; w < electron_size; w++){
				 
				 // create instance of electron
				 Electron *electron = (Electron*)branchElectron->At(w);
				 
				 if ((electron->PT > lepton_pt) && (abs(electron->Eta) < 2.4)){
						leptons_bool = 0;
					 }
				 }
			 
			// the lepton cuts
			if (leptons_bool == 0){continue;} // skips event if it doesn't make the lepton cut
			lepton_cuts++;

			Int_t index_of_first = 1;
			Int_t index_of_second = 0; //should this be 1 or 0????
			// create vector 
			vector <Photon*>luxons;
			
			Int_t numPhotons = 0;
			// loop through photons
			for (int r = 0; r < photon_size; r++){
				Photon *photon = (Photon*)branchPhoton->At(r);
				if (photon->PT < photon_pt || abs(photon->Eta > 2.4)){continue;}
				luxons.push_back(photon);
				numPhotons++;
			}
								
			if (numPhotons < 2){continue;}
			// cut for phtons < 2
			less_than_two_Gms_cuts++;
			
			if (numPhotons == 2){
				index_of_first = 0;
				index_of_second = 1;
				}
			
			int photon_PTs[numPhotons];	
			//keep only the 2 most energetic photons
			if (numPhotons > 2) {
				
				Double_t ptmax1 = 0.0;
				for (Int_t index = 0; index < numPhotons; index++) {
					Photon *phot = luxons[index];
					Double_t phot_pt = phot->P4().Pt();
					photon_PTs[index] = phot_pt;
					if (phot_pt > ptmax1) {
						ptmax1 = phot_pt;
						index_of_first = index;
					}
				}
			
				Double_t ptmax2 = 0.0;
				for (Int_t index2 = 0; index2 < numPhotons; index2++) {
					Photon *Phot = luxons[index2];
					Double_t Phot_pt = Phot->P4().Pt(); 
					if ((Phot_pt > ptmax2) && (index2 != index_of_first)) {
						ptmax2 = Phot_pt;
						index_of_second = index2;
					}
				}
			}
			
			// create the two photons to be used next
			Photon *photon1 = luxons[index_of_first];
			Photon *photon2 = luxons[index_of_second];

			// create vector for higgs_on_shell_from_BB
			TLorentzVector Higgs_gmgm_vec = photon1->P4() + photon2->P4();

			/*if (numPhotons >= 3){
				cout << "there are: " << numPhotons << " in this event" << endl;
				cout << "index_of_first: " << index_of_first << endl;
				cout << "index_of_second: " << index_of_second << endl;
				
				for (int i = 0; i < numPhotons; i++){

					cout << "this is photon with index: " << i << "and PT: " << photon_PTs[i] << endl;
				}
					
				//cout << "press enter to continute: " << endl;
				//cin.get();
				
				}*/
			
			 // find the mass of these photons's mothers
			 double nrg_sum = 0, px_sum = 0, py_sum = 0, pz_sum = 0;
					
			// create a lorentz vector for photon 1 
			Double_t phot_nrg1,phot_px1,phot_py1,phot_pz1;
			phot_nrg1=photon1->E;
			phot_px1=photon1->P4().Px(); 
			phot_py1=photon1->P4().Py();
			phot_pz1=photon1->P4().Pz();
					
			// create a lorentz vector fotr photon 2
			Double_t phot_nrg2,phot_px2,phot_py2,phot_pz2;
			phot_nrg2=photon2->E;
			phot_px2=photon2->P4().Px(); 
			phot_py2=photon2->P4().Py();
			phot_pz2=photon2->P4().Pz();
			nrg_sum = phot_nrg1+phot_nrg2, px_sum = phot_px1+phot_px2, py_sum = phot_py1+phot_py2, pz_sum = phot_pz1+phot_pz2;
			
			Double_t mass_parent = sqrt((nrg_sum)*(nrg_sum)-(px_sum)*(px_sum)-(py_sum)*(py_sum)-(pz_sum)*(pz_sum));
			if (mass_parent > ggmass_max || mass_parent < ggmass_min){continue;} // skip the event if it's not within range of mass
			gmgm_mass_cuts++;
			
			// find the deltaR between two photons
					
			Double_t phi_term = pow(photon1->Phi-photon2->Phi,2);
			Double_t eta_term = pow(photon1->Eta-photon2->Eta,2);
			Double_t delR = sqrt(phi_term+eta_term);
			//cout << "Eta_hist: "/ 
			if (delR < ggdeltaR){continue;}
			gmgm_isol_cuts++;
			
			//find the deltaR between each photon and each bjet
			Bool_t failed_the_isolation = false;
			for (Int_t j= 0; j < numPhotons; j++) {
				if (j != index_of_first && j != index_of_second) { continue; }
				Photon* disphoton = luxons[j];
				for (Int_t t = 0; t < numBJets; t++) {
					if (t != BJet_index_1 && t != BJet_index_2) { continue; }
					Jet* disjet = demJets[t];
					Double_t deltaPhi = pow(disphoton->Phi-disjet->Phi, 2);
					Double_t deltaEta = pow(disphoton->Eta-disjet->Eta, 2);
					Double_t deltaArrrr = sqrt(deltaPhi + deltaEta);
					if (deltaArrrr < bgdeltaR) { 
						failed_the_isolation = true;
						} //end of deltaArrrr check
					} //end of bjet loop
				} //end of photon loop
			if (failed_the_isolation) { continue; }
			bgm_isol_cuts++;		

		        TLorentzVector bbgg_vec = Higgs_bb_vec + Higgs_gmgm_vec;
			
			// find PT of the offshell higgs
			Double_t Higgs_off_shell_PT = bbgg_vec.Pt();	
			
			// make cuts on Higgs PT accordingly
			if (Higgs_off_shell_PT < bbgg_pt){continue;} // gets cut
			Higgs_PT_cuts++;
			
			// find the Eta of the offshell higgs
			Double_t Higgs_off_shell_Eta = bbgg_vec.Eta();
			
			// make cuts on Higgs Eta accordingly
			if (abs(Higgs_off_shell_Eta) > bbgg_eta){continue;}
			Higgs_Eta_cuts++;
			//cout << "higgs_off_shell_Eta: " << abs(Higgs_off_shell_Eta) << endl;
			
			// find the mass of the diHiggs
			/*Double_t higgs_nrg_bb,higgs_px_bb,higgs_py_bb,higgs_pz_bb;
			higgs_nrg_bb=Higgs_bb_vec.E();
			higgs_px_bb=Higgs_bb_vec.Px(); 
			higgs_py_bb=Higgs_bb_vec.Py();
			higgs_pz_bb=Higgs_bb_vec.Pz();
			
			Double_t higgs_nrg_gmgm,higgs_px_gmgm,higgs_py_gmgm,higgs_pz_gmgm;
			higgs_nrg_gmgm=Higgs_gmgm_vec.E();
			higgs_px_gmgm=Higgs_gmgm_vec.Px(); 
			higgs_py_gmgm=Higgs_gmgm_vec.Py();
			higgs_pz_gmgm=Higgs_gmgm_vec.Pz();*/
		
			Double_t nrg_sum2 = 0, px_sum2 = 0, py_sum2 = 0, pz_sum2 = 0;

			nrg_sum2 = bbgg_vec.E();
			px_sum2 = bbgg_vec.Px();
			py_sum2 = bbgg_vec.Py();
			pz_sum2 = bbgg_vec.Pz();
						
			Double_t diHiggs_mass = sqrt((nrg_sum2)*(nrg_sum2)-(px_sum2)*(px_sum2)-(py_sum2)*(py_sum2)-(pz_sum2)*(pz_sum2)) ; // parent_mass is from bb, mass_parent is from gmgm
			
			//cout << "dihiggs mass cuts: " << diHiggs_mass << endl;
			
			// make the cuts on diHiggs mass accordingly
			if (diHiggs_mass < bbgg_mass){continue;}
			diHiggs_mass_cuts++;
		
		parent_mass_hist_gmgm->Fill(mass_parent);

				
	} // end of event loop
			
			// make canvas, draw histograms of mass
			cout << "hi: " << endl;
			//TCanvas *canvas = new TCanvas("canvas","Mass of Higgs",10,10,900,600);
			cout << "hi2: " << endl;
			//canvas->Divide(2,1);
			//canvas->cd(1);
			parent_mass_hist_bb->Draw();
			//canvas->cd(2);
			parent_mass_hist_gmgm->Draw();
			
			Double_t cross_section = 0;
			Int_t luminosity = 3000;
			// print summary of cut results
			cout << "*******************************************************" << endl;
			cout << " SUMMARY FOR SAMPLE TYPE: " << endl;
			if (c == 'h'){cross_section = .0892; cout << "SIGNAL: gg -> HH --> bbgmgm" << endl;}
			if(c == 't'){cross_section = 1.39; cout << "BACKGROUND: gg --> ttH --> bbWWgmgm" << endl;}
			if(c == 'z'){cross_section = .333; cout << "BACKGROUND: gg --> ZH --> bbgmgm" << endl;}
			cout << "*******************************************************" << endl;
			cout << " Number of events: "  << '\t' << '\t'  << '\t' << numberOfEntries<< endl;
			cout << " Made the less than two b quarks cuts: " << '\t' <<  less_than_two_Bs_cuts << endl;
			cout << " Made the W-cuts: " << '\t' << '\t'  << '\t' << w_cuts << endl;
			cout << " Made the bb isolation cuts: " << '\t' << '\t' << bb_isol_cuts << endl;
			cout << " Made the bb parent mass cuts: "<< '\t'  << '\t' << bb_mass_cuts << endl;
			cout << " Made the lepton cuts: " << '\t'<< '\t'  << '\t' << lepton_cuts << endl;
			cout << " Made the less than two photons cuts: " << '\t' << less_than_two_Gms_cuts << endl;
			cout << " Made the gmgm parent mass cuts: " << '\t' << gmgm_mass_cuts << endl;
			cout << " Made the gmgm isolation cuts: " << '\t' << '\t' << gmgm_isol_cuts << endl;
			cout << " Made the b-gm isolation cuts: " << '\t' << '\t' <<  bgm_isol_cuts << endl;
			cout << "------------------------------------------------------" << endl;
			cout << " ADVANCED CUTS: " << endl;
			cout << "------------------------------------------------------" << endl;
			cout << " Made the Higgs PT cuts: " << '\t' << '\t'  << Higgs_PT_cuts << endl;
			cout << " Made the Higgs Eta cuts: " << '\t' << '\t' <<  Higgs_Eta_cuts << endl;
			cout << " Made the DiHiggs M_{HH} Mass cuts: " << '\t' << diHiggs_mass_cuts << endl;
			cout << "------------------------------------------------------" << endl;
			cout << " CALCULATING NUMBER OF EVENTS" << endl;
			cout << "------------------------------------------------------" << endl;
			Double_t e_times_a = (double)diHiggs_mass_cuts/numberOfEntries;
			cout << " With advanced cuts: " << endl;
			cout << '\t' << "Luminosity: " << '\t' << '\t' << '\t' <<luminosity << " fb"<< endl;			
			cout << '\t' << "Cross section: " << '\t' << '\t' << '\t' << cross_section << endl;
			cout << '\t' << "Efficieny * Acceptance: " << '\t'  << e_times_a << endl;
			if (c == 'h'){cout << '\t' << "S";}
			if (c == 't'){cout << '\t' << "B";}
			if (c == 'z'){cout << '\t' << "B";}

			cout << '\t' << '\t' << '\t'  << '\t' << e_times_a * luminosity * cross_section << endl << endl;
			
			e_times_a = (double)bgm_isol_cuts/numberOfEntries;
			cout << " Without advanced cuts: " << endl;
			cout << '\t' << "Luminosity: " << '\t' << '\t' << '\t' <<luminosity << " fb"<< endl;			
			cout << '\t' << "Cross section: " << '\t' << '\t' << '\t' << cross_section << endl;
			cout << '\t' << "Efficieny * Acceptance: " << '\t'  << e_times_a << endl;
			cout << '\t' << "B: " << '\t' << '\t' << '\t'  << '\t' << e_times_a * luminosity * cross_section << endl;
 			cout << "*******************************************************" << endl;

			return e_times_a * luminosity * cross_section;

} // end of function

/*
void optimalSignificance(){

  Double_t significance_max = 0;
  //using boost::assign::map_list_of;
  
  // create map of cuts and initialize all the cuts to zero (they will soon be tested at the appropriate ranges)
  //map<string,int> cuts = map_list_of("jet_pt",0)("wmass_width",0)("bbdeltaR",0)("bbmass_width",0)("lepton_pt",0)("photon_pt",0)("ggmass_width",0)("ggdeltaR",0)("bgdeltaR",0)("bbgg_pt",0)("bbgg_eta",0)("bbgg_mass",0); // change the last one here

  map<string,int> MapOfCuts;
  // add the cuts to the map
  MapOfCuts.insert(std::pair<string,int>("jet_pt",0)); 
  MapOfCuts.insert(std::pair<string,int>("wmass_width",0)); 
  MapOfCuts.insert(std::pair<string,int>("bbdeltaR",0)); 
  MapOfCuts.insert(std::pair<string,int>("bbmass_width",0)); 
  MapOfCuts.insert(std::pair<string,int>("lepton_pt",0)); 
  MapOfCuts.insert(std::pair<string,int>("photon_pt",0)); 
  MapOfCuts.insert(std::pair<string,int>("ggmass_width",0)); 
  MapOfCuts.insert(std::pair<string,int>("ggdeltaR",0)); 
  MapOfCuts.insert(std::pair<string,int>("bgdeltaR",0)); 
  MapOfCuts.insert(std::pair<string,int>("bbgg_pt",0)); 
  MapOfCuts.insert(std::pair<string,int>("bbgg_eta",0)); 
  MapOfCuts.insert(std::pair<string,int>("bbgg_mass",0)); 

// begin nested for loops of 12 cuts
for (Int_t a = 20; a <= 40; a += 5) // jet_pt > 30 +/- 10 ---------------------------------- (1)
  {
    for (Int_t b = -10; b <= 10; b+= 5) // wmass_width loop =  80 +/- 10 --------------------- (2)
    {
      for (Int_t c = 0.3; b <= 0.5; c+= 0.1) // bbdeltaR < 0.4 +/- 0.1 ----------------------- (3) 
	{
	  for (Int_t d = 0; d <= 20 ; d+=5) // bbmass_width = 130 +/- 20 --------------------- (4) 
	    {
	      for (Int_t e = 20; e <= 40 ; a += 5) // lepton_pt > 30 +/- 10 ------------------ (5)
		{
		  for (Int_t f = 20 ; f <= 40; f += 5) // photon_pt > 30 +/- 10 -------------- (6) 
		    {
		      for (Int_t g = 0; g <= 5; g +=5) // ggmass_width = 125 +/- 5 ----------- (7) 
			{
			  for (Int_t h = 0.3; g <= 0.5; h += 0.1) // ggDeltaR < 0.4 +/- 0.1 -- (8)
			    {
			      for (Int_t i = 0.3; i < 0.5; i += 0.1) // bgDeltaR ------------- (9)
				{
				  for (Int_t j = 85; j < 115 ; j += 5) //bbgg_pt = 100 +/- 15 (advanced cut) -- (10)
				    {
   				      for (Int_t k = 1.0; k < 3.0; k += 0.2) //bbgg_eta = 2 +/- 1 (advanced cut) -- (11)
					{	
					  for (Int_t l = 340; l <= 360 ; l += 5) //bbgg_mass_width = 350 +/- 10 (advanced cut) -- (12)
					   {
					    // calculate the new significance in this configuration of cuts
					     Double_t significance_new = numEvents('h',MapOfCuts)/(sqrt(numEvents('z',MapOfCuts) + numEvents('t',MapOfCuts) + 0.5)); // significance = s / sqrt(b + 0.5)
					    
					    // if the new significance is higher, than adjust the values of the cuts accordingly in the map and also update the significance_max
					    if (significance_new > significance_max)
					    {
					      MapOfCuts["jet_pt"] = a;
					      MapOfCuts["wmass_width"] = b;
					      MapOfCuts["bbdeltaR"] = c;
					      MapOfCuts["bbmass_width"] = d;
					      MapOfCuts["lepton_pt"] = e;
					      MapOfCuts["photon_pt"] = f;
					      MapOfCuts["ggmass_width"] = g;
					      MapOfCuts["ggdeltaR"] = h;
					      MapOfCuts["bgdeltaR"] = i;
					      MapOfCuts["bbgg_pt"] = j;
					      MapOfCuts["bbgg_eta"] = k;
					      MapOfCuts["bbgg_mass"] = l;
					      significance_max = significance_new;
					    }
					  }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  }

 cout << "The higest possible significance on these samples is: " << significance_max;

}*/
