#include <iostream>
#include "sigReader.h"
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
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "TMathBase.h"
#include "TLatex.h"
#include "TString.h"
#include "TChain.h"

// create the histograms

void cuts(){
	// create a chain of root files
	//TChain chain("Delphes");
	//chain.Add(s);
	//chain.Add(s);
	TChain *chain = new TChain("Delphes");
	chainMaker(chain);
	
	//create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	Long64_t numberOfEntries = treeReader->GetEntries();
	
	// get pointers to branches that will be used in this macro
	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	
	// create histogram to eventually plot the events that make it past the cuts
	TH1F *mass_phot_pair_sg = new TH1F("numEvents", "Number of Events that Made Cut", 140, 0, 500);
	mass_phot_pair_sg->GetXaxis()->SetTitle("M_{gmgm}");
	mass_phot_pair_sg->GetYaxis()->SetTitle("# Events");	
	
	TH1F *lepton_PT_hist = new TH1F("leptonPT", "PT of Leptons", 140, 0, 500);
	TH1F *lepton_Eta_hist = new TH1F("leptonEta", "Eta of Leptons", 140, -7, 7);
	TH1F *photon_Eta_hist = new TH1F("photonEta", "Eta of Photons", 140,-7,7);
	TH1F *photon_delR_hist = new TH1F("photondelR", "delR btwn photons",50,0,1);
	TH1F *num_photons = new TH1F("numPhoton", "num_photons in event", 100, 0, 5);
	
	TH1F *parent_mass_hist = new TH1F("parent_mass", "apparent parent mass", 200, 0, 2000);
	parent_mass_hist->GetXaxis()->SetTitle("M_{bb}");
	parent_mass_hist->GetYaxis()->SetTitle("# events");	

	Int_t num_that_made_matt_cuts = 0;	
	Int_t num_gotten_this_far = 0;
	
	// the event loop
	for (Int_t jentry = 0; jentry < numberOfEntries; jentry++){

		treeReader->ReadEntry(jentry);
		
		//Instatiate vector of non-btagged jets
		vector<Jet*> otrosjets;
		
		/*declare and initilaze the event sizes for every branch, 
		 * i.e. set of particles for this event*/
		Int_t photon_size = branchPhoton->GetEntriesFast();
		Int_t muon_size = branchMuon->GetEntriesFast();
		Int_t electron_size = branchElectron->GetEntriesFast();
        Int_t jet_size = branchJet->GetEntriesFast();
		
		/*THIS IS THE CODE FOR MAKING THE B-CUTS*/
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
            if (jet->P4().Pt() < 30.0) { continue;
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
        if (numBJets < 2) { continue;
			}
		cout << "this is event: " << jentry << "and there are at least two bjets in this event!" << endl;

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
                cout << parent_mass << endl;
                //check whether parent mass falls within 20GeV window of W
                if (parent_mass < 90.0 && parent_mass > 70.0) {
                    wcutfail = true;
                } //end of parent mass check
            } //end of jet2 loop
        } //end of jet1 loop
        //cut the event if it fails the W mass check
        if (wcutfail) {
            continue;
        }
        cout << "this is event: " << jentry << "and I made the w-cuts" << endl;
		//initialize bjet_index1 and 2
        Int_t BJet_index_1 = 0;
        Int_t BJet_index_2 = 0;
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
					if (index = BJet_index_1) {continue; }
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
        
        //restrict deltaR to be greater than 0.4
        Double_t d1 = pow(bjet1->P4().Eta()-bjet2->P4().Eta(), 2);
        Double_t d2 = pow(bjet1->P4().Phi()-bjet2->P4().Phi(), 2);
        Double_t deltaR = sqrt(d1+d2);
        if (deltaR < 0.4) { continue;
        }
        cout << "this is event: " << jentry << "and I made matt's deltaR (isolation)" << endl;
        
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
        parent_mass_hist->Fill(parent_mass);
        if (parent_mass > 137.5 || parent_mass < 112.5) { continue;
        }
						
				cout << "this is event: " << jentry << "and I made the Matt's mass cuts" << endl;
				cout << "this is event: " << jentry << "and I made ALL OF MATT'S CUTS" << endl << endl;
				num_that_made_matt_cuts++;		
				num_gotten_this_far++;
		cout << "num that have gotten this far: " << num_gotten_this_far << endl;
		cout << "prese enter tocontinute: " << endl;
		//cin.get();
		// if it makes it this far, then now it just has to go through the gm gm cuts	
			 
			 Int_t leptons_bool = 1;
			 
			 // loop through muons
			 for (int p=0; p < muon_size; p++){
				 
				 // create isntance of particle
				 Muon *muon = (Muon*)branchMuon->At(p);
				
				if ((muon->PT > 20) && (abs(muon->Eta) < 2.4)){
						leptons_bool = 0;
					}
					
			 }
			 
			 // loop through the electrons
			 for (int w = 0; w < electron_size; w++){
				 
				 // create instance of electron
				 Electron *electron = (Electron*)branchElectron->At(w);
				 
				 if ((electron->PT > 20) && (abs(electron->Eta) < 2.4)){
						leptons_bool = 0;
					 }
				 }
			 
				    
			 // the lepton cuts
			 cout << "leptons_bool: " << leptons_bool << endl;
			 if (leptons_bool == 0){continue;} // skips event if it doesn't make the lepton cut
			 cout << "this is event: " << jentry << "and I made Hussein's lepton cuts" << endl;
			 
				    
			Int_t index_of_first = 0;
			Int_t index_of_second = 0;
			// create vector 
			vector <Photon*>luxons;
			
			Int_t numPhotons = 0;
			// loop through photons
			for (int r = 0; r < photon_size; r++){
				Photon *photon = (Photon*)branchPhoton->At(r);
				if (photon->PT < 30.0 || abs(photon->Eta > 2.4)){continue;}
				luxons.push_back(photon);
				numPhotons++;
			}
			
			
			num_photons->Fill(numPhotons);
		
			if (numPhotons < 2){continue;}
		
		
			//keep only the 2 most energetic photons
			if (numPhotons > 2) {
				Double_t ptmax1 = 0.0;
				Double_t ptmax2 = 0.0;
				for (Int_t index = 0; index < numPhotons; index++) {
					Photon *phot = luxons[index];
					Double_t phot_pt = phot->P4().Pt();
					if (phot_pt > ptmax1) {
						ptmax1 = phot_pt;
						index_of_first = index;
					}
				}
				
				for (Int_t index = 0; index < numPhotons; index++) {
					Photon *Phot = luxons[index];
					Double_t Phot_pt = Phot->P4().Pt(); 
					if (Phot_pt > ptmax2) {
						if (index = index_of_first) {continue; }
						ptmax2 = Phot_pt;
						index_of_second = index;
					}
				}
			}
		
			// create the two photons to be used next
			//cout << "index_of_first: " << index_of_first << endl;
			//cout << "index_of_second: " << index_of_second << endl;
			Photon *photon1 = luxons[index_of_first];
			Photon *photon2 = luxons[index_of_second];
			cout << "I got here and the photons were created!" << endl;
			cout << "photon2-.P4().Px(): " << photon2->P4().Px() << endl;
			 // find the mass of these photons's mothers
			 double nrg_sum = 0, px_sum = 0, py_sum = 0, pz_sum = 0;
					
			
			// create a lorentz vector for photon 1 
			Double_t phot_nrg1,phot_px1,phot_py1,phot_pz1;
			phot_nrg1=photon1->E;
			phot_px1=photon1->P4().Px(); 
			phot_py1=photon1->P4().Py();
			phot_pz1=photon1->P4().Pz();
					
					cout << "the first loretnz vector was created" << endl;
			// create a lorentz vector fotr photon 2
			Double_t phot_nrg2,phot_px2,phot_py2,phot_pz2;
			phot_nrg2=photon2->E;
			phot_px2=photon2->P4().Px(); 
			cout << "im here ran out of messages" << endl;
			phot_py2=photon2->P4().Py();
			phot_pz2=photon2->P4().Pz();
					cout << "the second lorentz vectow was created" << endl;
			nrg_sum = phot_nrg1+phot_nrg2, px_sum = phot_px1+phot_px2, py_sum = phot_py1+phot_py2, pz_sum = phot_pz1+phot_pz2;
			
			cout << "I found the energy sum" << endl;
			
			Double_t mass_parent = sqrt((nrg_sum)*(nrg_sum)-(px_sum)*(px_sum)-(py_sum)*(py_sum)-(pz_sum)*(pz_sum));
			if (parent_mass > 130 || parent_mass < 120){continue;} // skip the event if it's not within range of mass
			cout << "this is event: " << jentry << "and I made Husseins invariant mass cuts and yes i believe mass changes" << endl;

				// find the momemntum
				
					Double_t phot_PT = photon1->PT;
					Double_t phot_PT2 = photon2-> PT;
					if ((phot_PT < 30) || (phot_PT2 < 30)){continue;} // skip the event if it's not within momentum
					        cout << "this is event: " << jentry << "and I made Husseins momentum cuts" << endl;

					cout << "phot_PT: " << phot_PT << endl;
					cout << "phot_PT2: " << phot_PT2 << endl;
				
				// find the eta
				
					Double_t phot_eta = photon1->Eta;
					Double_t phot_eta2 = photon2->Eta;
					if ((abs(phot_eta) > 2.4) || (abs(phot_eta2) > 2.4)){continue;} 
					        cout << "this is event: " << jentry << "and I made Husseins eta cuts" << endl;
					
				    
					cout << "phot_eta: " << phot_eta << endl;
					cout << "phot_eta2: " << phot_eta2 << endl;
				
				
				// find the deltaR between two photons
					
					Double_t phi_term = pow(photon1->Phi-photon2->Phi,2);
					Double_t eta_term = pow(photon1->Eta-photon2->Eta,2);
					Double_t delR = sqrt(phi_term+eta_term);
					if (delR < 0.4){continue;}
					        cout << "this is event: " << jentry << "and I made Husseins isolation cuts " << endl;
						
						//cout << "I've made all of Matt's cuts, hit enter to continue: " << endl;
				     
					
					cout << "delR: " << delR << endl;
					
						cout << "this is event: " << jentry << "and I made all of the cuts" << endl;
				        cout << "num_that_made_matt_cutts: " << num_that_made_matt_cuts << endl;
				        cout << "I've made all of the cuts cuts, hit enter to continue: " << endl;


					// if it's made it this far, then it's made all cuts. go ahead and plot the mass
				// if it passes the cuts, plot the invariant mass for that event
		cout << "----------------------------------------------------------" << endl << endl;

					mass_phot_pair_sg->Fill(mass_parent);

				
	} // enf of event loop
				cout << "num that made matts cuts: " << num_that_made_matt_cuts << endl;
				cout << "Draw me like one of your french girls" << endl;
				cout << "num that got this far: " << num_gotten_this_far << endl;
			//mass_phot_pair_sg->Draw();
			//lepton_PT_hist->Draw();
			//lepton_Eta_hist->Draw();
				//num_photons->Draw();
				TCanvas *canvas = new TCanvas("canvas","hi",10,10,1200,700);
				canvas->Divide(2,1);
				canvas->cd(1);
				parent_mass_hist->Draw();
				canvas->cd(2);
				mass_phot_pair_sg->Draw();
	
} // end of function

