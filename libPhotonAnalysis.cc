/*libPhotonAnalysis is a library of functions which can be used to analyze the performance of photons in samples. 
This library is designed to help you validate your Delphes samples.*/

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
#include <string>
#include <string.h>
#include <iostream>
#include <map>
#include <utility>
#include <stdbool.h>
using namespace std;

/*This function calculates the photon reconstruction effi*/
void photonReconstructionEfficiency(TString s){
	// create a chain of root files
	TChain *chain = new TChain("Delphes");

	cout << "Calculating Reconstruction Efficiency ";
	// if char c = s this means running cuts on signal
	if (s == "h"){
	  cout << "for HH signal sample." << endl;
	  chainMakerHH(chain);
	}
	
	// if char c = b this means runnign cuts on the background m
	else if (s == "t" ){
	  cout << "for ttH background sample"<< endl;
	  chainMakerttH(chain);
	}
		
	else if (s == "z"){
	  cout << "for ZH background sample "<< endl;
	  chainMakerZH(chain);
	}
		
	//create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	Long64_t numberOfEntries = treeReader->GetEntries();
	Long64_t numEvents = numberOfEntries;
	
	// get pointers to branches that will be used in this macro
	TClonesArray *branchRecPhoton = treeReader->UseBranch("Photon");
	TClonesArray *branchGenPhoton = treeReader->UseBranch("Particle");
	
	// declare variables holding amounts of filtered and unfiltered rec/gen photons
	Int_t numGenPhotonsFiltered = 0, numGenPhotonsUnfiltered = 0, numRecPhotonsFiltered = 0, numRecPhotonsUnfiltered = 0;

	// begin the event loop

	for (Int_t jentry = 0; jentry < numberOfEntries; jentry++){
		
		if (jentry % 5000 == 0){cout << "processing event no: " << jentry << endl;}
		
		// read the event
		treeReader->ReadEntry(jentry);

		// get the raw (unfiltered) number of rec and gen photons in this event
		Int_t rec_photon_size = branchRecPhoton->GetEntriesFast();
		Int_t gen_particle_size = branchGenPhoton->GetEntriesFast();

		numGenPhotonsUnfiltered += gen_particle_size; // add the num of gen photons to total unfiltered amount
		numRecPhotonsUnfiltered += rec_photon_size; // add the num of rec photons to total unfiltered amount
		
		// loop through the reconstructed photons and match them to generated photons
		
		// declare all the generated photons in this event to be unmatched
		bool MatchedGenPhotons[2] = {false};

		for (Int_t r = 0; r < rec_photon_size; r++){ // loop through the reconstructed photons in this event
 			
 			Photon *rec_photon = (Photon*)branchRecPhoton->At(r); // create instance of rec photon
 			
 			/*	don't considerthis photon if it has a low PT or is not within fiducial volume
 				the point of this is to eliminate photons which did not come from the Higgs
 				and ideally the cut would look at the actual parents using the M1 and M2 attributes 
 				of the photons, but until that is figured out this filter will take its place
 			*/ 
 			if (rec_photon->PT < 20 || abs(rec_photon->Eta) > 2.4 ){continue;} 
 			numRecPhotonsFiltered++; 
 			Int_t index_of_closest_gen_phot = 0;
			
			for (Int_t g = 0; g < gen_particle_size; g++){ // loop through the generated photons in this event

				GenParticle *gen_particle = (GenParticle*) branchGenPhoton->At(g); // create gen particle instance

				if (gen_particle->PID != 22){continue;} // skip this gen particle if it's not a gen photon
				/*	don't considerthis photon if it has a low PT or is not within fiducial volume
 				the point of this is to eliminate photons which did not come from the Higgs
 				and ideally the cut would look at the actual parents using the M1 and M2 attributes 
 				of the photons, but until that is figured out this filter will take its place
 				*/ 
				if (gen_particle->PT < 20 || abs(gen_particle->Eta) > 2.4 ){continue;}
				numGenPhotonsFiltered++;

				if (MatchedGenPhotons[g] == true){continue;} // skip this gen photon if it's already been matched in this event

				// calculate the distance in the eta-phi plane between this reconstructed and generated photon
				Double_t phiSquared = pow(rec_photon->Phi - gen_particle->Phi,2);
				Double_t etaSquared = pow(rec_photon->Eta - gen_particle->Eta,2);
				Float_t distance = sqrt(phiSquared+etaSquared);

				/* if the distance is less than 0.3, i.e. only check for the matching within the same cone on 
				a given reconstructed photon*/
				Int_t min_delta_PT = 5000; // set a low standard to meet for the first time (note: this is kind of hacky)
				if (distance < 0.3) {
						
					// now try to find the one with the closest Photon_PT value (concept of a cone)
					Double_t current_delta_PT = abs(rec_photon->PT - gen_particle->PT); // find this distance		

					if (current_delta_PT < min_delta_PT){
						// set the index of said one with closest PT value 
						index_of_closest_gen_phot = g;

						/*
						// if we've gotten this far then the reconstructed photon was matched
						// so add it to the total number of reconstructed photons observed in the sample
						// also important to note that we can only incremement the total here b/c we've
						// confirmed that the photon is legit (it's been matched to a generated photon
						// in the eta-phi plane and has a similar PT) and not just a random photon						// photon from 
						*/
						min_delta_PT = current_delta_PT;	
					
					} // end of delta_PT comparison
					
				} // end of eta if statement
			
			} // end of gen_particle loop
		
			// mark this gen photon as matched for this event						
			MatchedGenPhotons[index_of_closest_gen_phot] = true; 

		} // end of rec_photon loop
			
	} // end of photon loop
	
	// display results in console (eventually should display these results in a pad)
	cout << "*************************************************************************" << endl;
	cout << "PHOTON RECONSTRUCTION EFFICIENCY SUMMARY: " << setw(24) << numEvents << " events" << endl;
	cout << "-------------------------------------------------------------------------" << endl;
	cout << "Predicted Results (expected)" << endl;
	cout << "-------------------------------------------------------------------------" << endl;
	cout << "Expected # of generated photons for this sample (2 * numEvents): " << 2 * numEvents << endl;
	cout << "Expected # of reconstructed photons for this sample: " << '\t' << '\t' << 2 * numEvents << endl;
	cout << "-------------------------------------------------------------------------" << endl;
	cout << "Actual Results" << endl;	
	cout << "-------------------------------------------------------------------------" << endl;
	cout << "Unfiltered amount of generated photons: " << '\t' << '\t' << '\t' << numGenPhotonsUnfiltered << endl;
	cout << "Unfiltered amount of reconstructed photons: " << '\t'<< '\t' << '\t' << numRecPhotonsUnfiltered << endl;
	cout << "Filtered amount of generated photons for this sample: " << '\t' << '\t' << numGenPhotonsFiltered << endl;
	cout << "Filtered amount of reconstructed photons for this sample: "  << '\t' << numRecPhotonsFiltered << endl;
	cout << "-------------------------------------------------------------------------" << endl;
 	cout << "\033[1;34m";
	cout << "Photon Reconstruction Efficiency: " << '\t' << '\t' << '\t' << '\t' << ((float)numRecPhotonsFiltered/numGenPhotonsFiltered)*100<< "\%";	
 	cout << "\033[0m\n";
	cout << "-------------------------------------------------------------------------" << endl;
	cout << "*************************************************************************" << endl;
}

void photonResolution(TString s){


}

/*This function creates and displays two canvases each with two histograms with the distribution of the number of generated photons
per event. A couple of histograms will take into account ALL photons (including those that come from noisy processes
like pions) while the other two will filter for a certain PT and fiducial volume attributes to filter out the photons
which did not come from a Higgs.
CANVAS 1: (left) Unfiltered num of generated photons (right) filtered num of generated photon
CANVAS 2: (left) Unfiltered num of reconstructed photons (right) filtered num of reconstructed photon
 */
void photonCountDistribution(TString s){
// create a chain of root files
	TChain *chain = new TChain("Delphes");

	cout << "Finding the Photon Count Distribution ";
	// if char c = s this means running cuts on signal
	if (s == "h"){
	  cout << "for HH signal sample." << endl;
	  chainMakerHH(chain);
	}
	
	// if char c = b this means runnign cuts on the background m
	else if (s == "t" ){
	  cout << "for ttH background sample"<< endl;
	  chainMakerttH(chain);
	}
		
	else if (s == "z"){
	  cout << "for ZH background sample "<< endl;
	  chainMakerZH(chain);
	}
	
	// create histograms to display photon count distribution per event
	TH1F *gen_photon_distribution_hist = new TH1F("gen_photon_distribution_hist", "Number of Generated Photons in an Event", 5, 0, 5);
	TH1F *rec_photon_distribution_hist = new TH1F("rec_photon_distribution_hist", "Number of Reconstructed Photons in an Event", 5, 0, 5);

	gen_photon_distribution_hist->GetYaxis()->SetTitle("# events");	
	gen_photon_distribution_hist->GetXaxis()->SetTitle("# Gen Photons");

	rec_photon_distribution_hist->GetYaxis()->SetTitle("# events");	
	rec_photon_distribution_hist->GetXaxis()->SetTitle("# Rec Photons");

	//create object of class ExRootTreeReader
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	Long64_t numberOfEntries = treeReader->GetEntries();
	Long64_t numEvents = numberOfEntries;
	// get pointers to branches that will be used in this macro
	TClonesArray *branchGenPhoton = treeReader->UseBranch("Particle");
	TClonesArray *branchRecPhoton = treeReader->UseBranch("Photon");
	
	// declare variables
	Int_t actualNumGenPhotons = 0;
	Int_t actualNumRecPhotons = 0;

	// begin the event loop

	for (Int_t jentry = 0; jentry < numberOfEntries; jentry++){
		
		if (jentry % 5000 == 0){cout << "processing event no: " << jentry << endl;}

		// read the event
		treeReader->ReadEntry(jentry);

		// get number of rec and gen photons in this event
		Int_t rec_photon_size = branchRecPhoton->GetEntriesFast();
		Int_t gen_particle_size = branchGenPhoton->GetEntriesFast();

		Int_t numGenPhotons = 0;
		// loop through the generated photons
		for (Int_t g = 0; g < gen_particle_size; g++){
			
			// create generated particle object
			GenParticle *gen_particle = (GenParticle*) branchGenPhoton->At(g);
			
			// make sure it's a photon
			if (gen_particle->PID != 22){continue;}

			// make cuts on fiducial volume and on PT to eliminate photons due to pions and stuff
			// ideally you'd make a cut based on the parents, but until that is figured out this will have to do
			if (gen_particle->PT < 20 || abs(gen_particle->Eta) > 2.4){continue;} 

			// increment the number of generated photons in this event
			numGenPhotons++;
		}
		// after looping through the generated photons and counting them, 
		// fill gen_photon_distribution histogram with this value
		gen_photon_distribution_hist->Fill(numGenPhotons);

		Int_t numRecPhotons = 0;
		// loop through the generated photons
		for (Int_t r = 0; r < rec_photon_size; r++){
 			
			// create reconstructed photon object
 			Photon *rec_photon = (Photon*)branchRecPhoton->At(r);
 			
			// make cuts on fiducial volume and on PT to eliminate photons due to pions and stuff
			// ideally you'd make a cut based on the parents, but until that is figured out this will have to do
			if (rec_photon->PT < 20 || abs(rec_photon->Eta) > 2.4){continue;} 			

			// increment the number of generated photons in this event
			numRecPhotons++;
		}

		// after looping through the generated photons and counting them, 
		// fill rec_photon_distribution histogram with this value
		rec_photon_distribution_hist->Fill(numRecPhotons);

	} // end of event loop

	// draw the histograms on a canvas
	TCanvas *photonCountCanvas = new TCanvas("photonCountCanvas","Distribution of Photon Count Per Event",10,10,675,520);

	photonCountCanvas->Divide(2,1);

	photonCountCanvas->cd(1);
	gen_photon_distribution_hist->Draw();

	photonCountCanvas->cd(2);
	rec_photon_distribution_hist->Draw();

} // end of photonCountDistribution function

void massOfParents(TLorentzVector parent1, TLorentzVector parent2){

}
