//#if !defined(__CINT__) || defined(__MAKECINT__)
#ifndef CHAINMAKER_HH
#define CHAINMAKER_HH

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <TString.h>
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
#include "modules/Delphes.h"

void chainMaker(TChain *chain) {

  //change this for whatever you want
  //TString run = "s12-zeem20-test-v9";

  ifstream ifs;
  ifs.open("/home/habdalla/Desktop/Delphes-3.0.8/root_files.txt");
  assert(ifs.is_open());
  string line;
  while (getline(ifs, line)) {
    if (line[0] == ' ') continue;
    else {
      stringstream ss(line);

      TString thing1;
      ss >> thing1;

      //this option if you haven't downloaded it to hadoop
	  //chain->Add(TString("/mnt/hadoop/dkralph/powheg/production/HH/dihiggs-bbgg-14tev/")+thing1);
      chain->Add(TString("/scratch3/habdalla/hi/post-hepmc-files/")+thing1);

      //this option if you've downloaded it to hadoop
      //chain->Add(TString("/mnt/hadoop/cms/store/user/paus")+catalog+TString("/")+run+TString("/")+filename);
    }
  }
  ifs.close();
}

#endif
