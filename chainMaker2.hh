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
#include "/home/habdalla/Delphes-3.0.8/classes/DelphesClasses.h"
#include "/home/habdalla/Delphes-3.0.8/modules/Delphes.h"
#include "/home/habdalla/Delphes-3.0.8/external/ExRootAnalysis/ExRootTreeReader.h"

void chainMaker2(TChain *chain) {

  TString catalogDir = "/home/cmsprod/catalog/t2mit";
  TString catalog = "/delphes/309";

  //change this for whatever you want
  TString run = "s12-t14-h125gg-vh-ul1tdr-v1-aod";

  cout << catalogDir << catalog << "/" << run << "/Files" << endl;
  ifstream ifs;
  ifs.open(catalogDir+catalog+TString("/")+run+TString("/Files"));
  assert(ifs.is_open());
  string line;
  while (getline(ifs, line)) {
    if (line[0] == ' ') continue;
    else {
      stringstream ss(line);
      TString fileset;
      TString filename;
      Int_t thing3, thing4, thing5, thing6, thing7, thing8;
      ss >> fileset >> filename >> thing3 >> thing4 >> thing5 >> thing6 >> thing7 >> thing8;
      //cout << "root://xrootd.cmsaf.mit.edu//store/user/paus" << catalog << "/" << run << "/" << filename <<  endl;

      //this option if you haven't downloaded it to hadoop
      chain->Add(TString("root://xrootd.cmsaf.mit.edu//store/user/paus")+catalog+TString("/")+run+TString("/")+filename);

      //this option if you've downloaded it to hadoop
      //chain->Add(TString("/mnt/hadoop/cms/store/user/paus")+catalog+TString("/")+run+TString("/")+filename);
    }
  }
  ifs.close();
}

#endif
