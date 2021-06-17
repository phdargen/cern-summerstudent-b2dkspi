/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <math.h>

#include "TMath.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include <TChain.h>
#include "TStopwatch.h"

//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace std;
using namespace TMVA;

void TMVAClassificationApplication(TString decay = "B2DKspi", TString dataType = "Data", TString myMethod = "BDTG", TString trainedOn = "MC" ) 
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------
   TChain* theTree = new TChain("DecayTree");

   TString inDir = "/eos/lhcb/user/p/phdargen/summerStudents21/";
   TString outFileName = inDir+"BDT/";
    
   if(dataType == "Data"){ 	  
    theTree->Add(inDir+"Stripped/Data_"+decay+"_LL_12.root");
    theTree->Add(inDir+"Stripped/Data_"+decay+"_DD_12.root");       
    outFileName += decay+"_data.root";
   }

   else if(dataType == "MC"){ 	  
      theTree->Add(inDir+"Stripped/MC_"+decay+"_DD_12.root");
      theTree->Add(inDir+"Stripped/MC_"+decay+"_LL_12.root");
      outFileName += decay+"_mc.root";
   }
   else {
   	cout << "Unknown options, I'll crash now." << endl;
   	throw "ERROR";
   }

   // Ouput tree
   TFile *hFile = new TFile(outFileName,"RECREATE");
   TTree* tree = theTree->CloneTree(0);

   // This loads the library
   TMVA::Tools::Instance();

   // --- Create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader    
   Float_t r_PV_CHI2NDOF;

   // Use same names and order as in TMVAClassification.cpp !
   reader->AddVariable( "PV_CHI2NDOF", &r_PV_CHI2NDOF );

   // --- Book the MVA methods
   TString prefix = "weights/TMVAClassification_" + decay + "_" + trainedOn + "_";

   // Options used in TMVAClassification.cpp 
   std::vector<TString> weightFiles;
   weightFiles.push_back("all_LL");
   weightFiles.push_back("all_DD");
   //weightFiles.push_back("run1_LL");
   //weightFiles.push_back("run2_LL");
   //weightFiles.push_back("run1_DD");
   //weightFiles.push_back("run2_DD");

   for(int i= 0 ; i < weightFiles.size(); i++) 
        reader->BookMVA( myMethod + weightFiles[i], prefix + weightFiles[i] + "_" + myMethod + ".weights.xml" ); 

   Double_t PV_CHI2NDOF;
    
   theTree->SetBranchAddress( "PV_CHI2NDOF", &PV_CHI2NDOF );
   
   Int_t run, KsCat; 
   theTree->SetBranchAddress( "run", &run );
   theTree->SetBranchAddress( "KsCat", &KsCat );
   
   //output file--------------------------------------------------------------------------------------------------------------------------
   Float_t BDTG_response;
   double BDTG;
   tree->Branch("BDTG_response",&BDTG_response, "BDTG_response/F");
   tree->Branch("BDTG",&BDTG, "BDTG/D");

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;

   TStopwatch sw;
   sw.Start();
   int N = theTree->GetEntries();
   for (Long64_t ievt=0; ievt< N ;ievt++) {

      if (ievt%5000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

        theTree->GetEntry(ievt);
         
        r_PV_CHI2NDOF = float(PV_CHI2NDOF);

        TString methodName = myMethod + "run";
        methodName += run;
        if(KsCat == 0) methodName += "_LL" ;
        else methodName += "_DD";
        methodName += "_all" ;
       
        BDTG_response=reader->EvaluateMVA(methodName);
        BDTG = double(BDTG_response);
              
        tree->Fill();    
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   tree->Write();
   delete reader;
   hFile->Close();
    
   cout << "Wrote to file: " << outFileName << endl;
   cout << "==> TMVAClassificationApplication is done!" << endl << endl;
} 

void applyToAll(TString myMethod = "BDTG", TString trainedOn = "MC" ){
    TMVAClassificationApplication("B2DKspi", "Data", myMethod, trainedOn );
    TMVAClassificationApplication("B2DKspi", "MC", myMethod, trainedOn );   
}
