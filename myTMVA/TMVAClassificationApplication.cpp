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

   TString inDir = "/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/";
   TString outFileName = inDir+"BDTG/";
    
   // Change to your selected files 
   if(dataType == "Data"){ 	  
    theTree->Add(inDir+"Selection/Stripped/Data_"+decay+"_LL_11_s.root");
    theTree->Add(inDir+"Selection/Stripped/Data_"+decay+"_LL_12_s.root");
    theTree->Add(inDir+"Selection/Stripped/Data_"+decay+"_LL_15_s.root");
    theTree->Add(inDir+"Selection/Stripped/Data_"+decay+"_LL_16_s.root");
    theTree->Add(inDir+"Selection/Stripped/Data_"+decay+"_LL_17_s.root");
    theTree->Add(inDir+"Selection/Stripped/Data_"+decay+"_LL_18_s.root");
    theTree->Add(inDir+"Selection/Stripped/Data_"+decay+"_DD_11_s.root"); 
    theTree->Add(inDir+"Selection/Stripped/Data_"+decay+"_DD_12_s.root"); 
    theTree->Add(inDir+"Selection/Stripped/Data_"+decay+"_DD_15_s.root"); 
    theTree->Add(inDir+"Selection/Stripped/Data_"+decay+"_DD_16_s.root"); 
    theTree->Add(inDir+"Selection/Stripped/Data_"+decay+"_DD_17_s.root"); 
    theTree->Add(inDir+"Selection/Stripped/Data_"+decay+"_DD_18_s.root");       
    outFileName += decay+"_data.root";
   }
   else if(dataType == "MC"){ 	  
      theTree->Add(inDir+"Selection/Stripped/MC_"+decay+"_DD_12_s.root");
      theTree->Add(inDir+"Selection/Stripped/MC_"+decay+"_LL_12_s.root");
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
   Float_t r_B_PT;
   Float_t r_PV_CHI2NDOF;
   Float_t r_weight;
   Float_t r_run;   
   Float_t r_B_DTF_MM;
   Float_t r_B_BKGCAT;   
   Float_t r_B_IPCHI2_OWNPV;
   Float_t r_B_DIRA_OWNPV;   
   Float_t r_B_FDCHI2_OWNPV;
   Float_t r_B_TAU;   
   Float_t r_D_ENDVERTEX_Z;
   Float_t r_B_ENDVERTEX_Z;   
   Float_t r_D_FDCHI2_ORIVX;
   Float_t r_D_DIRA_OWNPV;   
   Float_t r_Ks_FDCHI2_ORIVX;
   Float_t r_Ks_PT;
   Float_t r_Ks_DIRA_OWNPV;   
   Float_t r_B_ENDVERTEX_CHI2;
   Float_t r_pi_ProbNNpi;
   Float_t r_K_D_ProbNNk;
   Float_t r_pi_PT;
   
   // Use same names and order as in TMVAClassification.cpp !
   reader->AddVariable( "B_PT", &r_B_PT );
   reader->AddVariable( "PV_CHI2NDOF", &r_PV_CHI2NDOF );
   
   //reader->AddVariable("weight",&r_weight);
   //reader->AddVariable("run",&r_run);
   //reader->AddVariable("B_DTF_MM",&r_B_DTF_MM);
   //reader->AddVariable("B_BKGCAT",&r_B_BKGCAT);
   //reader->AddVariable("B_IPCHI2_OWNPV",&r_B_IPCHI2_OWNPV);// M:
   reader->AddVariable("log(1-B_DIRA_OWNPV)",&r_B_DIRA_OWNPV);// M:
   //reader->AddVariable("B_FDCHI2_OWNPV",&r_B_FDCHI2_OWNPV);// M:
   reader->AddVariable("B_TAU",&r_B_TAU);// M:
   //reader->AddVariable("D_ENDVERTEX_Z",&r_D_ENDVERTEX_Z);// M:
   reader->AddVariable("B_ENDVERTEX_Z",&r_B_ENDVERTEX_Z);// M:
   reader->AddVariable("log(D_FDCHI2_ORIVX)",&r_D_FDCHI2_ORIVX);// M:
   //reader->AddVariable("D_DIRA_OWNPV",&r_D_DIRA_OWNPV);// M:
   reader->AddVariable("log(Ks_FDCHI2_ORIVX)",&r_Ks_FDCHI2_ORIVX);// M:
   reader->AddVariable("Ks_PT",&r_Ks_PT);// M:
   //reader->AddVariable("Ks_DIRA_OWNPV",&r_Ks_DIRA_OWNPV);// M:
   reader->AddVariable("B_ENDVERTEX_CHI2",&r_B_ENDVERTEX_CHI2);// M:
   reader->AddVariable("pi_ProbNNpi",&r_pi_ProbNNpi);// M:
   reader->AddVariable("K_D_ProbNNk",&r_K_D_ProbNNk);// M:
   reader->AddVariable("pi_PT",&r_pi_PT);// M:

   // --- Book the MVA methods
   TString prefix = "myTMVA/weights/TMVAClassification_" + decay + "_" + trainedOn + "_";

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

   Double_t B_PT, PV_CHI2NDOF;
   theTree->SetBranchAddress( "B_PT", &B_PT );
   theTree->SetBranchAddress( "PV_CHI2NDOF", &PV_CHI2NDOF );
   
   Int_t run, KsCat; 
   theTree->SetBranchAddress( "run", &run );
   theTree->SetBranchAddress( "KsCat", &KsCat );
   
   Double_t B_DTF_MM, B_IPCHI2_OWNPV, B_DIRA_OWNPV, B_FDCHI2_OWNPV, B_TAU, D_ENDVERTEX_Z, B_ENDVERTEX_Z, D_FDCHI2_ORIVX, D_DIRA_OWNPV, Ks_FDCHI2_ORIVX, Ks_PT, Ks_DIRA_OWNPV, B_ENDVERTEX_CHI2, pi_ProbNNpi, K_D_ProbNNk, pi_PT;
   //theTree->SetBranchAddress("B_DTF_MM",&B_DTF_MM);
   //theTree->SetBranchAddress("B_IPCHI2_OWNPV",&B_IPCHI2_OWNPV);// M:
   theTree->SetBranchAddress("B_DIRA_OWNPV",&B_DIRA_OWNPV);// M:
   //theTree->SetBranchAddress("B_FDCHI2_OWNPV",&B_FDCHI2_OWNPV);// M:
   theTree->SetBranchAddress("B_TAU",&B_TAU);// M:
   //theTree->SetBranchAddress("D_ENDVERTEX_Z",&D_ENDVERTEX_Z);// M:
   theTree->SetBranchAddress("B_ENDVERTEX_Z",&B_ENDVERTEX_Z);// M:
   theTree->SetBranchAddress("D_FDCHI2_ORIVX",&D_FDCHI2_ORIVX);// M:
   //theTree->SetBranchAddress("D_DIRA_OWNPV",&D_DIRA_OWNPV);// M:
   theTree->SetBranchAddress("Ks_FDCHI2_ORIVX",&Ks_FDCHI2_ORIVX);// M:
   theTree->SetBranchAddress("Ks_PT",&Ks_PT);// M:
   //theTree->SetBranchAddress("Ks_DIRA_OWNPV",&Ks_DIRA_OWNPV);// M:
   theTree->SetBranchAddress("B_ENDVERTEX_CHI2",&B_ENDVERTEX_CHI2);// M:
   theTree->SetBranchAddress("pi_ProbNNpi",&pi_ProbNNpi);// M:
   theTree->SetBranchAddress("K_D_ProbNNk",&K_D_ProbNNk);// M:
   theTree->SetBranchAddress("pi_PT",&pi_PT);// M:   
   
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
         
        r_B_PT = float(B_PT);
        r_PV_CHI2NDOF = float(PV_CHI2NDOF);
        //r_B_DTF_MM = float(B_DTF_MM); 
   		r_B_IPCHI2_OWNPV = float(B_IPCHI2_OWNPV);
   		r_B_DIRA_OWNPV = float(B_DIRA_OWNPV);  
   //	r_B_FDCHI2_OWNPV = float(B_FDCHI2_OWNPV);
   		r_B_TAU = float(B_TAU);   
   		r_D_ENDVERTEX_Z = float(D_ENDVERTEX_Z);
   		r_B_ENDVERTEX_Z = float(B_ENDVERTEX_Z);   
   		r_D_FDCHI2_ORIVX = float(D_FDCHI2_ORIVX);
   //	r_D_DIRA_OWNPV = float(D_DIRA_OWNPV);  
   		r_Ks_PT = float(Ks_PT);
   //	r_Ks_DIRA_OWNPV = float(Ks_DIRA_OWNPV);   
   		r_B_ENDVERTEX_CHI2 = float(B_ENDVERTEX_CHI2);
   		r_pi_ProbNNpi = float(pi_ProbNNpi);
   		r_K_D_ProbNNk = float(K_D_ProbNNk);
   		r_pi_PT = float(pi_PT);
        

        // Might need to change this depending on your options in TMVAClassification.cpp and your workflow
        TString methodName = myMethod; 
        methodName += "all";
        if(KsCat == 0) methodName += "_LL" ;
        else methodName += "_DD";
       
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
