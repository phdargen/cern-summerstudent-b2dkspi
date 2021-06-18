#include <cstdlib>
#include <iostream>
#include <map>
#include <math.h>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVAGui.C"
#include "../TMVA-v4.2.0/test/variables.C"
#include "../TMVA-v4.2.0/test/efficiencies.C"
#include "../TMVA-v4.2.0/test/mvas.C"
#include "../TMVA-v4.2.0/test/correlations.C"
#include "../TMVA-v4.2.0/test/mvaeffs.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

void TMVAClassification( TString myMethodList = "BDTG", TString trainOn = "MC", TString decay = "B2DKspi", TString run = "all", TString Ks = "all" )
{

   // Input files
   TChain* background = new TChain("DecayTree");
   TString inDir = "/eos/lhcb/user/p/phdargen/summerStudents21/Stripped/";

   // Change to your selected files
   if(Ks == "LL" || Ks == "all")background->Add(inDir+"Data_"+decay+"_LL_12.root");
   if(Ks == "DD" || Ks == "all")background->Add(inDir+"Data_"+decay+"_DD_12.root");

   TChain* signal = new TChain("DecayTree");
   if(trainOn == "MC"){
       if(Ks == "LL" || Ks == "all")signal->Add(inDir+"MC_"+decay+"_LL_12.root");
       if(Ks == "DD" || Ks == "all")signal->Add(inDir+"MC_"+decay+"_DD_12.root");
   }
   else signal->Add("");

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.

   //TString outDir = "plots";
   TString outDir = "figs/TMVA/";
   outDir +=  myMethodList+ "_" + decay + "_" + trainOn + "_" + run + "_" + Ks;
 
   TString outfileName = "TMVA_" + decay + "_" +  myMethodList + "_" + trainOn + "_" + run + "_" + Ks + ".root";
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. 
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification_" + decay + "_" + trainOn + "_" + run + "_" + Ks, outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   signal->SetBranchStatus("*",0);  // disable all branches
   signal->SetBranchStatus("*PT*",1); 
   signal->SetBranchStatus("*CHI2*",1); 
   signal->SetBranchStatus("weight",1);
   signal->SetBranchStatus("run",1);
   signal->SetBranchStatus("B_DTF_MM",1);
   signal->SetBranchStatus("B_BKGCAT",1);

   background->SetBranchStatus("*",0);  // disable all branches
   background->SetBranchStatus("*PT*",1); 
   background->SetBranchStatus("*CHI2*",1); 
   background->SetBranchStatus("weight",1);
   background->SetBranchStatus("run",1);
   background->SetBranchStatus("B_DTF_MM",1);
   background->SetBranchStatus("B_BKGCAT",1);

   // Define the input variables that shall be used for the MVA training
   factory->AddVariable( "B_PT", "B_PT", "MeV", 'F' );
   factory->AddVariable( "PV_CHI2NDOF", "#chi^{2}_{DTF}/ndf", "", 'F' );

   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts;
   //if(run != "all")mycuts += "run == " + run.ReplaceAll("run","");
   if(trainOn == "MC") mycuts += "B_DTF_MM > 5000 && B_DTF_MM < 6000 && B_BKGCAT < 30";

   TCut mycutb = "B_DTF_MM > 5500";
   //if(run != "all")mycutb += "run == " + run;
   
   factory->AddSignalTree    ( signal,     signalWeight     );
   factory->AddBackgroundTree( background, backgroundWeight );
   factory->PrepareTrainingAndTestTree( mycuts,mycutb,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   //factory->SetSignalWeightExpression("weight");

   // ---- Book MVA methods

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=500:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:nCuts=40:MaxDepth=3:NegWeightTreatment=Pray" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=500:MinNodeSize=3.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=80" );

   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

   if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   variables(outfileName,"InputVariables_Id", "TMVA Input Variables",kFALSE, kTRUE, outDir);
   correlations( outfileName,  kFALSE, kFALSE, kTRUE ,outDir);
   efficiencies( outfileName,  2, kTRUE ,outDir);
   mvas( outfileName, CompareType,  kTRUE , outDir, true);

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );
 
}

void trainAll( TString myMethodList = "BDTG", TString trainOn = "MC") {
	gROOT->SetBatch(true);
	TMVAClassification( myMethodList, trainOn , "B2DKspi", "all",  "LL" );
 	TMVAClassification( myMethodList, trainOn , "B2DKspi", "all",  "DD" );
}
