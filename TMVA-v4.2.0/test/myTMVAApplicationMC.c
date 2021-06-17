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

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;

void myTMVAApplication( TString myMethodList = "BDT" ) 
{   
   string sample = "2";
   //string sample = "2";
   string inFileName;
   string outFileName;

   inFileName="/auto/data/dargent/Bu2JpsiKpipi/MC/MC_reweighted_forBDT_PIDcorr_cat10.root";
   outFileName="/auto/data/dargent/Bu2JpsiKpipi/MC/MC_bdt_PIDcorr_cat10.root";

#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------

   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 1;
   Use["CutsD"]           = 1;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 1;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 1;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 1;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 1; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 1; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 1;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 0; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 1;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod 
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   /**Float_t var1, var2;
   Float_t var3, var4;
   reader->AddVariable( "myvar1 := var1+var2", &var1 );
   reader->AddVariable( "myvar2 := var1-var2", &var2 );
   reader->AddVariable( "var3",                &var3 );
   reader->AddVariable( "var4",                &var4 );
	*/

   Float_t var[15];
   Float_t var2[15];

   reader->AddVariable("Bplus_PT:=Bplus_PT/1000.",&var[0]);
   //reader->AddVariable("log_Bplus_TAU := log(Bplus_TAU)",&var[1]);
   reader->AddVariable("Bplus_IPCHI2_OWNPV",&var[2]);
   reader->AddVariable("logBplus_FDCHI2 := log(Bplus_FDCHI2_OWNPV)",&var[3]);
   reader->AddVariable("log_Bplus_DIRA := log(1-Bplus_DIRA_OWNPV)",&var[4]);
   reader->AddVariable("Bplus_ENDVERTEX_CHI2",&var[5]);

   //reader->AddVariable("Kplus_PT:=Kplus_PT/1000.",&var[6]);
   //reader->AddVariable("piplus_PT:=piplus_PT/1000.",&var[7]);
   //reader->AddVariable("piminus_PT:=piminus_PT/1000.",&var[8]);
   //reader->AddVariable("muplus_PT:=muplus_PT/1000.",&var[9]);
   //reader->AddVariable("muminus_PT:=muminus_PT/1000.",&var[10]);

   //reader->AddVariable( "log_Kplus_IPCHI2 := log(Kplus_IPCHI2_OWNPV)",&var[11]);
   //reader->AddVariable("log_piplus_IPCHI2:=log(piplus_IPCHI2_OWNPV)",&var[12]);
   //reader->AddVariable("log_piminus_IPCHI2:=log(piminus_IPCHI2_OWNPV)",&var[13]);
   //reader->AddVariable("log_muplus_IPCHI2:=log(muplus_IPCHI2_OWNPV)",&var[14]);
   //reader->AddVariable("log_muminus_IPCHI2:=log(muminus_IPCHI2_OWNPV)",&var[15]);
   reader->AddVariable( "log_min_IPCHI2:=log(min(min(min(min(muplus_IPCHI2_OWNPV,muminus_IPCHI2_OWNPV),Kplus_IPCHI2_OWNPV),piminus_IPCHI2_OWNPV),piplus_IPCHI2_OWNPV))", &var2[0] );
   reader->AddVariable( "cos := cos( max(max(angK,angPip),angPim))", &var2[1] );
   reader->AddVariable( "Kplus_PIDK", &var2[2] );

   // Spectator variables declared in the training have to be added to the reader, too
   ///Float_t spec1,spec2;
   ///reader->AddSpectator( "spec1 := var1*2",   &spec1 );
   ///reader->AddSpectator( "spec2 := var1*3",   &spec2 );

   Float_t spec[15];
   reader->AddSpectator( "Bplus_MM",   &spec[0] );
   reader->AddSpectator( "mKpi",  &spec[0] );
   reader->AddSpectator( "mpipi",  &spec[1] );
   reader->AddSpectator( "mKpipi", &spec[2] );
   reader->AddSpectator( "mJpsipi",  &spec[3] );
   reader->AddSpectator( "mJpsipipi",  &spec[4] );
   reader->AddSpectator( "DTF_mKpi",  &spec[5] );
   reader->AddSpectator( "DTF_mpipi",  &spec[6] );
   reader->AddSpectator( "DTF_mKpipi",  &spec[7]);
   reader->AddSpectator( "DTF_mJpsipi",  &spec[8] );
   reader->AddSpectator( "DTF_mJpsipipi",  &spec[9] );

/**
   Float_t Category_cat1, Category_cat2, Category_cat3;
   if (Use["Category"]){
      // Add artificial spectators for distinguishing categories
      reader->AddSpectator( "Category_cat1 := var3<=0",             &Category_cat1 );
      reader->AddSpectator( "Category_cat2 := (var3>0)&&(var4<0)",  &Category_cat2 );
      reader->AddSpectator( "Category_cat3 := (var3>0)&&(var4>=0)", &Category_cat3 );
   }
*/
   // --- Book the MVA methods

   TString dir    = "weights/";
   TString prefix = ("PIDK_TMVA_s"+sample).c_str();

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile ); 
      }
   }
   
   // Book output histograms
   UInt_t nbin = 100;
   TH1F   *histLk(0), *histLkD(0), *histLkPCA(0), *histLkKDE(0), *histLkMIX(0), *histPD(0), *histPDD(0);
   TH1F   *histPDPCA(0), *histPDEFoam(0), *histPDEFoamErr(0), *histPDEFoamSig(0), *histKNN(0), *histHm(0);
   TH1F   *histFi(0), *histFiG(0), *histFiB(0), *histLD(0), *histNn(0),*histNnbfgs(0),*histNnbnn(0);
   TH1F   *histNnC(0), *histNnT(0), *histBdt(0), *histBdtG(0), *histBdtD(0), *histRf(0), *histSVMG(0);
   TH1F   *histSVMP(0), *histSVML(0), *histFDAMT(0), *histFDAGA(0), *histCat(0), *histPBdt(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   TFile *input(0);
   input = TFile::Open( inFileName.c_str() ); // check if file in local directory exists
   
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;
   
   // --- Event loop

   // Prepare the event tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select data sample" << std::endl;

   TTree* dataTree = (TTree*)input->Get("DecayTree");

   ///variables used in bdt training
   Double_t mvar[15];
   Double_t mvar2[15];
   
   dataTree->SetBranchAddress( "Bplus_PT", &mvar[0] );
   dataTree->SetBranchAddress( "Bplus_TAU", &mvar[1] );
   dataTree->SetBranchAddress( "Bplus_IPCHI2_OWNPV", &mvar[2] );
   dataTree->SetBranchAddress( "Bplus_FDCHI2_OWNPV", &mvar[3] );
   dataTree->SetBranchAddress( "Bplus_DIRA_OWNPV", &mvar[4] );
   dataTree->SetBranchAddress( "Bplus_ENDVERTEX_CHI2", &mvar[5] );

   dataTree->SetBranchAddress( "Kplus_PT", &mvar[6] );
   dataTree->SetBranchAddress( "piplus_PT", &mvar[7] );
   dataTree->SetBranchAddress( "piminus_PT", &mvar[8] );
   dataTree->SetBranchAddress( "muplus_PT", &mvar[9] );
   dataTree->SetBranchAddress( "muminus_PT", &mvar[10] );

   dataTree->SetBranchAddress( "Kplus_IPCHI2_OWNPV", &mvar[11] );
   dataTree->SetBranchAddress( "piplus_IPCHI2_OWNPV", &mvar[12] );
   dataTree->SetBranchAddress( "piminus_IPCHI2_OWNPV", &mvar[13] );
   dataTree->SetBranchAddress( "muplus_IPCHI2_OWNPV", &mvar[14] );
   dataTree->SetBranchAddress( "muminus_IPCHI2_OWNPV", &mvar[15] );

   dataTree->SetBranchAddress( "angK", &mvar2[0] );
   dataTree->SetBranchAddress( "angPip", &mvar2[1] );
   dataTree->SetBranchAddress( "angPim", &mvar2[2] );
   dataTree->SetBranchAddress( "Kplus_PIDK_corr", &mvar2[3] );

/*
  ///Needed branches for analysis
  Double_t Bplus_MM, K_1_1270_plus_MM, J_psi_1S_MM;
  Double_t Bplus_TAU;
  Double_t Bplus_PX,Bplus_PY,Bplus_PZ,Bplus_PE;
  Double_t J_psi_1S_PE, J_psi_1S_PX, J_psi_1S_PY, J_psi_1S_PZ,J_psi_1S_PT;
  Double_t Kplus_P,Kplus_PE, Kplus_PX, Kplus_PY, Kplus_PZ;
  Double_t piplus_P,piplus_PE, piplus_PX, piplus_PY, piplus_PZ;
  Double_t piminus_P,piminus_PE, piminus_PX, piminus_PY, piminus_PZ;
  Double_t muplus_PE, muplus_PX, muplus_PY, muplus_PZ;
  Double_t muminus_PE, muminus_PX, muminus_PY, muminus_PZ;
  Double_t K1[5];
  Double_t ev[4];

  dataTree->SetBranchAddress("nCandidate",&ev[0]) ;
  dataTree->SetBranchAddress("nTracks",&ev[1]) ;
  dataTree->SetBranchAddress("nPV",&ev[2]) ;
  dataTree->SetBranchAddress("eventNumber",&ev[3]) ;

  dataTree->SetBranchAddress("Bplus_MM",&Bplus_MM) ;
  dataTree->SetBranchAddress("J_psi_1S_MM",&J_psi_1S_MM) ;
  dataTree->SetBranchAddress("K_1_1270_plus_MM",&K_1_1270_plus_MM) ;
  
  dataTree->SetBranchAddress( "Bplus_PX", &Bplus_PX );
  dataTree->SetBranchAddress( "Bplus_PY", & Bplus_PY);
  dataTree->SetBranchAddress( "Bplus_PZ", & Bplus_PZ);
  dataTree->SetBranchAddress( "Bplus_PE", &Bplus_PE);

  dataTree->SetBranchAddress("J_psi_1S_PE",&J_psi_1S_PE) ;
  dataTree->SetBranchAddress("J_psi_1S_PX",&J_psi_1S_PX) ;
  dataTree->SetBranchAddress("J_psi_1S_PY",&J_psi_1S_PY) ;
  dataTree->SetBranchAddress("J_psi_1S_PZ",&J_psi_1S_PZ) ;
  dataTree->SetBranchAddress("J_psi_1S_PT",&J_psi_1S_PT) ;

  dataTree->SetBranchAddress("K_1_1270_plus_PT",&K1[0]) ;
  dataTree->SetBranchAddress("K_1_1270_plus_PE",&K1[1]) ;
  dataTree->SetBranchAddress("K_1_1270_plus_PX",&K1[2]) ; 
  dataTree->SetBranchAddress("K_1_1270_plus_PY",&K1[3]) ;
  dataTree->SetBranchAddress("K_1_1270_plus_PZ",&K1[4]) ;

  dataTree->SetBranchAddress("Kplus_P",&Kplus_P) ;
  dataTree->SetBranchAddress("Kplus_PE",&Kplus_PE) ;
  dataTree->SetBranchAddress("Kplus_PX",&Kplus_PX) ;
  dataTree->SetBranchAddress("Kplus_PY",&Kplus_PY) ;
  dataTree->SetBranchAddress("Kplus_PZ",&Kplus_PZ) ;
  
  dataTree->SetBranchAddress("piplus_P",&piplus_P) ;
  dataTree->SetBranchAddress("piplus_PE",&piplus_PE) ;
  dataTree->SetBranchAddress("piplus_PX",&piplus_PX) ;
  dataTree->SetBranchAddress("piplus_PY",&piplus_PY) ;
  dataTree->SetBranchAddress("piplus_PZ",&piplus_PZ) ;
  
  dataTree->SetBranchAddress("piminus_P",&piminus_P) ;
  dataTree->SetBranchAddress("piminus_PE",&piminus_PE) ;
  dataTree->SetBranchAddress("piminus_PX",&piminus_PX) ;
  dataTree->SetBranchAddress("piminus_PY",&piminus_PY) ;
  dataTree->SetBranchAddress("piminus_PZ",&piminus_PZ) ;
    
  dataTree->SetBranchAddress("muplus_PE",&muplus_PE) ;
  dataTree->SetBranchAddress("muplus_PX",&muplus_PX) ;
  dataTree->SetBranchAddress("muplus_PY",&muplus_PY) ;
  dataTree->SetBranchAddress("muplus_PZ",&muplus_PZ) ;
  
  dataTree->SetBranchAddress("muminus_PE",&muminus_PE) ;
  dataTree->SetBranchAddress("muminus_PX",&muminus_PX) ;
  dataTree->SetBranchAddress("muminus_PY",&muminus_PY) ;
  dataTree->SetBranchAddress("muminus_PZ",&muminus_PZ) ;
*/

  ///output file
  TFile *target = new TFile(outFileName.c_str(),"RECREATE" );
  TTree* new_tree = dataTree->CloneTree();  
  Double_t BDT_response;
  TBranch* Bra_bdt =new_tree->Branch(("BDT_response"+sample).c_str(), &BDT_response , ("BDT_response"+sample+"/D").c_str());

//TTree *tree = new TTree("tree","tree");
//   
// 
//   tree->Branch("BDT_response", &BDT_response );
//   tree->Branch("Bplus_MM",&Bplus_MM) ;
//   tree->Branch("J_psi_1S_MM",&J_psi_1S_MM) ;
//   tree->Branch("K_1_1270_plus_MM",&K_1_1270_plus_MM) ;
//   
//   tree->Branch("Bplus_TAU",&mvar[1]) ;
//   tree->Branch( "Bplus_PX", &Bplus_PX );
//   tree->Branch( "Bplus_PY", & Bplus_PY);
//   tree->Branch( "Bplus_PZ", & Bplus_PZ);
//   tree->Branch( "Bplus_PE", &Bplus_PE);
//   tree->Branch( "Bplus_ENDVERTEX_CHI2", &mvar[5] );
// 
//   tree->Branch("J_psi_1S_PE",&J_psi_1S_PE) ;
//   tree->Branch("J_psi_1S_PX",&J_psi_1S_PX) ;
//   tree->Branch("J_psi_1S_PY",&J_psi_1S_PY) ;
//   tree->Branch("J_psi_1S_PZ",&J_psi_1S_PZ) ;
//   tree->Branch("J_psi_1S_PT",&J_psi_1S_PT) ;
// 
//   tree->Branch("K_1_1270_plus_PT",&K1[0]) ;
//   tree->Branch("K_1_1270_plus_PE",&K1[1]) ;
//   tree->Branch("K_1_1270_plus_PX",&K1[2]) ;
//   tree->Branch("K_1_1270_plus_PY",&K1[3]) ;
//   tree->Branch("K_1_1270_plus_PZ",&K1[4]) ;
// 
//   tree->Branch("Kplus_P",&Kplus_P) ;
//   tree->Branch("Kplus_PE",&Kplus_PE) ;
//   tree->Branch("Kplus_PX",&Kplus_PX) ;
//   tree->Branch("Kplus_PY",&Kplus_PY) ;
//   tree->Branch("Kplus_PZ",&Kplus_PZ) ;
//   
//   tree->Branch("piplus_P",&piplus_P) ;
//   tree->Branch("piplus_PE",&piplus_PE) ;
//   tree->Branch("piplus_PX",&piplus_PX) ;
//   tree->Branch("piplus_PY",&piplus_PY) ;
//   tree->Branch("piplus_PZ",&piplus_PZ) ;
//   
//   tree->Branch("piminus_P",&piminus_P) ;
//   tree->Branch("piminus_PE",&piminus_PE) ;
//   tree->Branch("piminus_PX",&piminus_PX) ;
//   tree->Branch("piminus_PY",&piminus_PY) ;
//   tree->Branch("piminus_PZ",&piminus_PZ) ;
//     
//   tree->Branch("muplus_PE",&muplus_PE) ;
//   tree->Branch("muplus_PX",&muplus_PX) ;
//   tree->Branch("muplus_PY",&muplus_PY) ;
//   tree->Branch("muplus_PZ",&muplus_PZ) ;
//   
//   tree->Branch("muminus_PE",&muminus_PE) ;
//   tree->Branch("muminus_PX",&muminus_PX) ;
//   tree->Branch("muminus_PY",&muminus_PY) ;
//   tree->Branch("muminus_PZ",&muminus_PZ) ;
// 
//   tree->Branch("nCandidate",&ev[0]) ;
//   tree->Branch("nTracks",&ev[1]) ;
//   tree->Branch("nPV",&ev[2]) ;
//   tree->Branch("eventNumber",&ev[3]) ;

   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << dataTree->GetEntries() << " events" << std::endl;
   dataTree->SetBranchStatus("*",0);
   dataTree->SetBranchStatus("*CHI*",1);
   dataTree->SetBranchStatus("*PT",1);
   dataTree->SetBranchStatus("ang*",1);
   dataTree->SetBranchStatus("*TAU",1);
   dataTree->SetBranchStatus("*DIRA*",1);
   dataTree->SetBranchStatus("*PID*",1);
   dataTree->SetBranchStatus( "Bplus_BKGCAT", 1 );

   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<dataTree->GetEntries();ievt++) {

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

       dataTree->GetEntry(ievt);

       var[0]=mvar[0]/1000.;
       var[1]=log(mvar[1]);
       var[2]=mvar[2];
       var[3]=log(mvar[3]);
       var[4]=log(1-mvar[4]);
       var[5]=mvar[5];

       var[6]=mvar[6]/1000.;
       var[7]=mvar[7]/1000.;
       var[8]=mvar[8]/1000.;
       var[9]=mvar[9]/1000.;
       var[10]=mvar[10]/1000.;
       
       var[11]=log(mvar[11]);
       var[12]=log(mvar[12]);
       var[13]=log(mvar[13]);
       var[14]=log(mvar[14]);
       var[15]=log(mvar[15]);

       var2[0]=log(min(min(min(min(mvar[11],mvar[12]),mvar[13]),mvar[14]),mvar[15]));
       var2[1]=cos(max(max(mvar2[0],mvar2[1]),mvar2[2]));
       var2[2]=mvar2[3];

       BDT_response=reader->EvaluateMVA("BDT method");

       Bra_bdt->Fill();

      ///var1 = userVar1 + userVar2;
      ///var2 = userVar1 - userVar2;

      // --- Return the MVA outputs and fill into histograms

	
      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
      if (Use["BDT"          ])   histBdt    ->Fill( BDT_response );
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );         
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }         

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
         probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
         rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
      }
	
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   //tree->AddFriend(fname,"merged"); 
   //tree->Write();

   // Get efficiency for cuts classifier
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer  
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {      
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( 0.7, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: " 
                      << cutsMin[ivar] 
                      << " < \"" 
                      << mcuts->GetInputVar(ivar)
                      << "\" <= " 
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }

   // --- Write histograms

   ///TFile *target  = new TFile( "myTMVApp.root","RECREATE" );
   if (Use["Likelihood"   ])   histLk     ->Write();
   if (Use["LikelihoodD"  ])   histLkD    ->Write();
   if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
   if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
   if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
   if (Use["PDERS"        ])   histPD     ->Write();
   if (Use["PDERSD"       ])   histPDD    ->Write();
   if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
   if (Use["KNN"          ])   histKNN    ->Write();
   if (Use["HMatrix"      ])   histHm     ->Write();
   if (Use["Fisher"       ])   histFi     ->Write();
   if (Use["FisherG"      ])   histFiG    ->Write();
   if (Use["BoostedFisher"])   histFiB    ->Write();
   if (Use["LD"           ])   histLD     ->Write();
   if (Use["MLP"          ])   histNn     ->Write();
   if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
   if (Use["MLPBNN"       ])   histNnbnn  ->Write();
   if (Use["CFMlpANN"     ])   histNnC    ->Write();
   if (Use["TMlpANN"      ])   histNnT    ->Write();
   if (Use["BDT"          ])   histBdt    ->Write();
   if (Use["BDTD"         ])   histBdtD   ->Write();
   if (Use["BDTG"         ])   histBdtG   ->Write(); 
   if (Use["RuleFit"      ])   histRf     ->Write();
   if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
   if (Use["SVM_Poly"     ])   histSVMP   ->Write();
   if (Use["SVM_Lin"      ])   histSVML   ->Write();
   if (Use["FDA_MT"       ])   histFDAMT  ->Write();
   if (Use["FDA_GA"       ])   histFDAGA  ->Write();
   if (Use["Category"     ])   histCat    ->Write();
   if (Use["Plugin"       ])   histPBdt   ->Write();

   // Write also error and significance histos
   if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }

   // Write also probability hists
   if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }

   //target->WriteTObject(tree):
   new_tree->Write();
   std::cout << "--- Created root file: \"myTMVApp.root\" containing the MVA output histograms" << std::endl;
   delete reader;
   target->Close();
   input->Close(); 
   //delete tree;
   //delete dataTree;
    
   std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;
} 
