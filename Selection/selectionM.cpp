//Select events and fill them into a new tree
//Philippe d'Argent
#include <TChain.h>
#include <THStack.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TNtuple.h>

#include <iostream>
using namespace std;

int main() {

	 	//P: Load file		
    	TChain* tree=new TChain("DecayTree","RECREATE");
		tree->Add("/eos/lhcb/user/p/phdargen/summerStudents21/Stripped/Data_B2DKspi_DD_11.root");
		//tree->Add("/home/maria/Work/cern-summerstudent-b2dkspi/Maryia/Data/Data_B2DKspi_DD_12.root");
		//tree->Add("/home/maria/Work/cern-summerstudent-b2dkspi/Maryia/Data/Data_B2DKspi_DD_15.root");
		//tree->Add("/home/maria/Work/cern-summerstudent-b2dkspi/Maryia/Data/Data_B2DKspi_DD_16.root");
		//tree->Add("/home/maria/Work/cern-summerstudent-b2dkspi/Maryia/Data/Data_B2DKspi_DD_17.root");
		//tree->Add("/home/maria/Work/cern-summerstudent-b2dkspi/Maryia/Data/Data_B2DKspi_DD_18.root");
		//tree->Add("/home/maria/Work/cern-summerstudent-b2dkspi/Maryia/Data/MC_B2DKspi_DD_12_BdDstKspi.root");
		//tree->Add("/home/maria/Work/cern-summerstudent-b2dkspi/Maryia/Data/MC_B2DKspi_DD_12_BsDstKsK.root");
		//tree->Add("/home/maria/Work/cern-summerstudent-b2dkspi/Maryia/Data/MC_B2DKspi_DD_12_BsDstKspi.root");
		//tree->Add("/home/maria/Work/cern-summerstudent-b2dkspi/Maryia/Data/DMC_B2DKspi_DD_12.root");


	  	//P: Needed branches
  		Double_t B_MM,B_DTF_MM,B_PT,B_IPCHI2_OWNPV,B_FDCHI2_OWNPV,B_TAU,D_ENDVERTEX_Z, B_ENDVERTEX_Z;
  		Double_t D_FDCHI2_ORIVX,D_DIRA_OWNPV,Ks_FDCHI2_ORIVX, D_MM, Ks_MM, ProbNNpi, ProbNNk;
  		Double_t Ks_PT, Ks_DIRA_OWNPV,B_ENDVERTEX_CHI2, TRACK_GhostProb;
  		Int_t B_BKGCAT, B_ENDVERTEX_NDOF, KsCat;
  		bool hasRich; 
  		
  		tree->SetBranchAddress("B_MM",&B_MM) ;
		tree->SetBranchAddress("D_MM",&D_MM) ;
		tree->SetBranchAddress("Ks_MM",&Ks_MM) ;
		tree->SetBranchAddress("B_DTF_MM",&B_DTF_MM) ;
    	tree->SetBranchAddress("B_PT",&B_PT) ;
    	tree->SetBranchAddress("B_IPCHI2_OWNPV",&B_IPCHI2_OWNPV) ;
	   	tree->SetBranchAddress("B_FDCHI2_OWNPV",&B_FDCHI2_OWNPV) ;
    	tree->SetBranchAddress("B_TAU",&B_TAU) ;
    	tree->SetBranchAddress("D_ENDVERTEX_Z",&D_ENDVERTEX_Z) ;
    	tree->SetBranchAddress("B_ENDVERTEX_Z",&B_ENDVERTEX_Z) ;
    	tree->SetBranchAddress("D_FDCHI2_ORIVX",&D_FDCHI2_ORIVX) ;
 	  	tree->SetBranchAddress("D_DIRA_OWNPV",&D_DIRA_OWNPV) ;
    	tree->SetBranchAddress("Ks_FDCHI2_ORIVX",&Ks_FDCHI2_ORIVX) ;
    	tree->SetBranchAddress("Ks_PT",&Ks_PT) ;
     	tree->SetBranchAddress("Ks_DIRA_OWNPV",&Ks_DIRA_OWNPV) ;
 		tree->SetBranchAddress("B_ENDVERTEX_CHI2",&B_ENDVERTEX_CHI2) ;
 		tree->SetBranchAddress("B_ENDVERTEX_NDOF",&B_ENDVERTEX_NDOF) ;
 		tree->SetBranchAddress("KsCat",&KsCat) ;
 		tree->SetBranchAddress("B_BKGCAT",&B_BKGCAT) ;

		//M: Open root file (https://root.cern/manual/storing_root_objects/)
		TFile* input = TFile::Open("/eos/lhcb/user/p/phdargen/summerStudents21/Stripped/Data_B2DKspi_DD_11.root");

		//P: Create output file 
		TFile* output = new TFile("output.root","RECREATE");
  		TTree* summary_tree = tree->CloneTree(0);

    	//P: Add new branches
    	summary_tree->Branch("B_DTF_MM",&B_DTF_MM);
    	summary_tree->Branch("D_MM",&D_MM);
    	summary_tree->Branch("Ks_MM",&Ks_MM);

    
		//P: Define some histograms
		//M: To draw histograms on the same plot (example https://root.cern.ch/doc/master/classTHStack.html)
		THStack *hisB_DTF_MM = new THStack("hisB_DTF_MM","B_DTF_MM 2"); 
		THStack *hisD_MM = new THStack("hisD_MM","D_MM 2");
		THStack *hisKs_MM = new THStack("hisKs_MM","Ks_MM 2");
		
		//M: Histograms after the selection
    	TH1F *hB_DTF_MM = new TH1F("hB_DTF_MM","B_DTF_MM after",100,4700,6100);
   		TH1F *hD_MM = new TH1F("hD_MM","D_MM after",100,1760,1980);
    	TH1F *hKs_MM = new TH1F("hKs_MM","Ks_MM after",100,460,530);
    			
		//M: Histograms before the selection 

		TH1F *hbB_DTF_MM = (TH1F*)input -> Get("B_DTF_MM");
		TH1F *hbD_MM = (TH1F*)input -> Get("D_MM");
		TH1F *hbKs_MM = (TH1F*)input -> Get("Ks_MM");

		//hbB_DTF_MM->SetFillColor(kRed);
		//hbD_MM->SetFillColor(kRed);
		//hbKs_MM->SetFillColor(kRed);
		
		//M: Fill THStack with hist before the selection
		hisB_DTF_MM -> Add(hbB_DTF_MM);
		hisD_MM -> Add(hbB_DTF_MM);
		hisKs_MM -> Add(hbB_DTF_MM);
    	    
    	//P: Loop over tree
    	int nEvents = tree->GetEntries();
  		for ( Int_t j = 0 ; j < nEvents ; j++ ) {
    		tree->GetEntry(j) ;
    		if (0ul == (j % 2000ul)) cout << "Read event " << j << "/" << nEvents << endl;
		
		    //P: Add your cuts
		    //M: Ks category
		    if (KsCat != 1) continue;
		    
		    //M: Kinematic		    
            if (B_PT<=2000) continue;
            if (B_IPCHI2_OWNPV>=20) continue;
			if (B_FDCHI2_OWNPV<=200) continue;
		    if (B_TAU<=0.0001) continue;
		   	if (D_ENDVERTEX_Z - B_ENDVERTEX_Z <= 0) continue;
		   	if (D_FDCHI2_ORIVX <= 0) continue;
		   	if (D_DIRA_OWNPV <= 0) continue;
		   	if (Ks_FDCHI2_ORIVX <= 0) continue; 
		   	if (Ks_DIRA_OWNPV <= 0) continue; 
		   	if (Ks_PT <= 200) continue;
		   	
		    //M: Quality
		   	if (TRACK_GhostProb >= 0.5) continue;
		    if ((B_ENDVERTEX_CHI2/B_ENDVERTEX_NDOF) >= 10) continue;
		    //if (B_BKGCAT >= 30) continue; // Only for MC!
		    		    
		    //M: PID
		    //if (ProbNNpi <= 0.2) continue;
		    //if (ProbNNk <= 0.3) continue;
		    //if (hasRich != true) continue;
		    
		    		       
		    //P: Fill histograms
            hB_DTF_MM -> Fill(B_DTF_MM);
            hB_DTF_MM -> SetFillColor(kBlue-7);
            
            hD_MM -> Fill(D_MM);
            hD_MM -> SetFillColor(kBlue-7);
            
            hKs_MM -> Fill(Ks_MM);
            hKs_MM -> SetFillColor(kBlue-7);

        
            //P: Fill output tree
            summary_tree->Fill();
		}

		cout << "before selection: " << tree->GetEntries() << endl;
		cout << "after selection: " << summary_tree->GetEntries() << endl;

		output->Write();
		
		//M: Fill THStack with hist after the selection
		hisB_DTF_MM -> Add(hB_DTF_MM);
		hisD_MM -> Add(hD_MM);
		hisKs_MM -> Add(hKs_MM);
			
		//M: Canvas
		TCanvas *cB_DTF_MM = new TCanvas("cB_DTF_MM","cB_DTF_MM Distribution");
		hisB_DTF_MM -> Draw();
		hisB_DTF_MM-> SaveAs("/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/Selection/hisB_DTF_MM.root");
		
		TCanvas *cD_MM = new TCanvas("cD_MM","cD_MM Distribution");
		hisD_MM -> Draw();
		hisD_MM-> SaveAs("/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/Selection/hisD_MM.root");
		
		TCanvas *cKs_MM = new TCanvas("cKs_MM","cKs_MM Distribution");
		hisKs_MM -> Draw();
		hisKs_MM-> SaveAs("/afs/cern.ch/work/m/mbuhayeu/public/cern-summerstudent-b2dkspi-master/Selection/hisKs_MM.root");

		output -> Close();
		input -> Close();

return 0;
}
