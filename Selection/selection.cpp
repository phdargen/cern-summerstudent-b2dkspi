//Select events and fill them into a new tree
//Philippe d'Argent
#include <TChain.h>
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

 	//Load file	
    TChain* tree=new TChain("DecayTree");
    tree->Add("/eos/lhcb/user/p/phdargen/summerStudents21/Stripped/Data_B2DKspi_LL_11.root");
    //tree->Add("/eos/lhcb/user/p/phdargen/summerStudents21/Stripped/Data_B2DKspi_DD_11.root");

  	//Needed branches (add some more)
  	Double_t B_MM,B_PT;
	tree->SetBranchAddress("B_MM",&B_MM) ;
    tree->SetBranchAddress("B_PT",&B_PT) ;

	//Create output file
	TFile* output = new TFile("output.root","RECREATE");
    TTree* summary_tree = tree->CloneTree(0);

    //Add new branches
    //...
    
    //Define some histograms
    //...
    
    //Loop over tree
    int nEvents = tree->GetEntries();
  	for ( Int_t j = 0 ; j < nEvents ; j++ ) {
    		tree->GetEntry(j) ;
    		if (0ul == (j % 1000ul)) cout << "Read event " << j << "/" << nEvents << endl;
		
            //Add your cuts
            if(!B_PT>2000)continue;
        
            //Fill histograms
            //...
        
            //Fill output tree
            summary_tree->Fill();
	}
	
	output->Write();
	output->Close();

return 0;
}

