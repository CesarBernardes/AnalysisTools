#include "TFile.h"
#include "TTree.h"

void print_list_average_APV_CMN(){

TFile *file = TFile::Open("output_CalculateCMNonTheFly.root");

TTree *tree = (TTree*)file->Get("hybridAna/Miscellanea/tracker_tree");

tree->SetBranchStatus("*",0);  // disable all branches
tree->SetBranchStatus("vec_APVCMN",1); //activate this
tree->SetBranchStatus("vec_fedID",1); //activate this
tree->SetBranchStatus("vec_lumiS",1); //activate this
tree->SetBranchStatus("vec_apvID",1); //activate this
tree->SetBranchStatus("vec_fedCH",1); //activate this

std::vector<Double_t> *vec_APVCMN=0;
std::vector<Int_t> *vec_fedID=0;
std::vector<Int_t> *vec_lumiS=0;
std::vector<Int_t> *vec_apvID=0;
std::vector<Int_t> *vec_fedCH=0;

tree->SetBranchAddress("vec_APVCMN", &vec_APVCMN);
tree->SetBranchAddress("vec_fedID", &vec_fedID);
tree->SetBranchAddress("vec_lumiS", &vec_lumiS);
tree->SetBranchAddress("vec_apvID", &vec_apvID);
tree->SetBranchAddress("vec_fedCH", &vec_fedCH);


TH1F *h_CMN[440][96][2];///FED Id:50-489, FED Ch:0-95, APV Id:0-1

for(unsigned int i_FEDid=50; i_FEDid<=489; i_FEDid++){
   for(unsigned int i_FEDch=0; i_FEDch<=95; i_FEDch++){
      for(unsigned int i_APVid=0; i_APVid<=1; i_APVid++){
         h_CMN[i_FEDid-50][i_FEDch][i_APVid] = new TH1F(Form("CMN_FEDid%d_FEDch%d_APVid%d",i_FEDid,i_FEDch,i_APVid),Form("CMN_FEDid%d_FEDch%d_APVid%d",i_FEDid,i_FEDch,i_APVid),2000,-1000,1000);	      
      }	      
   }	    
}	


Long64_t n_entries = tree->GetEntriesFast();
///Long64_t n_entries = 10;
for(Long64_t ii=0; ii<n_entries; ii++){

   tree->GetEntry(ii);
   if(ii % 100 == 0) printf("current entry = %lld out of %lld : %.3f %%\n", ii, n_entries, (Double_t)ii / n_entries * 100);

   for(unsigned int i=0; i<vec_APVCMN->size(); i++){
      if((*vec_lumiS)[i]<0)continue; ///IMPORTANT: LS selection (see Control Plots)	   
      for(unsigned int i_FEDid=50; i_FEDid<=489; i_FEDid++){
	 if((*vec_fedID)[i]!=i_FEDid)continue;     
         for(unsigned int i_FEDch=0; i_FEDch<=95; i_FEDch++){	   
            if((*vec_fedCH)[i]!=i_FEDch)continue;     
            for(unsigned int i_APVid=0; i_APVid<=1; i_APVid++){
               if((*vec_apvID)[i]!=i_APVid)continue; 		 
               h_CMN[i_FEDid-50][i_FEDch][i_APVid]->Fill((*vec_APVCMN)[i]);
            }
         }
      }
   }	   
	
}

gSystem->Exec("mkdir -p APV_CMNaverage-Minus127_PbPb2024_RunXXXXXX_2sigma");
std::ofstream APV_CMNaverage_txt("./APV_CMNaverage-Minus127_PbPb2024_RunXXXXXX_2sigma/FEDId_FEDCh_APVId_AverageAPVCMN-Minus127_PbPb2024_RunXXXXXX_2sigma.txt");

//To do control plots
///TCanvas *c1 = new TCanvas();
///c1->cd();

for(unsigned int i_FEDid=50; i_FEDid<=489; i_FEDid++){
   for(unsigned int i_FEDch=0; i_FEDch<=95; i_FEDch++){
      for(unsigned int i_APVid=0; i_APVid<=1; i_APVid++){
         
	 ///h_CMN[0][i_FEDch][i_APVid]->Draw("hist"); //To do control plots

	 //if channel is too low you put a positive number indicating by how much it has to be raised
	 //to that number you add the sigma value
	 /////if(h_CMN[i_FEDid-50][i_FEDch][i_APVid]->GetMean()<127.){
	 if(h_CMN[i_FEDid-50][i_FEDch][i_APVid]->GetMean()<(127.+2.*9.2)){ //checking 2-sigmas - the sigma should be decided by looking directly into the overall spread of the CMN distribution (look directly into the input tree)
	    if(h_CMN[i_FEDid-50][i_FEDch][i_APVid]->GetMean()==0){}//not used APVs
            else{ 	    
	       APV_CMNaverage_txt<<i_FEDid<<" "<<i_FEDch<<" "<<i_APVid<<" "<< (127 - h_CMN[i_FEDid-50][i_FEDch][i_APVid]->GetMean()) + 2.*9.2 <<endl; //checking 2-sigmas
	    }   
	 }else{//no need to shift
	    //APV_CMNaverage_txt<<i_FEDid<<" "<<i_FEDch<<" "<<i_APVid<<" "<<0.0<<endl;//Erik asked toremove these cases	 
	 }	 

	 //To control plots the CMN shifts and check the sigmas
         ///if(i_FEDid==50 && i_FEDch==0 && i_APVid==0)c1->Print("APV_CMNaverage_RunXXXXXX/hist_test.pdf(");
         ///if(!(i_FEDid==50 && i_FEDch==0 && i_APVid==0) && !(i_FEDid==489 && i_FEDch==95 && i_APVid==1))c1->Print("APV_CMNaverage_RunXXXXXX/hist_test.pdf");
         ///if(i_FEDid==489 && i_FEDch==95 && i_APVid==1)c1->Print("APV_CMNaverage_RunXXXXXX/hist_test.pdf)");
         ///c1->Update();
      }	
   }
}
//tree->ResetBranchAddresses(); // detach "tree" from local variables
//delete file; // automatically deletes "file", too

}	
