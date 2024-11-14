{

//TFile *fout = new TFile("output.root","RECREATE");

auto froot = TFile::Open("output_CalculateCMNonTheFly_387892.root");
TTree *t0 = (TTree *)froot->Get("hybridAna/Miscellanea/tracker_tree");

std::vector<double> *vec_APVCMN=0;
std::vector<int> *vec_fedID=0;
std::vector<int> *vec_fedCH=0;
std::vector<int> *vec_apvID=0;

t0->SetBranchAddress("vec_APVCMN",&vec_APVCMN);
t0->SetBranchAddress("vec_fedID",&vec_fedID);
t0->SetBranchAddress("vec_fedCH",&vec_fedCH);
t0->SetBranchAddress("vec_apvID",&vec_apvID);


auto h = new TH1F("hist","hist",2400,-1200,1200);
//h->SetTitle("FED 106 - Selected APVs (10% higher rates of bad APVs)");
h->SetTitle("FED 112 - Selected APVs (10% higher rates of bad APVs)");
h->GetXaxis()->SetTitle("CM");
h->GetYaxis()->SetTitle("Entries / bin");

TTree *t1 = new TTree("ntuple","data from ascii file");

//Long64_t nlines = t1->ReadFile("APV_BAD_APV_Rates/FEDId_FEDCh_APVId_APVsWithTwoSigmasOutInMisbehaving_FED106.txt","fedid:fedch:apvid:cmn"); //FED106
Long64_t nlines = t1->ReadFile("APV_BAD_APV_Rates/FEDId_FEDCh_APVId_APVsWithTwoSigmasOutInMisbehaving_FED112.txt","fedid:fedch:apvid:cmn"); //FED112

printf(" found %lld points\n",nlines);
//t1->Draw("fedid","");

Float_t fedid,fedch,apvid;
t1->SetBranchAddress("fedid",&fedid);
t1->SetBranchAddress("fedch",&fedch);
t1->SetBranchAddress("apvid",&apvid);


Int_t nevents0 = t0->GetEntries();

for (Int_t i=0;i<nevents0;i++) {

   t0->GetEntry(i);
   for(Int_t ii=0; ii<(*vec_APVCMN).size(); ii++){

      Int_t nevents1 = t1->GetEntries();
      for (Int_t j=0;j<nevents1;j++){
         t1->GetEntry(j);

         if((*vec_fedID)[ii]==fedid && (*vec_fedCH)[ii]==fedch && (*vec_apvID)[ii]==apvid){
            float aux_cmn = (*vec_APVCMN)[ii];
            h->Fill(aux_cmn);
            break;
         }

      }
   }

}


TFile *fout = new TFile("output.root","RECREATE");

h->Write();

}
