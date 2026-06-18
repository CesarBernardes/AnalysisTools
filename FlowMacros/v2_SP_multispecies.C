#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TRegexp.h"
#include "TString.h"

#include "TH1D.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"

using namespace std;

// ================= SETTINGS =================
const int nCent = 9;
const int nPtBins = 10;
const double ptMin = 0.5;
const double ptMax = 5.0;
const double dpt   = (ptMax-ptMin)/nPtBins;

enum Species { PROT=0, PION=1, KAON=2, CH=3, NSPEC=4 };
const char* spName[NSPEC] = {"prot","pion","kaon","ch"};

// ================= ETA =================
double compute_eta(double px, double py, double pz){
  double pt = sqrt(px*px + py*py);
  if(pt==0) return 0.0;
  double theta = atan2(pt,pz);
  if(theta<0) theta += TMath::Pi();
  return -log(tan(theta/2.0));
}

// ================= RAPIDITY =================
double compute_y(double E, double pz){
  double num = E + pz;
  double den = E - pz;
  if(num<=0 || den<=0) return 999.0;
  return 0.5 * log(num/den);
}

// ================= MAIN =================
void v2_SP_multispecies(const char* infile="OO_ampt_100K.dat",const char* sysTag = "OO"){
  // ---------- accumulators ----------
  double num[NSPEC][nCent][nPtBins] = {};
  double num_sq[NSPEC][nCent][nPtBins] = {};
  int    count[NSPEC][nCent][nPtBins] = {};
  double den[nCent] = {};
  int nevts[nCent] = {};//Cesar
  // ---------- file handling ----------
  TString pattern(infile);
  TString dir = gSystem->DirName(pattern);
  TString base = gSystem->BaseName(pattern);

  vector<TString> files;
  TSystemDirectory directory("dir",dir);
  TList* fileList = directory.GetListOfFiles();
  if(!fileList){ cout<<"Directory error\n"; return; }

  TIter next(fileList);
  TSystemFile* file;
  TRegexp re(base,true);

  while((file=(TSystemFile*)next())){
    TString fname=file->GetName();
    if(!file->IsDirectory() && fname.Index(re)!=kNPOS)files.push_back(dir+"/"+fname);
  }
  
  cout<<"Found "<<files.size()<<" files\n";

  // ---------- variables ----------
  int eventID,dummy,n_particles;
  int npartP,npartT,nElaP,nInElaP,nElaT,nInElaT;
  float b,unused;

  int pid;
  float px,py,pz,E,x,y,z,t;

  // ---------- counters ----------
  long long eventCounter = 0;

  // ================= LOOP =================
  for(auto &fname:files){
    ifstream in(fname.Data());
    if(!in.is_open()){cout<<"Cannot open "<<fname<<endl; continue;}
    cout<<"========== Reading file: "<<fname<<" =========="<<endl;
    while(true){
      if(!(in >> eventID >> dummy >> n_particles >> b >> npartP >> npartT >> nElaP >> nInElaP >> nElaT >> nInElaT >> unused)) break;
      // ===== event counter =====
      eventCounter++;
      if (eventCounter % 1000 == 0) cout << "[INFO] Events processed: " << eventCounter << "\r" << flush;
        
      int nchRef = 0;
      double QxA=0.0, QyA=0.0;
      double QxB=0.0, QyB=0.0;
      double QAnorm=0.0;
      double QBnorm=0.0;

      struct Track{
        double px,py,pz,E,phi,y;
        bool is[NSPEC];
      };
      vector<Track> tracks;
      // -------- PARTICLES --------
      for(int i=0;i<n_particles;i++){
        if(!(in >> pid >> px >> py >> pz >> E >> x >> y >> z >> t))break;
        double pt = sqrt(px*px + py*py);
        if(pt < 0.02) continue;
        double eta = compute_eta(px,py,pz);

	if(fabs(eta)>3.0 && fabs(eta)<5.0 && pt>0.02){
          if (abs(pid) == 211 || abs(pid) == 321 || abs(pid) == 2212)nchRef++;
        }
        double phi = atan2(py,px);

        // subevents
	if(eta > 3.0 && eta < 5.0){ QxA += pt*cos(2.*phi); QyA += pt*sin(2.*phi); QAnorm += pt; }
	if(eta < -3.0 && eta > -5.0){ QxB += pt*cos(2.*phi); QyB += pt*sin(2.*phi); QBnorm += pt; }

        // analysis region
        if(fabs(eta)>1.5) continue;

        // rapidity
        double mass = 0.139;
        if(abs(pid)==2212) mass = 0.938;
        if(abs(pid)==321) mass = 0.494;

        double Ecalc = sqrt(px*px + py*py + pz*pz + mass*mass);
        double yrap = compute_y(Ecalc,pz);

        if(pt<ptMin || pt>=ptMax) continue;

        Track tr;
        tr.px=px; tr.py=py; tr.pz=pz;
        tr.E=Ecalc; tr.phi=phi; tr.y=yrap;

        for(int s=0;s<NSPEC;s++) tr.is[s]=false;

        if(abs(pid)==2212) tr.is[PROT]=true;
        if(abs(pid)==211) tr.is[PION]=true;
        if(abs(pid)==321) tr.is[KAON]=true;
        if(tr.is[PROT] || tr.is[PION] || tr.is[KAON]) tr.is[CH]=true;
        tracks.push_back(tr);
      }

      // -------- centrality --------
      //int centEdges[nCent+1]={0,13,22,28,37,47,62,85,105,100000}; //HeHe
      int centEdges[nCent+1]={0, 22, 54, 82, 120, 172, 243, 341, 411, 100000}; //OO
      int centBin=-1;
      for(int i=0;i<nCent;i++){
        if(nchRef>=centEdges[i] && nchRef<centEdges[i+1]){
          centBin=8-i;
          break;
        }
      }

      if(centBin<0) continue;

      // -------- denominator --------
      if(QAnorm<=0.0 || QBnorm<=0.0) continue;
      double AB = (QxA*QxB + QyA*QyB)/(QAnorm*QBnorm);
      if(AB<=0) continue;
      den[centBin] += AB;
      nevts[centBin]++;

      // -------- numerator --------
      for(auto &t:tracks){
	double aux_pt = sqrt(t.px*t.px + t.py*t.py);
	int ptBin = (int)((aux_pt - ptMin)/dpt);
	if(ptBin<0 || ptBin>=nPtBins) continue;
	double ux = cos(2.*t.phi);
        double uy = sin(2.*t.phi);
        // Use opposite subevent based on rapidity sign
        // forward track (y>0) → correlate with backward subevent B
        // backward track (y<0) → correlate with forward subevent A
        double uQ = 0.0;
        if(t.y > 0)uQ = (ux*QxB + uy*QyB)/QBnorm;
        else uQ = (ux*QxA + uy*QyA)/QAnorm;
        for(int s=0;s<NSPEC;s++){
          if(t.is[s]){
            num[s][centBin][ptBin] += uQ;
            num_sq[s][centBin][ptBin] += uQ*uQ;
            count[s][centBin][ptBin]++;
          }
        }
      }
    }
    in.close();
  }
  
  cout<<"========== ALL FILES DONE =========="<<endl;
  cout<<"Total processed events = "<<eventCounter<<endl;
  
  // ================= OUTPUT =================
  TString outName = Form("v2_SP_multispecies_%s.root", sysTag);
  TFile* fout = TFile::Open(outName,"RECREATE");

  for(int c=0;c<nCent;c++){
    double denom = (den[c]>0) ? sqrt(den[c]/nevts[c]) : 0;
    for(int s=0;s<NSPEC;s++){
      TH1D *h = new TH1D(Form("v2_%s_cent%d",spName[s],c),"", nPtBins,ptMin,ptMax);
      for(int i=0;i<nPtBins;i++){
        if(count[s][c][i]>0 && denom>0){
          double mean = num[s][c][i] / count[s][c][i];
          double var  = num_sq[s][c][i]/count[s][c][i] - mean*mean;
          if(var<0) var=0;
          double v2 = mean / denom;
          double err = sqrt(var/count[s][c][i]) / denom;
          h->SetBinContent(i+1,v2);
          h->SetBinError(i+1,err);
        }
      }

      h->Write();
    }
  }

  fout->Close();
  cout<<"DONE → v2_SP_multispecies.root"<<endl;
}
