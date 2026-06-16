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
const int nYBins = 20;
const double yMin = -2.0;
const double yMax =  2.0;
const double dy   = (yMax-yMin)/nYBins;

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
void v1_SP_multispecies_HeHe(const char* infile="/data/sdogra/ampt_helium/inPutFiles/*.dat",const char* sysTag = "HeHe"){
  // ---------- accumulators ----------
  double num[NSPEC][nCent][nYBins] = {};
  double num_sq[NSPEC][nCent][nYBins] = {};
  int    count[NSPEC][nCent][nYBins] = {};
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
      double QAnorm=0.0;//Cesar
      double QBnorm=0.0;//Cesar

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
        // BUG 2 FIX — centrality from 3<|η|<5, fully separated from EP subevents
	if(fabs(eta)>3.0 && fabs(eta)<5.0 && pt>0.02){
          if (abs(pid) == 211 || abs(pid) == 321 || abs(pid) == 2212)nchRef++;
        }
        double phi = atan2(py,px);

        // subevents — unchanged
        ///if(eta > 2.4){ QxA += cos(phi); QyA += sin(phi); }
        ///if(eta < -2.4){ QxB += cos(phi); QyB += sin(phi); }
	///Cesar - begin
        if(eta > 3.0 && eta < 5.0){ QxA += pt*cos(phi); QyA += pt*sin(phi); QAnorm += pt; }
        if(eta < -3.0 && eta > -5.0){ QxB += pt*cos(phi); QyB += pt*sin(phi); QBnorm += pt; }
	//if(eta > 3.0 && eta < 5.0){ QxA += pt*cos(2.*phi); QyA += pt*sin(2.*phi); QAnorm += pt; }//test v2
	//if(eta < -3.0 && eta > -5.0){ QxB += pt*cos(2.*phi); QyB += pt*sin(2.*phi); QBnorm += pt; }//test v2
        ///Cesar - end

        // analysis region
        if(fabs(eta)>2.4) continue;
	//if(pt<1.0)continue;//Cesar: test v2

        // rapidity
        double mass = 0.139;
        if(pid==2212) mass = 0.938;
        if(abs(pid)==321) mass = 0.494;

        double Ecalc = sqrt(px*px + py*py + pz*pz + mass*mass);
        double yrap = compute_y(Ecalc,pz);

        if(yrap<yMin || yrap>=yMax) continue;

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
      int centEdges[10]={0,13,22,28,37,47,62,85,105,100000};
      int centBin=-1;
      for(int i=0;i<9;i++){
        if(nchRef>=centEdges[i] && nchRef<centEdges[i+1]){
          centBin=8-i;
          break;
        }
      }

      if(centBin<0) continue;

      // -------- denominator --------
      ///double AB = QxA*QxB + QyA*QyB;
      if(QAnorm<=0.0 || QBnorm<=0.0) continue;
      double AB = (QxA*QxB + QyA*QyB)/(QAnorm*QBnorm); //Cesar
      if(AB<=0) continue;
      den[centBin] += AB;
      nevts[centBin]++;//Cesar

      // -------- numerator --------
      for(auto &t:tracks){
        int yBin = (int)((t.y - yMin)/dy);
        if(yBin<0 || yBin>=nYBins) continue;
	double ux = cos(t.phi);
	double uy = sin(t.phi);
	//	double ux = cos(2.*t.phi);//test v2
	//      double uy = sin(2.*t.phi);//test v2
        // BUG 1 FIX — use opposite subevent based on rapidity sign
        // forward track (y>0) → correlate with backward subevent B
        // backward track (y<0) → correlate with forward subevent A
        double uQ = 0.0;
        if(t.y > 0)uQ = (ux*QxB + uy*QyB)/QBnorm; ///Cesar
        else uQ = (ux*QxA + uy*QyA)/QAnorm; ///Cesar
        for(int s=0;s<NSPEC;s++){
          if(t.is[s]){
            num[s][centBin][yBin] += uQ;
            num_sq[s][centBin][yBin] += uQ*uQ;
            count[s][centBin][yBin]++;
          }
        }
      }
    }
    in.close();
  }
  
  cout<<"========== ALL FILES DONE =========="<<endl;
  cout<<"Total processed events = "<<eventCounter<<endl;
  
  // ================= OUTPUT =================
  TString outName = Form("v1_SP_multispecies_%s.root", sysTag);
  TFile* fout = TFile::Open(outName,"RECREATE");

  TH1D *hSlope[NSPEC];
  for(int s=0;s<NSPEC;s++)hSlope[s] = new TH1D(Form("dv1dy_%s",spName[s]),"",nCent,0,nCent);
  for(int c=0;c<nCent;c++){
    double denom = (den[c]>0) ? sqrt(den[c]/nevts[c]) : 0; //Cesar
    for(int s=0;s<NSPEC;s++){
      TH1D *h = new TH1D(Form("v1_%s_cent%d",spName[s],c),"", nYBins,yMin,yMax);
      for(int i=0;i<nYBins;i++){
        if(count[s][c][i]>0 && denom>0){
          double mean = num[s][c][i] / count[s][c][i];
          double var  = num_sq[s][c][i]/count[s][c][i] - mean*mean;
          if(var<0) var=0;
          double v1 = mean / denom;
          double err = sqrt(var/count[s][c][i]) / denom;
          h->SetBinContent(i+1,v1);
          h->SetBinError(i+1,err);
        }
      }

      TH1D *hOdd = (TH1D*)h->Clone(Form("v1odd_%s_cent%d",spName[s],c));

      // BUG 3 FIX — propagate errors into hOdd
      for(int i=1;i<=nYBins;i++){
        int j = nYBins-i+1;
        double odd = 0.5*(h->GetBinContent(i) - h->GetBinContent(j));
        double err = 0.5*sqrt(pow(h->GetBinError(i),2) + pow(h->GetBinError(j),2));
        hOdd->SetBinContent(i,odd);
        hOdd->SetBinError(i,err);
      }

      TF1 *f = new TF1(Form("fit_%s_cent%d",spName[s],c),"[0]*x",-0.5,0.5);
      hOdd->Fit(f,"Q");

      double slope = f->GetParameter(0);
      double slopeErr = f->GetParError(0);

      hSlope[s]->SetBinContent(c+1,slope);
      hSlope[s]->SetBinError(c+1,slopeErr);

      h->Write();
      hOdd->Write();
      f->Write();
      cout<<"Cent "<<c<<" "<<spName[s] <<" dv1/dy = "<<slope<<" ± "<<slopeErr<<endl;
    }
  }

  for(int s=0;s<NSPEC;s++) hSlope[s]->Write();
  fout->Close();
  cout<<"DONE → v1_SP_multispecies.root"<<endl;
}
