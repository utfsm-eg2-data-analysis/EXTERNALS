#include "TROOT.h"
#include "TH2F.h"
#include "THStack.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TTree.h"
#include "TCut.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TBenchmark.h"
#include "TGraphErrors.h"
#include "TColor.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sys/stat.h>
#include <getopt.h>
#include <fstream>
#include <string>
#include <iomanip>
#include <map>
#include "TH1.h"
#include "TAxis.h"
#include "TEventList.h"
#include <TLatex.h>
#include <stdlib.h>
#include <TFile.h>
#include <TH1F.h>
#include <Riostream.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TFrame.h>
#include <TChain.h>
#include <TMarker.h>
#include <TLegend.h>
#include "TLine.h"
#include <TPaveText.h>

int main(int argc , char **argv){

const Int_t NRGBs = 5;const Int_t NCont = 255;
Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
gStyle->SetNumberContours(NCont);


  TString Metal        = (TString)argv[1];    //"C";
  TString Liq_or_Solid = (TString)argv[2];    //"liq", "sol";
  TString VC           = (TString)argv[3];    //"H" , "RD", "KPA"
  TString dataLoc, dataLocSim, file, fileSim;
  dataLoc    = "/work/smoran/data/externalFiles/";
  dataLocSim = "/work/smoran/simul/externalFiles/";

  file    = Form( "%s_external.root" , (const char*)Metal);
  if (Liq_or_Solid=="sol")  fileSim = Form( "%s.root" , (const char*)Metal);
  else if (Liq_or_Solid=="liq")  fileSim = "D2.root";

  Double_t  delta_q2 , delta_xb ;//, *v_q2 , *v_xb;
  
  Double_t limQ2=1.0;
  Double_t limW2=4.0; //1.8*1.8=3.24
  
  
  const Double_t  q2_min = 1.0; 
  const Double_t  q2_max = 4.; 
  const Double_t  xb_min = 0.2; 
  const Double_t  xb_max = 0.6;
  const Int_t     n_q2   = 60 ;  
  const Int_t     n_xb   = 44 ;

  delta_q2  = (q2_max-q2_min)/n_q2;   
  delta_xb  = (xb_max-xb_min)/n_xb;

//  Double_t v_q2[13]={1.0,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4};
//  Double_t v_xb[12]={0.12,0.16,0.2,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56};

Double_t v_q2[61]={1.00 ,1.05 , 1.10 , 1.15 , 1.20 , 1.25 , 1.30 , 1.35 , 1.40 , 1.45 , 1.50 , 1.55 , 1.60 , 1.65 , 1.70 , 1.75 , 1.80 , 1.85 , 1.90, 1.95 , 2.00 , 2.05 , 2.10 , 2.15 , 2.20 , 2.25 , 2.30 , 2.35 , 2.40 , 2.45 , 2.50 , 2.55 , 2.60 , 2.65 , 2.70 , 2.75 , 2.80, 2.85 , 2.90 , 2.95 , 3.00 , 3.05 , 3.10 , 3.15 , 3.20 , 3.25 , 3.30 , 3.35 , 3.40 , 3.45 , 3.50 , 3.55 , 3.60 , 3.65 , 3.70, 3.75 , 3.80 , 3.85 , 3.90 , 3.95 , 4.00 };

Double_t v_xb[45] = {0.12 , 0.13 , 0.14 , 0.15 , 0.16 , 0.17 , 0.18 , 0.19 , 0.20 , 0.21 , 0.22 , 0.23 , 0.24 , 0.25 , 0.26 , 0.27 , 0.28 , 0.29 , 0.30 , 0.31 , 0.32 , 0.33 , 0.34 , 0.35 , 0.36 , 0.37 , 0.38 , 0.39 , 0.40 , 0.41 , 0.42 , 0.43 , 0.44 , 0.45 , 0.46 , 0.47 , 0.48 , 0.49 , 0.50 , 0.51 , 0.52 , 0.53 , 0.54 , 0.55 , 0.56 };

  std::ofstream ofs (Form("centroids_%s_%s_%s_vc.txt", (const char*)Metal , (const char*)Liq_or_Solid, (const char*)VC), std::ofstream::out);
  std::ofstream ofs2(Form("centroids_%s_%s_%s_vc_for_plot.txt",(const char*)Metal,(const char*)Liq_or_Solid,(const char*)VC),std::ofstream::out);
  
  std::ofstream ofs3(Form("clas_kin_%s_%s_%s_vc.txt",(const char*)Metal,(const char*)Liq_or_Solid,(const char*)VC),std::ofstream::out);

  TCut cut_1, cut_2, cut_1Sim, cut_2Sim,liq , sol,solSim,liqSim ,  dis_cut, cut_w, cut_q2 , cut_Yb, target_cut, target_cutSim;

  if (VC=="HH")     { liq = Form("TargType==1 && Q2>%f && (0.938272*0.938272 + 2*0.938272*Nu - Q2)>%f", limQ2, limW2); 
  		      sol = Form("TargType==2 && Q2>%f && (0.938272*0.938272 + 2*0.938272*Nu - Q2)>%f", limQ2, limW2);}
		      
  else if (VC=="RD"){ liq = Form("VC_RD==1 && Q2>%f && (0.938272*0.938272 + 2*0.938272*Nu - Q2)>%f", limQ2, limW2); 
  		      sol = Form("VC_RD==2 && Q2>%f && (0.938272*0.938272 + 2*0.938272*Nu - Q2)>%f", limQ2, limW2);
		      liqSim = Form("TargType==1 && Q2>%f && (0.938272*0.938272 + 2*0.938272*Nu - Q2)>%f", limQ2, limW2);
		      solSim = Form("TargType==2 && Q2>%f && (0.938272*0.938272 + 2*0.938272*Nu - Q2)>%f", limQ2, limW2); }


  if (Liq_or_Solid == "liq" ){ target_cut = liq; target_cutSim = liqSim; } else if (Liq_or_Solid == "sol") {target_cut = sol;target_cutSim = solSim;}
  
  TFile   *f ,*fSim, *fout; TH2F* h;
  Double_t meanx , meany;
  fout = new TFile(  Form("out_%s_%s_%s_vc.root", (const char*)Metal ,(const char*)Liq_or_Solid , (const char*)VC)  ,  "recreate");
  f       = new TFile(dataLoc + file, "READ");
  fSim    = new TFile(dataLocSim + fileSim, "READ");
  
  TChain *t_elec    = new TChain();t_elec->Add(Form(dataLoc + file + "/ntuple_data_electrons" ,(const char*)Metal));
  TChain *t_elecRec = new TChain();t_elecRec->Add(Form(dataLocSim + fileSim + "/ntuple_accept_electrons" ,(const char*)Metal));
  TChain *t_elecThr = new TChain();t_elecThr->Add(Form(dataLocSim + fileSim + "/ntuple_thrown_electrons" ,(const char*)Metal));

  ofs<<"('Xb','Q2') (";
  
  ofs3<<"RUNPLAN for CLAS hadronization experiment\n";
  ofs3<<" inclusive:\n\n\n\n"; 
  ofs3<<" E     Ep    theta     W     y     x      Q2\n";
    
    for (int i = 0 ; i < n_q2 ; i++){
    cut_1 = Form("Q2>%f && Q2<%f && Q2>%f && (0.938272*0.938272 + 2*0.938272*Nu - Q2)>%f && Nphe>(Sector==0||Sector==1)*25+(Sector==2)*26+(Sector==3)*21+(Sector==4||Sector==5)*28  ", v_q2[i], v_q2[i+1], limQ2, limW2);
    cut_1Sim = Form("Q2>%f && Q2<%f && Q2>%f && (0.938272*0.938272 + 2*0.938272*Nu - Q2)>%f", v_q2[i], v_q2[i+1], limQ2, limW2);  
    for (int j = 0 ; j < n_xb  ; j++){
    
      cut_2 = Form("Xb>%f && Xb<%f", v_xb[j], v_xb[j+1]); 
      TCut cut, cutSim, cutSimThr;
      cut       = cut_1 && cut_2 && target_cut ;
      cutSim    = cut_1Sim && cut_2 && target_cutSim ;
      cutSimThr = cut_1Sim && cut_2;
      
t_elec->Draw((const char*)Form("Q2:Xb>>h(30,%f,%f,30,%f,%f)",v_xb[j],v_xb[j+1],v_q2[i],v_q2[i+1]),cut,"goff"); 
t_elecRec->Draw((const char*)Form("Q2:Xb>>hrec(30,%f,%f,30,%f,%f)",v_xb[j], v_xb[j+1], v_q2[i], v_q2[i+1]),cutSim ,"goff"); 
t_elecThr->Draw((const char*)Form("Q2:Xb>>hthr(30,%f,%f,30,%f,%f)",v_xb[j], v_xb[j+1], v_q2[i], v_q2[i+1]),cutSimThr,"goff"); 
  
   TH2F* h    = (TH2F*)gDirectory->GetList()->FindObject("h");
   TH2F* hrec = (TH2F*)gDirectory->GetList()->FindObject("hrec");
   TH2F* hthr = (TH2F*)gDirectory->GetList()->FindObject("hthr");
      
   TH2F *htmpData;  htmpData  = (TH2F*)h->Clone("tmpData");
   TH2F *htmpRec;   htmpRec   = (TH2F*)hrec->Clone("tmpRec");
   TH2F *htmpGen;   htmpGen   = (TH2F*)hthr->Clone("tmpGen");

htmpRec->Divide(htmpGen);
htmpData->Divide(htmpRec);

      meanx = htmpData->GetMean(1); 
      meany = htmpData->GetMean(2);
      
      if(j == n_xb-1 && i== n_q2-1){ofs<< "("<<meanx<<","<<meany<<")";} 
      else {ofs<< "("<<meanx<<","<<meany<<")"<< ",";}
      
      ofs2<<meanx<<" "<<meany<< "\n";
      
      
      if (meanx!=0 && meany!=0){
      Double_t xp1 = 5.014 - (meany/(2*meanx*0.93827)) ;
      Double_t yp1 = 2*(asin(sqrt(meany/( 4*5.014*xp1)))*(180./3.141592))   ;
      
      //ofs3<<"5.014 "<< "0"<< xp1 << " " << yp1<< "\n";
      ofs3<<"5.014 "<< "0"<< std::fixed << std::setprecision(3)  <<xp1 << " " <<std::fixed << std::setprecision(4)<< yp1<< "\n";
      }
      h->Delete();
      hrec->Delete();
      hthr->Delete();
      htmpData->Delete();
      htmpRec->Delete();
      htmpGen->Delete();
      
    } 
  }
  ofs<<")";
  ofs.close(); 
  ofs2.close();
  ofs3.close();
  TCanvas *c = new TCanvas("c", "canvas", 500, 500);
  t_elec->Draw("Q2:Xb>>h2(200,0.1,0.6,200,0.8,4.2)",  target_cut , "goff");
  TH2F *h2d = (TH2F*) gDirectory->GetList()->FindObject("h2");
  h2d->SetMarkerColor(5);
  h2d->SetStats(0);
  h2d->SetTitle(Form("%s , data %s Target, %s Vertex Cuts", (const char*)Metal ,(const char*) Liq_or_Solid , (const char*)VC ));
  gStyle->SetTitleFont(22,"t"); 
  h2d->SetTitleFont(22, "X");
  h2d->GetXaxis()->SetLabelFont(22);
  h2d->SetTitleFont(22, "Y");
  h2d->GetYaxis()->SetLabelFont(22);  
  h2d->GetXaxis()->SetTitle("Xb");
  h2d->GetYaxis()->SetTitle("Q2");
  h2d->GetYaxis()->CenterTitle(true) ; 
  h2d->GetXaxis()->CenterTitle(true);
  h2d->Draw("colz");
  TLine *lin_H = new TLine();  
  lin_H->SetLineWidth(2); 
  lin_H->SetLineStyle(2);
  TLine *lin_V = new TLine();  
  lin_V->SetLineWidth(2); 
  lin_V->SetLineStyle(2);
  Double_t ll, ii; 
  lin_H->SetLineColor(kGreen+4);
  lin_V->SetLineColor(kGreen+4);
  for (int i =0 ; i<n_q2+1 ; i++){ll = v_q2[i];lin_H->DrawLine(0.1, ll,0.59, ll);}
  for (int i =0 ; i<n_xb+1 ; i++){ii = v_xb[i];lin_V->DrawLine(ii, 0.9, ii,4.1) ;}
TGraph *g=new TGraph(Form("centroids_%s_%s_%s_vc_for_plot.txt",(const char*)Metal,(const char*)Liq_or_Solid,(const char*)VC)  ); 
  g->SetMarkerStyle(kFullDotLarge); 
  g->SetMarkerSize(0.9);
  g->Draw("psame");
  fout->cd();
  c->Write("centroid_plot", TObject::kOverwrite);
  delete c; 
  delete h2d;
  f->Close(); 
  fout->Close(); 
  delete fout ; 
  delete  f; 
  ofs.close();
  //delete v_q2; 
  //delete v_xb; 
  return 0;
}




