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
  TString dataLoc, file;
  dataLoc = "/work/smoran/data/externalFiles/";
  file = Form( "%s_external.root" , (const char*)Metal);
  Double_t  delta_q2 , delta_xb ;//, *v_q2 , *v_xb;
  
  Double_t limQ2=1.;
  Double_t limW2=4.; //1.8*1.8=3.24
  //HERE I'M USING NU INSTEAD OF XB BUT I'M NOT CHANGING THE LABELS.
  
  const Double_t  q2_min = 1.0; 
  const Double_t  q2_max = 4.5; 
  const Double_t  xb_min = 2.0; 
  const Double_t  xb_max = 4.26;
  const Int_t     n_q2   = 3;  
  const Int_t     n_xb   = 3;

  delta_q2  = (q2_max-q2_min)/n_q2;   
  delta_xb  = (xb_max-xb_min)/n_xb;

   Double_t v_q2[4]={1.0,1.33,1.76,4.1};
   Double_t v_xb[4]={2.2,3.2,3.73,4.25};


 // Double_t v_q2[7]={1., 1.17 , 1.33, 1.51 , 1.75 , 2.12, 4.};
 // Double_t v_xb[6]={ 0.12, 0.19, 0.23, 0.27 , 0.33, 0.57 };
  
  //Double_t v_q2[2]={1.0,4};

  std::ofstream ofs (Form("centroids_%s_%s_%s_vc.txt", (const char*)Metal , (const char*)Liq_or_Solid, (const char*)VC), std::ofstream::out);
  std::ofstream ofs2(Form("centroids_%s_%s_%s_vc_for_plot.txt",(const char*)Metal,(const char*)Liq_or_Solid,(const char*)VC),std::ofstream::out);
    std::ofstream ofs22(Form("centroids_%s_%s_%s_vc_for_plot2.txt",(const char*)Metal,(const char*)Liq_or_Solid,(const char*)VC),std::ofstream::out);
  std::ofstream ofs3(Form("clas_kin_%s_%s_%s_vc.txt",(const char*)Metal,(const char*)Liq_or_Solid,(const char*)VC),std::ofstream::out);
  
  //for(Int_t i = 0; i < n_q2+1; i++){  if(i == 0) v_q2[i] = q2_min;  else v_q2[i] = v_q2[i-1] + delta_q2; }
  //for(Int_t i = 0; i < n_xb+1; i++){  if(i == 0) v_xb[i] = xb_min;  else v_xb[i] = v_xb[i-1] + delta_xb; }

  TCut cut_1, cut_2 , liq , sol ,  dis_cut, cut_w, cut_q2 , cut_Yb, target_cut;
  /*
  if (VC=="HH")     { liq = "TargType==1"; sol = "TargType==2";}
  else if (VC=="RD"){ liq = "VC_RD==1"   ; sol = "VC_RD==2"   ;}
*/

  if (VC=="HH")     { liq = Form("TargType==1 && Q2>%f && (0.938272*0.938272 + 2*0.938272*Nu - Q2)>%f", limQ2, limW2); 
  		      sol = Form("TargType==2 && Q2>%f && (0.938272*0.938272 + 2*0.938272*Nu - Q2)>%f", limQ2, limW2);}
  else if (VC=="RD"){ liq = Form("VC_RD==1 && Q2>%f && (0.938272*0.938272 + 2*0.938272*Nu - Q2)>%f", limQ2, limW2); 
  		      sol = Form("VC_RD==2 && Q2>%f && (0.938272*0.938272 + 2*0.938272*Nu - Q2)>%f", limQ2, limW2); }



  if (Liq_or_Solid == "liq" ){ target_cut = liq; } else if (Liq_or_Solid == "sol") {target_cut = sol;}
  TFile   *f , *fout; ;TH2F* h;
  Double_t meanx , meany;
  fout = new TFile(  Form("out_%s_%s_%s_vc.root", (const char*)Metal ,(const char*)Liq_or_Solid , (const char*)VC)  ,  "recreate");
  f    = new TFile(dataLoc + file, "READ");
  TChain *t_elec = new TChain();t_elec->Add(Form(dataLoc + file + "/ntuple_data_electrons" ,(const char*)Metal));
 // t_elec->SetMaxEntryLoop(2000);
  //t_elec->Draw(">>list_dis","", "goff"); t_elec->SetEventList((TEventList*) gDirectory->Get("list_dis"));
  
  ofs<<"('Xb','Q2') (";
  
  ofs3<<"RUNPLAN for CLAS hadronization experiment\n";
  ofs3<<" inclusive:\n\n\n\n"; 
  ofs3<<" E     Ep    theta     W     y     x      Q2\n";
  
  
    for (int i = 0 ; i < n_q2 ; i++){
    cut_1 = Form("Q2>%f && Q2<%f && Q2>%f && (0.938272*0.938272 + 2*0.938272*Nu - Q2)>%f && Nphe>(Sector==0||Sector==1)*25+(Sector==2)*26+(Sector==3)*21+(Sector==4||Sector==5)*28 && TMath::Abs(YC)<1.4  ", v_q2[i], v_q2[i+1], limQ2, limW2);
    for (int j = 0 ; j < n_xb  ; j++){
    
      cut_2 = Form("Nu>%f && Nu<%f", v_xb[j], v_xb[j+1]); 
      TCut cut;
      cut = cut_1 && cut_2 && target_cut ;
      t_elec->Draw("Q2:Nu>>h(5,2.0,4.6,5,0.8,4.2)" , cut , "goff"); 
      TH2F* h = (TH2F*)gDirectory->GetList()->FindObject("h");
      //h->SetName((const char*) Form("h_%s_elec_%s_%d_%d", (const char*)Metal, (const char*)Liq_or_Solid , i,j ));
      //fout->cd();h->Write("" ,TObject::kOverwrite);
      
      meanx = h->GetMean(1); 
      meany = h->GetMean(2);
      
      // if (meanx!=0 && meany!=0){
      if(j == n_xb-1 && i== n_q2-1){ofs<< "("<<meanx<<","<<meany<<")";} 
      else {ofs<< "("<<meanx<<","<<meany<<")"<< ",";}
      //}
      ofs2<<meanx<<" "<<meany<< "\n";
      
      
      
      
      
      if (meanx!=0 && meany!=0){
      Double_t xp1 = 5.014 - (meany/(2*meanx*0.938272)) ;
      Double_t yp1 = 2*(asin(sqrt(meany/( 4*5.014*xp1)))*(180./3.141592))   ;
      
      //ofs3<<"5.014 "<< "0"<< xp1 << " " << yp1<< "\n";
      ofs22<<meanx<<" "<<yp1<< "\n";
      
      
      ofs3<<"5.014 "<< "0"<< std::fixed << std::setprecision(3)  <<xp1 << " " <<std::fixed << std::setprecision(4)<< yp1<< "\n";
      }else {ofs22<< "0 0 \n";}

      h->Delete();
    } 
  }
  ofs<<")";
  ofs.close(); 
  ofs2.close();
  ofs22.close();
  ofs3.close();
  TCanvas *c = new TCanvas("c", "canvas", 500, 500);
  t_elec->Draw("Q2:Nu>>h2(200,2.0,4.46,200,0.8,4.2)",  target_cut , "goff");
  TH2F *h2d = (TH2F*) gDirectory->GetList()->FindObject("h2");
  h2d->SetMarkerColor(5);
  h2d->SetStats(0);
  h2d->SetTitle(Form("%s , data %s Target, %s Vertex Cuts", (const char*)Metal ,(const char*) Liq_or_Solid , (const char*)VC ));
  gStyle->SetTitleFont(22,"t"); 
  h2d->SetTitleFont(22, "X");
  h2d->GetXaxis()->SetLabelFont(22);
  h2d->SetTitleFont(22, "Y");
  h2d->GetYaxis()->SetLabelFont(22);  
  h2d->GetXaxis()->SetTitle("Nu");
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
  for (int i =0 ; i<n_q2+1 ; i++){ll = v_q2[i];lin_H->DrawLine(2.12, ll,4.4, ll);}
  for (int i =0 ; i<n_xb+1 ; i++){ii = v_xb[i];lin_V->DrawLine(ii, 0.9, ii,4.1) ;}
  TGraph *g = new TGraph(Form("centroids_%s_%s_%s_vc_for_plot.txt",(const char*)Metal,(const char*)Liq_or_Solid,(const char*)VC)  ); 
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




