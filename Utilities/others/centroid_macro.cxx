
{


TStyle *seba=new TStyle("seba", "Seba's Style");
seba->SetHistLineWidth(3);
seba->SetTitleFont(22,"");seba->SetLabelFont(22,"XYZ");seba->SetTitleFont(22,"XYZ");seba->SetCanvasBorderMode(0);
seba->SetTitleBorderSize(0);seba->SetTitleOffset(0.82,"yx");
seba->SetStatFont(22);seba->SetLabelSize(0.05,"x");
seba->SetGridColor(kBlue-1);
seba->SetFrameBorderMode(0);seba->SetFrameFillColor(0);//seba->SetCanvasFillColor(0);
//seba->SetGridStyle(7);
seba->SetLabelSize(0.05,"y");seba->SetTitleXSize(0.05);
seba->SetTitleYSize(0.05);
seba->SetOptStat(kFALSE);gROOT->SetStyle("seba");
std::cout<<"Seba's Style Implemented."<< std::endl;

const Int_t NRGBs = 5;const Int_t NCont = 255;
Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
gStyle->SetNumberContours(NCont);

  TString Metal        = "C";    //"C";
  TString Liq_or_Solid = "liq";    //"liq", "sol";
  TString VC           = "RD";    //"H" , "RD", "KPA"
  TString dataLoc, file;
  dataLoc = "/work/smoran/data/";
  //file = "C.root";
  file = Form( "%s.root" , (const char*)Metal);
  Double_t  delta_q2 , delta_xb , *v_q2 , *v_xb;
  
  const Double_t q2_min  = 1.; 
  const Double_t  q2_max = 4.; 
  const Double_t xb_min = 0.12; 
  const Double_t xb_max = 0.56;
  const Int_t n_q2 = 10 ;  
  const Int_t n_xb = 10 ;
  
  
  delta_q2  = (q2_max-q2_min)/n_q2;   delta_xb  = (xb_max-xb_min)/n_xb;
  v_q2  = new Double_t[n_q2+1];  v_xb  = new Double_t[n_xb+1];
  std::ofstream ofs (Form("centroids_%s_%s_%s_vc.txt", (const char*)Metal , (const char*)Liq_or_Solid, (const char*)VC), std::ofstream::out);
  std::ofstream ofs2(Form("centroids_%s_%s_%s_vc_for_plot.txt",(const char*)Metal,(const char*)Liq_or_Solid,(const char*)VC),std::ofstream::out);
  for(Int_t i = 0; i < n_q2+1; i++){  if(i == 0) v_q2[i] = q2_min;  else v_q2[i] = v_q2[i-1] + delta_q2; }
  for(Int_t i = 0; i < n_xb+1; i++){  if(i == 0) v_xb[i] = xb_min;  else v_xb[i] = v_xb[i-1] + delta_xb; }
  TCut cut_1, cut_2 , liq , sol ,  dis_cut, cut_w, cut_q2 , cut_Yb, target_cut;
  cut_w="W>2"; cut_q2="Q2>1"; cut_Yb="Nu/5.015<0.85"; dis_cut= cut_w && cut_q2 && cut_Yb;
  if (VC=="H"){liq = "TargType==1"; sol = "TargType==2";}
  else if (VC=="KPA"){ liq = "Vert_cut_KPA==1"  ; sol = "Vert_cut_KPA==2";}
  else if (VC=="OS"){ liq = "Vert_cut_OS==1"  ; sol = "Vert_cut_OS==2";}
  else if (VC=="TM"){ liq = "Vert_cut_TM==1"  ; sol = "Vert_cut_TM==2";}
  else if (VC=="RD"){ liq = "VC_RD==1"  ; sol = "VC_RD==2";}
  else if (VC=="RD_f"){ liq = "Vert_cut_RD_f==1"  ; sol = "Vert_cut_RD_f==2";}
  else if (VC=="OS_f"){ liq = "Vert_cut_OS_f==1"  ; sol = "Vert_cut_OS_f==2";}
  else if (VC=="TM_f"){ liq = "Vert_cut_TM_f==1"  ; sol = "Vert_cut_TM_f==2";}


  if (Liq_or_Solid == "liq" ){ target_cut = liq; } else if (Liq_or_Solid == "sol") {target_cut = sol;}
  TFile   *f , *fout; ;TH2F* h;
  Double_t meanx , meany;
  fout = new TFile(  Form("out_%s_%s_%s_vc.root", (const char*)Metal ,(const char*)Liq_or_Solid , (const char*)VC)  ,  "recreate");
  f    = new TFile(dataLoc + file, "READ");
  TChain *t_elec = new TChain();t_elec->Add(Form(dataLoc + file + "/ntuple_data_electrons" ,(const char*)Metal));
  t_elec->SetMaxEntryLoop(20000000);
  t_elec->Draw(">>list_dis","", "goff"); t_elec->SetEventList((TEventList*) gDirectory->Get("list_dis"));
  ofs<<"('Xb','Q2') (";
  
  for (int i = 0 ; i < n_q2 ; i++){
    cut_1 = Form("Q2>%f && Q2<%f", v_q2[i], v_q2[i+1]);
    for (int j = 0 ; j < n_xb  ; j++){
    
      cut_2 = Form("Xb>%f && Xb<%f", v_xb[j], v_xb[j+1]); 
      TCut cut;
      cut = cut_1 && cut_2 && target_cut ;
      t_elec->Draw("Q2:Xb>>h(5,0.0,0.6,5,0.8,4.2)" , cut , "goff"); 
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
      h->Delete();
    } 
  }
  ofs<<")";
  ofs.close(); 
  ofs2.close();
  TCanvas *c = new TCanvas("c", "canvas");c->SetFillColor(0);
  t_elec->Draw("Q2:Xb>>h2(250,0.1,0.6,250,0.8,4.2)",  target_cut , "goff");
  TH2F *h2d = (TH2F*) gDirectory->GetList()->FindObject("h2");
  h2d->GetYaxis()->SetDecimals(kTRUE);
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
  lin_H->SetLineColor(kGray);
  lin_V->SetLineColor(kGray);
  for (int i =0 ; i<n_q2+1 ; i++){ll = v_q2[i];lin_H->DrawLine(0.1, ll,0.59, ll);}
  for (int i =0 ; i<n_xb+1 ; i++){ii = v_xb[i];lin_V->DrawLine(ii, 0.9, ii,4.1) ;}
  TGraph *g = new TGraph(Form("centroids_%s_%s_%s_vc_for_plot.txt",(const char*)Metal,(const char*)Liq_or_Solid,(const char*)VC)  ); 
  g->SetMarkerStyle(kFullDotLarge); 
  g->SetMarkerSize(0.9);g->SetMarkerColor(kRed);
  g->Draw("psame");
  fout->cd();
  c->Write("centroid_plot", TObject::kOverwrite);
  c->SaveAs("out.pdf");
  delete c; 
  delete h2d;
  f->Close(); 
  fout->Close(); 
  delete fout ; 
  delete  f; 
  ofs.close();
  delete v_q2; 
  delete v_xb; 
}

