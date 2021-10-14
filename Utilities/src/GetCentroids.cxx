/*****************************************/
/*  GetCentroids.cxx                     */
/*                                       */
/*  Andrés Bórquez                       */
/*                                       */
/*****************************************/

// This program gets the centroids of (Nu, Q2) bins based on a (Nu, Q2) binning given by a CSV file
// October 2021

#include "Binning.hxx"
#include "Headers.hxx"
#include "UX.hxx"

int main(int argc, char **argv) {

  gProgram = "GetCentroids";

  /*** INPUT ***/

  parseCommandLine(argc, argv);
  printOptions();

  TChain *DataChain = new TChain();

  if (gTargetOption == "D") {  // unified D
    DataChain->Add(gWorkDir + "/out/GetSimpleTuple/data/C/*.root/ntuple_e");
    DataChain->Add(gWorkDir + "/out/GetSimpleTuple/data/Fe/*.root/ntuple_e");
    DataChain->Add(gWorkDir + "/out/GetSimpleTuple/data/Pb/*.root/ntuple_e");
  } else {
    DataChain->Add(gWorkDir + "/out/GetSimpleTuple/data/" + gTargetOption + "/*.root/ntuple_e");
  }

  TChain *SimChain = new TChain();
  // assign dirs
  if (gTargetOption == "D") {
    SimChain->Add(gWorkDir + "/out/GetSimpleTuple/" + gParticleOption + "-sim/" + gTargetOption + "/00/*.root/ntuple_e");
    SimChain->Add(gWorkDir + "/out/GetSimpleTuple/" + gParticleOption + "-sim/" + gTargetOption + "/01/*.root/ntuple_e");
    SimChain->Add(gWorkDir + "/out/GetSimpleTuple/" + gParticleOption + "-sim/" + gTargetOption + "/02/*.root/ntuple_e");
  } else if (gTargetOption == "Fe") {
    SimChain->Add(gWorkDir + "/out/GetSimpleTuple/" + gParticleOption + "-sim/" + gTargetOption + "/00/*.root/ntuple_e");
    SimChain->Add(gWorkDir + "/out/GetSimpleTuple/" + gParticleOption + "-sim/" + gTargetOption + "/01/*.root/ntuple_e");
  } else if (gTargetOption == "C" || gTargetOption == "Pb") {
    SimChain->Add(gWorkDir + "/out/GetSimpleTuple/" + gParticleOption + "-sim/" + gTargetOption + "/00/*.root/ntuple_e");
  }

  // define vertex cuts for data
  TCut CutVertex;
  if (gTargetOption == "D") {
    CutVertex = "TargType == 1 && vyec > -1.4 && vyec < 1.4";
  } else {  // in case of solid targets: C, Fe, Pb
    CutVertex = "TargType == 2 && vyec > -1.4 && vyec < 1.4";
  }

  // define cuts for data and sim. rec.
  TCut CutDIS = "Q2 > 1. && W > 2. && Yb < 0.85";

  // define cuts for MC
  TCut CutDIS_MC = "mc_Q2 > 1. && mc_W > 2. && mc_Yb < 0.85";

  /*** READ CSV FILE ***/

  std::ifstream BinningFile;
  BinningFile.open("binning_" + gParticleOption + ".csv", std::ios::in);
  std::vector<Double_t> Nu_low;
  std::vector<Double_t> Nu_up;
  std::vector<Double_t> Q2_low;
  std::vector<Double_t> Q2_up;

  // exit program if ifstream could not open file
  if (!BinningFile) {
    std::cerr << "ERROR: File could not be opened" << std::endl;
    return 1;
  }

  TString auxLine;
  TString auxString[4];
  while (BinningFile >> auxLine) {
    for (Int_t i = 0; i < 4; i++) {
      auxString[i] = ((TObjString *)(auxLine.Tokenize(",")->At(i)))->String();
    }
    Nu_low.push_back(auxString[0].Atof());
    Nu_up.push_back(auxString[1].Atof());
    Q2_low.push_back(auxString[2].Atof());
    Q2_up.push_back(auxString[3].Atof());
  }

  /*** GET EDGES & NUMBER OF BINS ***/

  // related to particle option - now, get edges
  Double_t EdgesQ2[NbinsQ2 + 1];
  Double_t EdgesNu[NbinsNu + 1];
  if (gParticleOption == "omega") {
    for (Int_t q = 0; q < NbinsQ2 + 1; q++) {
      EdgesQ2[q] = kEdgesQ2_Omega[q];
    }
    for (Int_t n = 0; n < NbinsNu + 1; n++) {
      EdgesNu[n] = kEdgesNu_Omega[n];
    }
  } else {  // "eta"
    for (Int_t q = 0; q < NbinsQ2 + 1; q++) {
      EdgesQ2[q] = kEdgesQ2_Eta[q];
    }
    for (Int_t n = 0; n < NbinsNu + 1; n++) {
      EdgesNu[n] = kEdgesNu_Eta[n];
    }
  }

  Int_t NTotalBins = (Int_t)Nu_low.size();

  /*** HISTOGRAMS ***/

  // define data, mc, and sim. rec. histogram
  TH2D *Hist_Data[NTotalBins];
  TH2D *Hist_MC[NTotalBins];
  TH2D *Hist_Sim[NTotalBins];

  // define acceptance and corr hists
  TH2D *Hist_Acceptance[NTotalBins];
  TH2D *Hist_Corr[NTotalBins];

  Double_t meanX[NTotalBins];
  Double_t meanY[NTotalBins];

  // define same hist properties for every hist
  TString HistProperties = Form("(200, %.2f, %.2f, 200, %.2f, %.2f)", EdgesNu[0], EdgesNu[NbinsNu], EdgesQ2[0], EdgesQ2[NbinsQ2]);

  TCut CutBin;
  TCut CutBin_MC;

  // loop over listed bins in file
  for (Int_t i = 0; i < NTotalBins; i++) {
    // update cut
    CutBin = Form("Q2 > %.2f && Q2 < %.2f && Nu > %.2f && Nu < %.2f", Q2_low[i], Q2_up[i], Nu_low[i], Nu_up[i]);
    CutBin_MC = Form("mc_Q2 > %.2f && mc_Q2 < %.2f && mc_Nu > %.2f && mc_Nu < %.2f", Q2_low[i], Q2_up[i], Nu_low[i], Nu_up[i]);

    // make data histogram
    DataChain->Draw(Form("Q2:Nu>>data_%i", i) + HistProperties, CutDIS && CutVertex && CutBin, "goff");
    Hist_Data[i] = (TH2D *)gROOT->FindObject(Form("data_%i", i));

    // make sim. rec. histogram
    SimChain->Draw(Form("Q2:Nu>>sim_%i", i) + HistProperties, CutDIS && CutBin, "goff");
    Hist_Sim[i] = (TH2D *)gROOT->FindObject(Form("sim_%i", i));

    // make mc histogram
    SimChain->Draw(Form("mc_Q2:mc_Nu>>mc_%i", i) + HistProperties, CutDIS_MC && CutBin_MC, "goff");
    Hist_MC[i] = (TH2D *)gROOT->FindObject(Form("mc_%i", i));

    // make acceptance
    Hist_Acceptance[i] = new TH2D(Form("acc_%i", i), "", 200, EdgesNu[0], EdgesNu[NbinsNu], 200, EdgesQ2[0], EdgesQ2[NbinsQ2]);
    Hist_Acceptance[i]->Divide(Hist_Sim[i], Hist_MC[i], 1, 1, "B");

    // correct data
    Hist_Corr[i] = new TH2D(Form("corr_%i", i), "", 200, EdgesNu[0], EdgesNu[NbinsNu], 200, EdgesQ2[0], EdgesQ2[NbinsQ2]);
    Hist_Corr[i]->Divide(Hist_Data[i], Hist_Acceptance[i], 1, 1);

    // get 2D centroids
    meanX[i] = Hist_Corr[i]->GetMean(1);
    meanY[i] = Hist_Corr[i]->GetMean(2);
  }

  /*** OUTPUT FILE ***/

  TString OutputFilename = "centroids_" + gParticleOption + "_" + gTargetOption + ".txt";
  std::ofstream OutFile(OutputFilename, std::ios::out);

  OutFile << "('Nu','Q2') (";
  for (Int_t i = 0; i < NTotalBins; i++) {
    OutFile << "(" << meanX[i] << "," << meanY[i] << ")";
    if (i + 1 == NTotalBins)
      OutFile << ")";
    else
      OutFile << ",";
  }
  OutFile.close();

  std::cout << "The following file has been created: " << OutputFilename << std::endl;

  return 0;
}
