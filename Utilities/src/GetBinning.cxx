/*****************************************/
/*  GetBinning.cxx                       */
/*                                       */
/*  Andrés Bórquez                       */
/*                                       */
/*****************************************/

// This program prints the binning of (Nu, Q2) in a CSV file
// June 2021

#include "Binning.hxx"
#include "Headers.hxx"
#include "UX.hxx"

int main() {

  gProgram = "GetBinning";

  /*** WRITE OUTPUT FILE ***/

  // Arrays kEdgesNu and kEdgesQ2 are defined in include/Headers.hxx

  TString OutFileName = "binning.csv";
  std::ofstream OutFile(OutFileName, std::ios::out);

  Int_t NbinsNu = (Int_t)(sizeof(kEdgesNu) / sizeof(kEdgesNu[0]));
  Int_t NbinsQ2 = (Int_t)(sizeof(kEdgesQ2) / sizeof(kEdgesQ2[0]));

  std::cout << "Binning : (Nu, Q2) = " << NbinsNu - 1 << " x " << NbinsQ2 - 1 << std::endl;

  for (Int_t n = 0; n < NbinsNu - 1; n++) {
    for (Int_t q = 0; q < NbinsQ2 - 1; q++) {
      OutFile << kEdgesNu[n] << "," << kEdgesNu[n + 1] << "," << kEdgesQ2[q] << "," << kEdgesQ2[q + 1] << std::endl;
      // std::cout << kEdgesNu[n] << "," << kEdgesNu[n + 1] << "," << kEdgesQ2[q] << "," << kEdgesQ2[q + 1] << std::endl;
    }
  }

  OutFile.close();
  std::cout << "The following file has been created: " << OutFileName << std::endl;
}