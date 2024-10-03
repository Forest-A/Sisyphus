#include <cmath>
#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TFile.h>
#include <THStack.h>
#include <Math/Integrator.h>
#include <Math/IntegratorOptions.h>

ROOT::Math::IntegratorOneDim *gIntegrator = nullptr;

double noscil(const Double_t xx){

  // Electron antineutrino flux from the neutrino
  double flux = 0.58 * exp(0.870 - 0.160 * xx - 0.091 * xx * xx)
              + 0.30 * exp(0.896 - 0.239 * xx - 0.0981 * xx * xx)
              + 0.07 * exp(0.976 - 0.162 * xx - 0.0790 * xx * xx)
              + 0.05 * exp(0.793 - 0.080 * xx - 0.1085 * xx * xx);
  
  // Positron energy (MeV)
  double Ee = xx - ( 939.565 - 938.272);

  // Positron momentum (Mev)
  double Pe = sqrt(Ee * Ee - 0.511 * 0.511);

  // Inverse beta-decay cross section (leading-order expression) (cm^2)
  double crossSec = 0.0952e-42 * Ee * Pe;

  return flux * crossSec;

}

// double oscil(Double_t *xx){

//   // Electron antineutrino flux from the neutrino
//   double flux = 0.58 * exp(0.870 - 0.160 * xx - 0.091 * xx * xx)
//               + 0.30 * exp(0.896 - 0.239 * xx - 0.0981 * xx * xx)
//               + 0.07 * exp(0.976 - 0.162 * xx - 0.0790 * xx * xx)
//               + 0.05 * exp(0.793 - 0.080 * xx - 0.1085 * xx * xx);
  
//   // Positron energy (MeV)
//   double Ee = xx - (939.565 - 938.272);
//   // Positron momentum (Mev)
//   double Pe = sqrt(Ee * Ee - 0.511 * 0.511);
//   // Inverse beta-decay cross section (leading-order expression) (cm^2)
//   double crossSec = 0.0952e-42 * Ee * Pe;

//   /////////////////////////////////////////////////
//   //https://pdg.lbl.gov/2024/listings/rpp2024-list-neutrino-mixing.pdf
//   const double theta12 = 0.587; // radian
//   const double theta13 = 0.1485;

//   const double m21 = 7.53e-5; //eV2
//   const double m32 = 2.455e-3; //NO
//   const double m31 = 2.463e-3;

//   const double delta21 = 1.27 * m21 * 52.5e3 / xx; //52.5 km baseline 
//   const double delta32 = 1.27 * m32 * 52.5e3 / xx;
//   const double delta31 = 1.27 * m31 * 52.5e3 / xx;
//   ////////////////////////////////////////////////
  
//   double P21 = pow(cos(theta13), 4)
//              * pow(sin(2 * theta12), 2)
//              * pow(sin(delta21), 2);

//   double P31 = pow(cos(theta12), 2)
//              * pow(sin(2 * theta13), 2)
//              * pow(sin(delta31), 2);

//   double P32 = pow(sin(theta12), 2)
//              * pow(sin(2 * theta13), 2)
//              * pow(sin(delta32), 2);

//   double prob = 1 - P21 - P31 - P32;

//   return flux * crossSec * prob;

// }

double oscil(const double xx) {

    // Ensure xx is not zero to avoid division by zero
    if (xx <= 0) {
        std::cerr << "Warning: xx must be greater than 0. Received: " << xx << std::endl;
        return 0; // Return zero or handle the error as needed
    }
    
    // Electron antineutrino flux from the neutrino
    double flux = 0.58 * exp(0.870 - 0.160 * xx - 0.091 * xx * xx)
                + 0.30 * exp(0.896 - 0.239 * xx - 0.0981 * xx * xx)
                + 0.07 * exp(0.976 - 0.162 * xx - 0.0790 * xx * xx)
                + 0.05 * exp(0.793 - 0.080 * xx - 0.1085 * xx * xx);
  
    // Positron energy (MeV)
    double Ee = xx - (939.565 - 938.272);
    // Positron momentum (MeV)
    double Pe = sqrt(Ee * Ee - 0.511 * 0.511);
    // Inverse beta-decay cross section (leading-order expression) (cm^2)
    double crossSec = 0.0952e-42 * Ee * Pe;

    // Neutrino mixing parameters
    const double theta12 = 0.587; // radian
    const double theta13 = 0.1485;

    const double m21 = 7.53e-5; // eV^2
    const double m32 = 2.455e-3; // Normal ordering
    const double m31 = 2.463e-3;

    const double delta21 = 1.27 * m21 * 52.5e3 / xx; // 52.5 km baseline 
    const double delta32 = 1.27 * m32 * 52.5e3 / xx;
    const double delta31 = 1.27 * m31 * 52.5e3 / xx;

    // Probability calculations
    double P21 = pow(cos(theta13), 4) * pow(sin(2 * theta12), 2) * pow(sin(delta21), 2);
    double P31 = pow(cos(theta12), 2) * pow(sin(2 * theta13), 2) * pow(sin(delta31), 2);
    double P32 = pow(sin(theta12), 2) * pow(sin(2 * theta13), 2) * pow(sin(delta32), 2);

    double prob = 1 - P21 - P31 - P32;

    return flux * crossSec * prob;
}

// Initialize the integrator and set the function once
void InitializeIntegrator() {
    if (!gIntegrator) {
        gIntegrator = new ROOT::Math::IntegratorOneDim(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
        if (gIntegrator) {
            gIntegrator->SetRelTolerance(1e-4);
            gIntegrator->SetFunction(oscil);  // Set the function only once
        }
        else {
            std::cerr << "Error: Failed to initialize the integrator!" << std::endl;
        }
    }
}

void plotSpectrum(){

  double EE1 = 0;
  double EE2 = 0;

  if (!gIntegrator) {
        InitializeIntegrator();
    }
  
  // double integral = gIntegrator->Integral(1.8, 10);
  // std::cout<< integral << std::endl;

  auto hs = new THStack("hs","Reactor Neutrino Spectrum");
  TH1F *hist1 = new TH1F("hist", "Non-Oscillating Reactor Neutrino Spectrum", 200, 0, 10); // MeV
  TH1F *hist2 = new TH1F("hist", "Non-Oscillating Reactor Neutrino Spectrum", 200, 0, 10); // MeV

  for(int ii = 0; ii < 200; ii++){

    double xx = hist1->GetXaxis()->GetBinCenter(ii);
    if(xx > 1.8){
      EE1 = noscil(xx);
      EE2 = gIntegrator->Integral(1.8, 10);
    }
    else{
      EE1 = 0;
      EE2 = 0;
    }

    // hist->Fill(xx, EE);
    hist1->SetBinContent(ii, EE1);
    hist2->SetBinContent(ii, EE2);
  
  }

  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
  hs->Add(hist1);
  hs->Add(hist2);
  hs->Draw();
  c1->SaveAs("spectrum.png");
  
}


int main(){

  plotSpectrum();
  
  return 0;
}










