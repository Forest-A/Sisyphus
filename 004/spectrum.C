#include <cmath>
#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TFile.h>
#include <THStack.h>
#include <Math/Integrator.h>
#include <Math/IntegratorOptions.h>
#include <TLegend.h>

ROOT::Math::IntegratorOneDim *gIntegrator = nullptr;

double noscil(const double xx) {
    // Ensure xx is not zero to avoid division by zero
    if (xx <= 0) {
        std::cerr << "Warning: xx must be greater than 0. Received: " << xx << std::endl;
        return 0; // Return zero or handle the error as needed
    }

    // Electron antineutrino flux from the reactor
    double flux = 0.58 * exp(0.870 - 0.160 * xx - 0.091 * xx * xx)
                + 0.30 * exp(0.896 - 0.239 * xx - 0.0981 * xx * xx)
                + 0.07 * exp(0.976 - 0.162 * xx - 0.0790 * xx * xx)
                + 0.05 * exp(0.793 - 0.080 * xx - 0.1085 * xx * xx);

    // Positron energy (MeV)
    double Ee = xx - (939.565 - 938.272);
    // Positron momentum (MeV)
    double Pe = sqrt(Ee * Ee - 0.511 * 0.511);
    // Inverse beta-decay cross section (cm^2)
    double crossSec = 0.0952e-42 * Ee * Pe;

    return flux * crossSec; // Return non-oscillating result
}

double oscil(const double xx) {
    // Ensure xx is not zero to avoid division by zero
    if (xx <= 0) {
        std::cerr << "Warning: xx must be greater than 0. Received: " << xx << std::endl;
        return 0; // Return zero or handle the error as needed
    }

    // Electron antineutrino flux from the reactor
    double flux = 0.58 * exp(0.870 - 0.160 * xx - 0.091 * xx * xx)
                + 0.30 * exp(0.896 - 0.239 * xx - 0.0981 * xx * xx)
                + 0.07 * exp(0.976 - 0.162 * xx - 0.0790 * xx * xx)
                + 0.05 * exp(0.793 - 0.080 * xx - 0.1085 * xx * xx);

    // Positron energy (MeV)
    double Ee = xx - (939.565 - 938.272);
    // Positron momentum (MeV)
    double Pe = sqrt(Ee * Ee - 0.511 * 0.511);
    // Inverse beta-decay cross section (cm^2)
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
        gIntegrator->SetRelTolerance(1e-4);
    }
}

void plotSpectrum() {
    // Initialize the integrator
    InitializeIntegrator();

    // Create histograms for the non-oscillating and oscillating spectra
    TH1F *hist1 = new TH1F("hist1", "Non-Oscillating Reactor Neutrino Spectrum", 200, 0, 10); // MeV
    TH1F *hist2 = new TH1F("hist2", "Oscillating Reactor Neutrino Spectrum", 200, 0, 10); // MeV

    // Normalisation factor
    double norm = 0;

    // Fill the histograms
    for (int ii = 1; ii <= hist1->GetNbinsX(); ii++) {
        double xx = hist1->GetXaxis()->GetBinCenter(ii);
        if (xx > 1.8) {
            hist1->SetBinContent(ii, noscil(xx));
            hist2->SetBinContent(ii, oscil(xx));
        } else {
            hist1->SetBinContent(ii, 0);
            hist2->SetBinContent(ii, 0);
        }
    }

    // Calculate the integrals for normalization
    norm = hist1->Integral();

    // Normalize the histograms
    if (norm > 0){
     hist1->Scale(1.0 / norm);
     hist2->Scale(1.0 / norm);
    }

    // Create a canvas and draw the histograms
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    hist1->SetLineColor(kBlue);
    hist1->SetLineWidth(2);
    hist1->Draw("HIST");
    
    hist2->SetLineColor(kRed);
    hist2->SetLineWidth(2);
    hist2->Draw("HIST SAME");

    // Add a legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hist1, "Non-Oscillating", "f");
    legend->AddEntry(hist2, "NO Oscillating", "f");
    legend->Draw();


    hist1->SetTitle("Normalized Reactor Neutrino Spectrum; Energy (MeV); Normalized Count");
    c1->SaveAs("spectrum.png");

    delete c1; 
}

int main() {
    plotSpectrum();
    return 0;
}











