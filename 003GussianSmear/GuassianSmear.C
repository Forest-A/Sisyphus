#include <iostream>
#include <cmath>
#include <TRandom3.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include "style.h"

// Global objects
TRandom3 gRand(123);  // Global random number generator
TH1D *gHist = nullptr;  // Global histogram
TCanvas *gCanvas = nullptr;  // Global canvas

// Function to generate noisy sine values and fill the histogram
void SmearSignal(const double params[5]) {
    if (gHist) delete gHist;  // Clean up previous histogram if it exists

    // Parameters for sine function and Gaussian noise
    const double amp = params[0];
    const double freq = params[1];
    const double phase = params[2];
    const double offset = params[3];
    const double noiseSigma = params[4];

    // Create histogram with appropriate range and binning
    gHist = new TH1D("gHist", "Generated Noisy Sine Signal", 200, 0, 100);

    const double nEvents = 1e5;
    const double xMin = 0;
    const double xMax = 100;  // Ensure this matches your range of interest
    const double maxY = amp + offset + 3 * noiseSigma; // Maximum possible y-value used for rejection sampling

    int fillEvents = 0;  // Counter for filled events

    while (fillEvents < nEvents) {
        // Generate a random x-value in the range [xMin, xMax]
        const double xVal = gRand.Uniform(xMin, xMax);

        // Compute the sine value at this x-value
        const double yVal = amp * TMath::Sin(freq * xVal + phase) + offset;

        // Add Gaussian noise to the y-value
        const double noisyYVal = yVal + gRand.Gaus(0, noiseSigma);

        // Generate a random number for comparison
        const double rr = gRand.Uniform(0, maxY); 

        // Acceptance criterion: accept the value if rr < noisyYVal
        if (rr < noisyYVal) {
            gHist->Fill(xVal);
            fillEvents++;  
        }
    }
}

// Function to plot and save the histogram
void Plot() {
    // Ensure canvas exists
    if (!gCanvas) {
        gCanvas = new TCanvas("gCanvas", "Smeared Sine Signal", 800, 600);
    }
    gCanvas->cd();
    gCanvas->Clear();  // Clear the canvas for a new plot

    gHist->SetTitle("Smeared Sine Signal");

    // Draw histogram
    style::ResetStyle(gHist);
    gHist->SetLineColor(kBlue);
    gHist->SetLineWidth(2);
    gHist->SetFillStyle(0);
    gHist->Draw();

    // Legend
    const double x0 = 0.2;
    const double y0 = 0.6;
    const double x1 = x0 + 0.34;
    TLegend *lg = new TLegend(x0, y0, x1, 0.92);
    style::ResetStyle(lg, 0.18, 0.68);
    lg->SetTextAlign(12);

    // lg->AddEntry(gHist, "Noisy Sine Signal", "l");

    lg->Draw();

    // Axis labels
    gHist->GetXaxis()->SetTitle("x-axis");
    gHist->GetYaxis()->SetTitle("Frequency");

    // Print out the plots
    gCanvas->Print("outplot/smeared_signal.png");

    delete lg;  // Clean up local object
}

int main() {
    // Canvas Setup
    gStyle->SetOptStat(0);
    style::SetGlobalStyle();
    style::fgkYTitleOffset = 1.6;
    style::IniColorCB();

    const Double_t currentLeft = 0.17;
    const Double_t currentTop = 0.06;
    const Double_t currentRight = 0.035;
    const Double_t currentBottom = 0.14;

    gCanvas = new TCanvas("gCanvas", "Smeared Sine Signal", 800, 600);
    gPad->SetLeftMargin(currentLeft);
    gPad->SetRightMargin(currentRight);
    gPad->SetTopMargin(currentTop);
    gPad->SetBottomMargin(currentBottom);

    // Title settings
    gStyle->SetOptTitle(1);
    gStyle->SetTitleX(0.6);
    gStyle->SetTitleW(0.8);

    // Parameters for sine function and Gaussian noise
    const double params[5] = {3, 1, 0, 15, 1};  // A, B, C, D, sigma

    // Perform smearing of the sine signal with Gaussian noise
    SmearSignal(params);

    // Plot and save the sampling results
    Plot();

    // Cleanup global objects
    delete gHist;
    delete gCanvas;

    return 0;
}
