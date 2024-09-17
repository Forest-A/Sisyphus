#include <iostream>
#include <cmath>
#include <TRandom3.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include "style.h"

// Global objects
TRandom3 gRand(123);
TH1D *gHist = nullptr;
TCanvas *gCanvas = nullptr;

// Function to generate sine values and Gaussian noise, then combine them into a single histogram
void SmearSignal(TH1D *const hist, const double *params) {
    if (!hist) {
        std::cerr << "Error: Histogram not properly initialized!" << std::endl;
        return;
    }

    // Parameters for sine function and Gaussian noise
    const double amp = params[0];
    const double freq = params[1];
    const double phase = params[2];
    const double offset = params[3];
    const double noiseSigma = params[4];
    
    const double nEvents = 1e5;
    const double xMin = 0;
    const double xMax = 40;  // Adjust the x-range as necessary
    const double smaxY = amp + offset; // Maximum possible y-value used for rejection sampling

    const double mean = 20.0;
    const double gmaxY = TMath::Gaus(mean, mean, noiseSigma, noiseSigma);
    
    int sfillEvents = 0;  // Counter for filled events
    int gfillEvents = 0;  // Counter for filled events

    // Create temporary histograms
    TH1D *histSine = new TH1D("histSine", "Sine Signal", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    TH1D *histNoise = new TH1D("histNoise", "Gaussian Noise", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    // Acceptance-Rejection sampling to generate sine values
    while (sfillEvents < nEvents) {
        const double xVal = gRand.Uniform(xMin, xMax);
        const double yVal = amp * TMath::Sin(freq * xVal + phase) + offset;
        const double rr = gRand.Uniform(0, smaxY);
        if (rr < yVal) {
            histSine->Fill(xVal);
            sfillEvents++;
        }
    }

   // Generate Gaussian noise and fill the noise histogram
     while (gfillEvents < nEvents) {
        // Generate a random x-value in the range [xMin, xMax]
        const double xVal = gRand.Uniform(xMin, xMax);
        const double yVal = TMath::Gaus(xVal, mean, noiseSigma);
        const double rr = gRand.Uniform(0, gmaxY);

        // Acceptance criterion: accept the value if rr < yVal
        if (rr < yVal) {
            histNoise->Fill(xVal);
            gfillEvents++;
        }
    }

    // Combine sine values and Gaussian noise into the final histogram
    for (int ii = 1; ii <= hist->GetNbinsX(); ii++) {
        const double sineContent = histSine->GetBinContent(ii);
        const double noiseContent = histNoise->GetBinContent(ii);
        hist->SetBinContent(ii, sineContent + noiseContent);
    }

    // Plot the noise histogram for debugging
    TCanvas *c1 = new TCanvas("c1", "Noise Histogram", 800, 600);
    histNoise->SetLineColor(kRed); // Set color for visibility
    histNoise->SetLineWidth(2);    // Set line width for visibility
    histNoise->Draw();             // Draw histogram

    // Save the plot
    c1->SaveAs("outplot/histNoise.png");

    // Clean up
    delete c1;
    delete histSine;
    delete histNoise;
}

// Function to plot and save the histogram
void Plot() {
    if (!gCanvas) {
        gCanvas = new TCanvas("gCanvas", "Smeared Sine Signal", 800, 600);
    }
    gCanvas->cd();
    gCanvas->Clear();  // Clear the canvas for a new plot

    gHist->SetTitle("Smeared Sine Signal");

    // Draw histogram
    style::ResetStyle(gHist);
    gHist->SetLineColor(kBlue);
    gHist->SetLineWidth(3);
    gHist->SetFillStyle(0);
    gHist->Draw();

    // Legend
    const double x0 = 0.2;
    const double y0 = 0.6;
    const double x1 = x0 + 0.34;
    TLegend *lg = new TLegend(x0, y0, x1, 0.92);
    style::ResetStyle(lg, 0.18, 0.68);
    lg->SetTextAlign(12);
    lg->AddEntry(gHist, "Noisy Sine Signal", "l");
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

    // Initialize histogram for sine values
    gHist = new TH1D("gHist", "Sine Histogram", 200, 0, 40);

    // Parameters for sine function and Gaussian noise
    const double params[5] = {1, 1, 0, 15, 1};  // A, B, C, D, sigma

    // Perform smearing of the sine signal with Gaussian noise
    SmearSignal(gHist, params);

    // Plot and save the sampling results
    Plot();

    // Cleanup global objects
    delete gHist;
    delete gCanvas;

    return 0;
}
