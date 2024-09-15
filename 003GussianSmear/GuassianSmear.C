#include <iostream>
#include <cmath>
#include <TMath.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include <TH1D.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>
#include "style.h"

// Global objects
TRandom3 gRnd(0);  // Global random number generator
TH1D *gHist = nullptr;  // Global histogram
TCanvas *gCanvas = nullptr;  // Global canvas

// Define the sine function
double Sine(const double xx, const double *params) {
    const double A = params[0];  // Amplitude
    const double B = params[1];  // Frequency
    const double C = params[2];  // Phase shift
    const double D = params[3];  // Offset
    return A * sin(B * xx + C) + D;
}

// Define the Gaussian function
double Gaussian(const double xx, const double sigma) {
    return exp(-0.5 * (xx * xx) / (sigma * sigma));
}

// Convolution of the sine function with the Gaussian
double convolvedFunction(const double xx, const double *params) {
    // Integration range for convolution 
    const double x_min = 0.0;
    const double x_max = 100.0;
    const int steps = 1e3;  // Number of steps for the numerical integration
    const double dx = (x_max - x_min) / steps;

    double result = 0;
    
    // Numerical convolution
    for (int ii = 0; ii < steps; ii++) {
        double x_prime = x_min + ii * dx;
        double sine_val = Sine(x_prime, params);  
        double gauss_val = Gaussian(xx - x_prime, params[4]);  
        result += sine_val * gauss_val * dx;
    }

    return result;
}

// Acceptance-rejection sampling
void ARSampling(const double *params) {
    const double x_min = 0;
    const double x_max = 100;
    const double maxVal = 1.0;  // Assumed maximum value of the convolved function
    int fillEvents = 0;  // Counter for filled events

    // Ensure histogram exists
    if (!gHist) {
        gHist = new TH1D("hist", "Samples from Convolved Function", 200, 0, 100);
    }

    for (int ii = 0; ii < 1e5; ii++) {
        double xx = gRnd.Uniform(x_min, x_max);  // Random uniform number between x_min and x_max
        double yy = gRnd.Uniform(0, maxVal);  // Random uniform number between 0 and maxVal

        // Check if we accept the sample based on the convolved function value
        double f_x = convolvedFunction(xx, params);
        if (yy < f_x) {
            gHist->Fill(xx);  // Accept the sample and fill the histogram
            fillEvents++; 
        }
    }
}

// Function to plot and save the histogram
void Plot() {
    // Ensure canvas exists
    if (!gCanvas) {
        gCanvas = new TCanvas("gCanvas", "Sampled Convolved Function", 800, 600);
    } 
    gCanvas->cd();
    gCanvas->Clear();  // Clear the canvas for new plot

    gHist->SetTitle("Convolution");

    // Draw histogram and fit function
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

    // lg->AddEntry(gHist, "Generated Convolution", "l");
    // lg->AddEntry(fitFunc, "Fitted Sine", "l");

    lg->Draw();
    
    // Axis labels
    gHist->GetXaxis()->SetTitle("x-axis");
    gHist->GetYaxis()->SetTitle("y-axis");
    
    // Print out the plots
    gCanvas->Print("outplot/convolved_sampling.png");

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

    gCanvas = new TCanvas("gCanvas", "Sampled Convolved Function", 800, 600);
    gPad->SetLeftMargin(currentLeft);
    gPad->SetRightMargin(currentRight);
    gPad->SetTopMargin(currentTop);
    gPad->SetBottomMargin(currentBottom);


    // Title settings
    gStyle->SetOptTitle(1);
    gStyle->SetTitleX(0.6);
    gStyle->SetTitleW(0.8);

    // Parameters for sine and Gaussian functions
    const double params[5] = {0.5, 2, 0, 1, 1};  // A, B, C, D, sigma = 1

    // Perform acceptance-rejection sampling
    ARSampling(params);

    // Plot and save the sampling results
    Plot();

    // Cleanup global objects
    delete gHist;
    delete gCanvas;

    return 0;
}
