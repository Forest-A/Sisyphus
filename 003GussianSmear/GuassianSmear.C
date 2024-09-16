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
TList *plotList = nullptr; 

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

// Chi-squared minimization function
void fcn(int &npar, double *gin, double &ff, double *par, int iflag) {
    if (!gHist) {
        std::cerr << "Error: Histogram not found!" << std::endl;
        return;
    }

    const int nBins = gHist->GetNbinsX();
    double chi2 = 0.0;

    for (int ii = 1; ii <= nBins; ii++) {
        const double xx = gHist->GetBinCenter(ii);
        const double data = gHist->GetBinContent(ii);
        const double fitValue = par[0] * sin(par[1] * xx + par[2]) + par[3];

        if (data > 0) {
            chi2 += TMath::Power((fitValue - data), 2) / TMath::Power(gHist->GetBinError(ii), 2);
        }
    }

    ff = chi2;
}

// Perform the fit using Minuit
int SingleFit(TF1* const fitFunc, double* const params, double* const errors, double& chi2) {
    TMinuit Minuit(4);
    Minuit.SetFCN(fcn);

    Minuit.DefineParameter(0, "Amplitude", params[0], 0.1, 0.1, 1e4);
    Minuit.DefineParameter(1, "Frequency", params[1], 0.01, 0.5, 5.0);
    Minuit.DefineParameter(2, "Phase", params[2], 0.1, -5, 5);
    Minuit.DefineParameter(3, "Offset", params[3], 0.1, -10, 1e4);

    Minuit.Migrad();  // Perform the minimization

    // Retrieve fitted parameters
    for (int ii = 0; ii < 4; ii++) {
        Minuit.GetParameter(ii, params[ii], errors[ii]);
    }

    // Update the fit function parameters
    fitFunc->SetParameters(params[0], params[1], params[2], params[3]);

    // Retrieve fitting status
    double amin, edm, errdef;
    int nvpar, nparx, icstat;
    Minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    chi2 = amin;

    return Minuit.GetStatus();  // Return fit status
}

// Plot the fit at each iteration and store it in the TList
void PlotFit(const int iteration, TF1 * const fitFunc, const double *params, const double *errors) {
    if (!gCanvas) {
        gCanvas = new TCanvas(Form("c1_iter_%d", iteration), "Smear Fit", 800, 600);
    } else {
        gCanvas->cd();  // Make sure to use the existing canvas
    }

    gCanvas->Clear();  // Clear the canvas for new plot

    gHist->SetTitle("SmearSine Fitting");

    // Draw histogram and fit function
    style::ResetStyle(gHist);
    gHist->SetLineColor(kBlue);
    gHist->SetLineWidth(3);
    gHist->SetFillStyle(0);
    gHist->Draw();
    
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(3);
    fitFunc->SetNpx(1e5);
    fitFunc->Draw("same");

    // Legend
    const double x0 = 0.2;
    const double y0 = 0.6;
    const double x1 = x0 + 0.34;
    TLegend *lg = new TLegend(x0, y0, x1, 0.92);
    style::ResetStyle(lg, 0.18, 0.68);
    lg->SetTextAlign(12);

    lg->AddEntry(gHist, "Generated Signal", "l");
    lg->AddEntry(fitFunc, "Fitted Line", "l");

    lg->Draw();
    
    // Axis labels
    gHist->GetXaxis()->SetTitle("x-axis");
    gHist->GetYaxis()->SetTitle("y-axis");

    // Optimise scales
    gHist->GetYaxis()->SetRangeUser(200, 1000);
    
    // Print out the plots
    gCanvas->Print(Form("outplot/SmearFit_%d.png", iteration));

    delete lg; 
}

// Perform iterative fitting with convergence check
void IterativeFit(TF1* const fitFunc, double* const params, double* const errors, const int maxIterations, const double convergenceThreshold) {
    double previousChi2 = 1e10;
    double currentChi2 = 0.0;

    for (int ii = 0; ii < maxIterations; ii++) {
        // Perform a single fit
        const int status = SingleFit(fitFunc, params, errors, currentChi2);

	// Output the chi-squared value for this iteration
        std::cout << "////////////////////////////Iteration " << ii + 1 << ": Chi2 = " << currentChi2 << "////////////////////////////" << std::endl;

        // Plot and store the results for this iteration
        PlotFit(ii + 1, fitFunc, params, errors);

        // Check for convergence
        if (status != 0) {
            std::cerr << "Fit did not converge. Status: " << status << std::endl;
            break;
        }

        // Check if the change in chi-squared is below the threshold
        if (TMath::Abs(previousChi2 - currentChi2) < convergenceThreshold) {
            std::cout << "Converged after " << ii + 1 << " iterations." << std::endl;
            break;
        }

        previousChi2 = currentChi2;
    }
}

int main() {

     // Redirect ROOT output to the log file
    const Int_t status = gSystem->RedirectOutput("outplot/see.log", "w");  // "w" to overwrite the file
    if (status != 0) {
        std::cerr << "Error: Unable to redirect output to file." << std::endl;
        return 1;
    }

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
    const double generate[5] = {3, 1, 0, 15, 1};  // A, B, C, D, sigma

    // Perform smearing of the sine signal with Gaussian noise
    SmearSignal(generate);

    // Initial parameter guesses (Amplitude, Frequency, Phase, Offset)
    double params[4] = {1.0, 1.0, 0.0, 2.0}; 
    double errors[4];

    // Create a list to store plots
    plotList = new TList();
    plotList->Add(gHist);
    
    // Create a ROOT file to store results
    TFile *ff = new TFile("outplot/SmearFits.root", "RECREATE");
    if (!ff || ff->IsZombie()) {
        std::cerr << "Error: Cannot create ROOT file." << std::endl;
        return 1;
    }

    // Perform the fit
    TF1 *fitFunc = new TF1("fitFunc", "[0] * sin([1] * x + [2]) + [3]", 0, 100);
    
    IterativeFit(fitFunc, params, errors, 10, 1e-6);

    // Write the list of plots to the ROOT file
    plotList->Write();

    // Save everything and close the file
    ff->Write(); 
    ff->Close();

    // Restore output to the default streams
    gSystem->RedirectOutput(0);

    // Cleanup
    delete ff;
    delete fitFunc;
    delete gHist;
    delete gCanvas;

    return 0;
}
