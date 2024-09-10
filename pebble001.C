#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TSystem.h>
#include "/home/motoko/TKI_Comparison/style/style.h"

// Global pointers to the random number generator and histogram
TRandom3 gRand(123);
TH1D *gHist = nullptr;
TCanvas *gCanvas = nullptr;  // Global canvas pointer

// Generate Gaussian values and fill histogram
void GenerateGaussianValues(TH1D *hist, const double mean, const double stdDev) {
  
    if (!hist) {
      std::cerr << "Error: Histogram not properly initialised!" << std::endl;
      return;
    }
    
    for (int ii = 0; ii < 1e5; ii++) {
        hist->Fill(gRand.Gaus(mean, stdDev));
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
        double fit = par[0] * TMath::Gaus(xx, par[1], par[2]);

        if (data > 0) {
	  chi2 += TMath::Power((fit - data), 2) / TMath::Power(gHist->GetBinError(ii), 2);
        }
    }
    ff = chi2;
}

// Perform the fit using Minuit
int SingleFit(TF1* fitFunc, double* params, double* errors, double& chi2) {
  
    TMinuit Minuit(3);
    Minuit.SetFCN(fcn);

    Minuit.DefineParameter(0, "Amplitude", params[0], 0.1, 0, 1e4);
    Minuit.DefineParameter(1, "Mean", params[1], 0.1, -10, 10);
    Minuit.DefineParameter(2, "Sigma", params[2], 0.1, 0, 10);

    Minuit.Migrad();  // Perform the minimization

    // Retrieve fitted parameters
    for (int ii = 0; ii < 3; ii++) {
        Minuit.GetParameter(ii, params[ii], errors[ii]);
    }

    // Update the fit function parameters
    fitFunc->SetParameters(params[0], params[1], params[2]);

    // Retrieve fitting status
    double amin, edm, errdef;
    int nvpar, nparx, icstat;
    Minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    chi2 = amin;

    return Minuit.GetStatus(); // Return fit status;
}

// Plot the Gaussian fit at each iteration and save it
void PlotGaussianFit(const int iteration, TF1 *fitFunc, const double *params, const double *errors, TFile *outputFile) {
    // Use the global canvas
    if (!gCanvas) {
      gCanvas = new TCanvas(Form("c1_iter_%d", iteration), "Gaussian Fit", 800, 600);
    }
    else {
        gCanvas->cd();  // Make sure to use the existing canvas
    }

    gCanvas->Clear();  // Clear the canvas for new plot

    gHist->SetTitle("Gaussian Fitting");

    // Draw histogram and fit function
    style::ResetStyle(gHist);
    gHist->SetLineColor(kBlue);
    gHist->SetLineWidth(3);
    gHist->SetFillStyle(0);
    gHist->Draw();
    
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(3);
    fitFunc->SetNpx(1e4);
    fitFunc->Draw("same");

    // Legend
    const double x0 = 0.2;
    const double y0 = 0.6;
    const double x1 = x0 + 0.34;
    TLegend *lg = new TLegend(x0, y0, x1, 0.92);
    style::ResetStyle(lg, 0.18, 0.68);
    lg->SetTextAlign(12);

    lg->AddEntry(gHist, "Generated Gaussian", "l");
    lg->AddEntry(fitFunc, "Fitted Gaussian", "l");

    lg->Draw();
    
    // Axis labels
    gHist->GetXaxis()->SetTitle("x-axis");
    gHist->GetYaxis()->SetTitle("y-axis");

    // Print out the plots
    gCanvas->Print(Form("GaussianFit_Iteration_%d.png", iteration));

    // Save the canvas to the ROOT file
    outputFile->cd();
    gCanvas->Write(Form("GaussianFit_Iteration_%d", iteration), TObject::kOverwrite);
    
    delete lg; 
}

// Perform iterative fitting with convergence check
void IterativeFit(TF1* fitFunc, double* params, double* errors, const int maxIterations, const double convergenceThreshold, TFile *outputFile) {
    double previousChi2 = 1e10;
    double currentChi2 = 0.0;

    for (int ii = 0; ii < maxIterations; ii++) {
        // Perform a single fit
        double status = SingleFit(fitFunc, params, errors, currentChi2);

        // Plot and save the results for this iteration
        PlotGaussianFit(ii + 1, fitFunc, params, errors, outputFile);

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
    Int_t status = gSystem->RedirectOutput("see.log", "w"); // "w" to overwrite the file
    if (status != 0) {
        std::cerr << "Error: Unable to redirect output to file." << std::endl;
        return 1;
    }
    
    // Create histogram
    gHist = new TH1D("hist", "Generated Gaussian Data", 100, -10, 10);

    // Generate Gaussian data
    GenerateGaussianValues(gHist, 0.0, 1.0);

    // Initial parameter guesses (Amplitude, Mean, Sigma)
    double params[3] = {
        gHist->GetMaximum(),  // Amplitude
        gHist->GetMean(),     // Mean
        gHist->GetStdDev()    // Sigma
    };
    double errors[3];

    // Create a ROOT file to store results
    TFile *ff = new TFile("GaussianFits.root", "RECREATE");
    if (!ff || ff->IsZombie()) {
        std::cerr << "Error: Cannot create ROOT file." << std::endl;
        return 1;
    }

    // Canvas Setup
    gStyle->SetOptStat(0);
    style::SetGlobalStyle();
    style::fgkYTitleOffset = 1.6;
    style::IniColorCB();

    Double_t currentLeft = 0.17; 
    Double_t currentTop = 0.06; 
    Double_t currentRight = 0.035;
    Double_t currentBottom = 0.14;

    // Ensure the global canvas is created with appropriate margins
    gCanvas = new TCanvas("gCanvas", "Gaussian Fit", 800, 600);
    gPad->SetLeftMargin(currentLeft);
    gPad->SetRightMargin(currentRight);
    gPad->SetTopMargin(currentTop);
    gPad->SetBottomMargin(currentBottom);

    // Title
    gStyle->SetOptTitle(1);
    gStyle->SetTitleX(0.6);
    gStyle->SetTitleW(0.8);

    // Perform the fit
    TF1 *fitFunc = new TF1("fitFunc", "[0] * TMath::Gaus(x,[1],[2])", -1e3, 1e3);
    
    IterativeFit(fitFunc, params, errors, 100, 1e-6, ff);

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
