#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TRandom3.h>
#include <TStyle.h>
#include "/home/motoko/TKI_Comparison/style/style.h"

// Global pointer to the random number generator and the histogram
TRandom3 *gRand = new TRandom3(123);
TH1D *gHist = nullptr;

// Generate Gaussian values and fill histogram
void GenerateGaussianValues(TH1D *hist, const double mean, const double stdDev) {
	for (int ii = 0; ii < 1e5; ii++) {
		 hist->Fill(gRand->Gaus(mean, stdDev));
	}
}

// Chi-squared minimization function
void fcn(int &npar, double *gin, double &ff, double *par, int iflag) {
  
    if (!gHist) {
        std::cerr << "Error: Histogram not found!" << std::endl;
        return;
    }

    int nBins = gHist->GetNbinsX();
    double chi2 = 0.0;

    for (int ii = 1; ii <= nBins; ii++) {
        double x = gHist->GetBinCenter(ii);
        double data = gHist->GetBinContent(ii);
        double fit = par[0] * TMath::Gaus(x, par[1], par[2]);

        if (data > 0) {
            chi2 += TMath::Power((fit - data), 2) / data;
        }
    }
    ff = chi2;
}

// Perform the fit using Minuit
void PerformSingleFit(TF1* fitFunc, double* params, double* errors) {
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
}

// Plot the Gaussian fit and save it
void PlotGaussianFit(TH1D *hist, TF1 *fitFunc, const double *params,  const double *errors ) {
    TCanvas *c1 = new TCanvas("c1", "Gaussian Fit", 800, 600);
    gStyle->SetOptStat(0);

    // Canvas Setup
    style::SetGlobalStyle();
    style::fgkYTitleOffset=1.6;//1.2;//square
    style::IniColorCB();

    Double_t currentLeft=0.17; 
    Double_t currentTop=0.06; 
    Double_t currentRight=0.035;
    Double_t currentBottom=0.14;

    style::PadSetup(c1);
    //cc->GetCanvas()->SetGrayscale();

    gPad->SetLeftMargin(currentLeft);//square
    gPad->SetRightMargin(currentRight);//0.03//square
    gPad->SetTopMargin(currentTop);
    gPad->SetBottomMargin(currentBottom);

    // Title
    gStyle->SetOptTitle(1);
    gStyle->SetTitleX(0.6);
    gStyle->SetTitleW(0.8);
    hist->SetTitle("Gaussian Fitting");

    // Draw lines 
    style::ResetStyle(hist);
    hist->SetLineColor(kBlue);
    hist->SetLineWidth(3);
    hist->SetFillStyle(0);
    hist->Draw();
    
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(3);
    fitFunc->SetNpx(1e4);
    fitFunc->Draw("same");

    // Legend
    const double x0 = 0.2;
    const double y0 = 0.6;
    const double x1 = x0 + 0.34;
    TLegend *lg = new TLegend(x0, y0, x1, 0.92);
    style::ResetStyle(lg,0.18,0.68);
    lg->SetTextAlign(12);

    lg->AddEntry(hist, "Generated Gaussian", "l");
    lg->AddEntry(fitFunc, "Fitted Gaussian", "l");

    lg->Draw();
    
    // Axis label
    hist->GetXaxis()->SetTitle("x-axis");
    hist->GetYaxis()->SetTitle("y-axis");

    // Save
    // Create a new ROOT file
    TFile *file = new TFile("GaussianFit.root", "RECREATE");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot create ROOT file." << std::endl;
        return;
    }

    // Create a TList
    TList *list = new TList();

    // Add objects to the TList
    list->Add(hist);
    list->Add(fitFunc);
    list->Add(c1);

    list->Write("aList");

    
    // c1->SaveAs("gaussian_fit.png");
    file->Close();

    delete c1;
    delete list;
    delete file;
}

int main() {
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

    // Perform the fit
    TF1 *fitFunc = new TF1("fitFunc", "[0] * TMath::Gaus(x,[1],[2])", -1e3, 1e3);
    PerformSingleFit(fitFunc, params, errors);

    // Plot and save the results
    PlotGaussianFit(gHist, fitFunc, params, errors);
    

    delete fitFunc;
    delete gHist;
    return 0;
}
