#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TRandom3.h>
#include <TStyle.h>

// Global pointer to the random number generator
TRandom3 *gRand = new TRandom3(123);

// Generate Gaussian values and fill histogram
void GenerateGaussianValues(TH1D *hist, double mean, double stdDev) {
	for (int ii = 0; ii < 1e5; ii++) {
		 hist->Fill(gRand->Gaus(mean, stdDev));
	}
}

// Chi-squared minimization function
void fcn(int &npar, double *gin, double &ff, double *par, int iflag) {
    TH1D *hist = static_cast<TH1D*>(gDirectory->Get("hist"));
    if (!hist) {
        std::cerr << "Error: Histogram not found!" << std::endl;
        return;
    }

    int nBins = hist->GetNbinsX();
    double chi2 = 0.0;

    for (int ii = 1; ii <= nBins; ii++) {
        double x = hist->GetBinCenter(ii);
        double data = hist->GetBinContent(ii);
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

// Plot the Gaussian fit
void PlotGaussianFit(TH1D *hist, TF1 *fitFunc) {
    TCanvas *c1 = new TCanvas("c1", "Gaussian Fit", 800, 600);
    gStyle->SetOptStat(0);

    hist->Draw();
    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("same");

    c1->SaveAs("gaussian_fit.png");
    delete c1;
}

int main() {
    // Create histogram
    TH1D *hist = new TH1D("hist", "Generated Gaussian Data", 100, -10, 10);

    // Generate Gaussian data
    GenerateGaussianValues(hist, 0.0, 1.0);

    // Initial parameter guesses (Amplitude, Mean, Sigma)
    double params[3] = {
        hist->GetMaximum(),  // Amplitude
        hist->GetMean(),     // Mean
        hist->GetStdDev()    // Sigma
    };
    double errors[3];

    // Perform the fit
    TF1 *fitFunc = new TF1("fitFunc", "[0]*TMath::Gaus(x,[1],[2])", -1e3, 1e3);
    PerformSingleFit(fitFunc, params, errors);

    // Plot the results
    PlotGaussianFit(hist, fitFunc);

    // Output results

    std::cout << "Fitted parameters:" << std::endl;
    std::cout << "Amplitude: " << params[0] << " ± " << errors[0] << std::endl;
    std::cout << "Mean: " << params[1] << " ± " << errors[1] << std::endl;
    std::cout << "Sigma: " << params[2] << " ± " << errors[2] << std::endl;

    delete fitFunc;
    delete hist;
    return 0;
}
