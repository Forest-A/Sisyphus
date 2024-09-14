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
#include <TList.h>
#include <TVirtualFFT.h>
#include "style.h"

// Global pointers to the random number generator and histogram
TRandom3 gRand(123);
TH1D *gHist = nullptr;
TCanvas *gCanvas = nullptr;
TList *plotList = nullptr;  

// Generate sine values and fill histogram with exactly nEvents entries
void GenerateSineValue(TH1D *const hist, const double amp, const double freq, const double phase, const double offset) {
    if (!hist) {
        std::cerr << "Error: Histogram not properly initialized!" << std::endl;
        return;
    }

    const double nEvents = 1e5;
    const double xMin = 0;
    const double xMax = 40;  // Adjust the x-range as necessary
    const double maxY = amp + offset; // Maximum possible y-value used for rejection sampling
    
    int fillEvents = 0;  // Counter for filled events

    // Acceptance-Rejection sampling to generate exactly nEvents
    while (fillEvents < nEvents) {
        // Generate a random x-value in the range [xMin, xMax]
        const double xVal = gRand.Uniform(xMin, xMax);

        // Compute the sine value at this x-value
        const double  yVal = amp * TMath::Sin(freq * xVal + phase) + offset;

        // Generate a random number for comparison
        const double rr = gRand.Uniform(0, maxY); 

        // Acceptance criterion: accept the value if rr < yVal
        if (rr < yVal) {
            hist->Fill(xVal);
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
void PlotSineFit(const int iteration, TF1 * const fitFunc, const double *params, const double *errors) {
    if (!gCanvas) {
        gCanvas = new TCanvas(Form("c1_iter_%d", iteration), "Sine Fit", 800, 600);
    } else {
        gCanvas->cd();  // Make sure to use the existing canvas
    }

    gCanvas->Clear();  // Clear the canvas for new plot

    gHist->SetTitle("Sine Fitting");

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

    lg->AddEntry(gHist, "Generated Sine", "l");
    lg->AddEntry(fitFunc, "Fitted Sine", "l");

    lg->Draw();
    
    // Axis labels
    gHist->GetXaxis()->SetTitle("x-axis");
    gHist->GetYaxis()->SetTitle("y-axis");

    // Optimise scales
    gHist->GetYaxis()->SetRangeUser(200, 1000);
    
    // Print out the plots
    gCanvas->Print(Form("outplot/SineFit_%d.png", iteration));

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
        PlotSineFit(ii + 1, fitFunc, params, errors);

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

// Fit with ROOT's internal fitting function and print the results
void FitWithROOT(TF1* const fitFunc, double* const params, double* const errors) {
    // Perform the fit using ROOT's internal fitting function
    gHist->Fit(fitFunc, "Q");

    std::cout << "------The internal fitting function results:" << std::endl;

    // Retrieve and print the fit results
    for (int ii = 0; ii < 4; ii++) {
        params[ii] = fitFunc->GetParameter(ii);
        errors[ii] = fitFunc->GetParError(ii);
        std::cout << "Parameter " << ii << ": " << params[ii] << " ± " << errors[ii] << std::endl;
    }

    // Print the chi-square value
    const double chi2 = fitFunc->GetChisquare();
    std::cout << "Chi-square: " << chi2  << std::endl;
}

//Extract data from the histogram directly with FFT
void Extraction(TH1D *const hist) {
    if (!hist) {
        std::cerr << "Error: Histogram not properly initialized!" << std::endl;
        return;
    }

    // Get histogram properties
    const int binmax = hist->GetMaximumBin();
    const int binmin = hist->GetMinimumBin();
    const double maxValue = hist->GetMaximum();
    const double minValue = hist->GetMinimum();
    const double maxError = hist->GetBinError(binmax);
    const double minError = hist->GetBinError(binmin);
    
    // Calculate vertical offset
    const double offset = (maxValue + minValue) / 2.0;
    const double offsetErr = std::sqrt(maxError * maxError + minError * minError);

    // Calculate amplitude
    const double amp = (maxValue - minValue) / 2.0;
    const double ampErr = std::sqrt(maxError * maxError + minError * minError);

    // Number of bins in the histogram
    const Int_t nBins = hist->GetNbinsX();
    const Int_t n_size = nBins;

    // Create an array to hold the size for TVirtualFFT
    Int_t fftSize = n_size;
    Int_t* sizeArray = &fftSize;

    // Create a TVirtualFFT object for real-to-complex FFT
    TVirtualFFT *fft = TVirtualFFT::FFT(1, sizeArray, "R2C ES K");
    if (!fft) {
        std::cerr << "Error: FFT object creation failed!" << std::endl;
        return;
    }

    // Create an array to hold the histogram data
    Double_t *in = new Double_t[n_size];
    for (int ii = 0; ii < n_size; ii++) {
        in[ii] = hist->GetBinContent(ii + 1);
    }

    // Perform transform
    fft->SetPoints(in);
    fft->Transform();

    // Allocate arrays for real and imaginary parts
    Double_t *real = new Double_t[n_size];
    Double_t *imaginary = new Double_t[n_size];
    fft->GetPointsComplex(real, imaginary);

    // Compute the magnitude
    double magnitudes[n_size / 2];
    for (int ii = 0; ii < n_size / 2; ii++) {
        magnitudes[ii] = TMath::Sqrt(real[ii] * real[ii] + imaginary[ii] * imaginary[ii]);
    }

    // Find the bin with the highest magnitude using a loop
    int maxBin = 1;
    double maxMagnitude = magnitudes[1];
    for (int ii = 2; ii < n_size / 2; ii++) {
        if (magnitudes[ii] > maxMagnitude) {
            maxMagnitude = magnitudes[ii];
            maxBin = ii;
        }
    }

    // Calculate the dominant frequency
    const double binWidth = hist->GetXaxis()->GetBinWidth(1);  
    const double dominantFreq = (maxBin * 1.0) / (n_size * binWidth); 

    // Extract the phase for the dominant frequency
    const double phase = TMath::ATan2(imaginary[maxBin], real[maxBin]);
    
    // Frequency error (resolution-based)
    const double freqResolution = 1.0 / (nBins * binWidth);  // Frequency resolution used as error estimation
    
    // Phase error (SNR)
    double noiseFloor = 0.0;
    int noiseBins = 0;
    for (int ii = nBins / 4; ii < nBins / 2; ii++) {  // Estimate noise floor from a region with no dominant peaks
        noiseFloor += magnitudes[ii];
        noiseBins++;
    }
    if (noiseBins > 0) {
        noiseFloor /= noiseBins;  // Average noise floor
    }
    else {
        noiseFloor = 0.0;  // Avoid division by zero
    }

    const double SNR = maxMagnitude / noiseFloor;  // Signal-to-noise ratio
    double phaseError;
    if (SNR > 0) {
        phaseError = 1.0 / SNR;
    } else {
        phaseError = 0.0;
    }

    // Output the results
    std::cout << "------The FFT extraction results:" << std::endl;
    std::cout << "Amplitude: " << amp << " ± " << ampErr << std::endl;
    std::cout << "Dominant Frequency: " << dominantFreq << " ± " << freqResolution << std::endl;
    std::cout << "Phase at Dominant Frequency: " << phase << " ± " << phaseError << std::endl;
    std::cout << "Offset: " << offset << " ± " << offsetErr << std::endl;

    // Clean up
    delete[] in;
    delete[] real;
    delete[] imaginary;
    delete fft;
}




int main() {
    // Redirect ROOT output to the log file
    const Int_t status = gSystem->RedirectOutput("outplot/see.log", "w");  // "w" to overwrite the file
    if (status != 0) {
        std::cerr << "Error: Unable to redirect output to file." << std::endl;
        return 1;
    }

    // Create histogram
    gHist = new TH1D("hist", "Generated Sine Data", 200, 0, 40);

    // Generate Sine data
    GenerateSineValue(gHist, 1.0, 1.0, 1.0, 3.0);

    // Initial parameter guesses (Amplitude, Frequency, Phase, Offset)
    double params[4] = {1.0, 1.0, 0.0, 2.0}; 
    double errors[4];

    // Create a ROOT file to store results
    TFile *ff = new TFile("outplot/SineFits.root", "RECREATE");
    if (!ff || ff->IsZombie()) {
        std::cerr << "Error: Cannot create ROOT file." << std::endl;
        return 1;
    }

    // Create a list to store plots
    plotList = new TList();
    plotList->Add(gHist);

    // Canvas Setup
    gStyle->SetOptStat(0);
    style::SetGlobalStyle();
    style::fgkYTitleOffset = 1.6;
    style::IniColorCB();

    const Double_t currentLeft = 0.17; 
    const Double_t currentTop = 0.06; 
    const Double_t currentRight = 0.035;
    const Double_t currentBottom = 0.14;

    // Ensure the global canvas is created with appropriate margins
    gCanvas = new TCanvas("gCanvas", "Sine Fit", 800, 600);
    gPad->SetLeftMargin(currentLeft);
    gPad->SetRightMargin(currentRight);
    gPad->SetTopMargin(currentTop);
    gPad->SetBottomMargin(currentBottom);

    // Title
    gStyle->SetOptTitle(1);
    gStyle->SetTitleX(0.6);
    gStyle->SetTitleW(0.8);

    // Perform the fit
    TF1 *fitFunc = new TF1("fitFunc", "[0] * sin([1] * x + [2]) + [3]", 0, 40);
    
    IterativeFit(fitFunc, params, errors, 10, 1e-6);

    // Perform fit with ROOT's internal function and print results
    FitWithROOT(fitFunc, params, errors);

    // Extract parameters from the histogram
    Extraction(gHist);
    
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
