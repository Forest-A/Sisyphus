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
#include <Math/Integrator.h>
#include <Math/IntegratorOptions.h>
#include <TStopwatch.h>

#include "style.h"

// Global pointers to the random number generator and histogram
TRandom3 gRand(123);
TH1D *gHist = nullptr;
TCanvas *gCanvas = nullptr;
TList *plotList = nullptr;
TMinuit *gMinuit;

// Global constants for the integrator and axis limits
const double gkxMax = 50.0;
const double gkxMin = 0.0;
ROOT::Math::IntegratorOneDim *gIntegrator = nullptr;

// Global variable for params and xx
const double* gParams = nullptr;
double gXX = 0;

// Generate sine values
void SmearSignal(TH1D *const hist, const double *params) {
    if (!hist) {
        std::cerr << "Error: Histogram not properly initialized!" << std::endl;
        return;
    }

    const double amp = params[0];
    const double freq = params[1];
    const double phase = params[2];
    const double offset = params[3];
    const double sigma = params[4];
    
    const double xMin = hist->GetXaxis()->GetXmin();
    const double xMax = hist->GetXaxis()->GetXmax(); 
    const double maxY = amp + offset;
    
    for (int ii = 0; ii < 1e6; ii++) {
        double EE = 0;
        // Acceptance-Rejection sampling
        while (true) {
            const double xVal = gRand.Uniform(xMin, xMax);
            const double yVal = amp * TMath::Sin(freq * xVal + phase) + offset;
            const double rr = gRand.Uniform(0, maxY); 

            // Acceptance criterion: accept the value if rr < yVal
            if (rr < yVal) {
                EE = xVal;
                break;
            }
        }

        double energy = gRand.Gaus(EE, sigma);
        gHist->Fill(energy);
        // gHist->SetBinError(ii, sigma);
    }
}

// Sine function
double Sine(double xx, const double *params) {
  double p1 = params[1]; // Amplitude
  double p2 = params[2]; // frequency
  double p3 = params[3]; // phase
  double uu = 1 / (gkxMax - gkxMin) *
    (1 + p1 / p2 * ((TMath::Cos(p2 * gkxMax + p3) - 
                     TMath::Cos(p2 * gkxMin + p3)))); 

  //substitute p4 for params[3] in the sine function expression
  return p1 * TMath::Sin(p2 * xx + p3) + uu;
}

// Gaussian function
double Gaussian(double xx, double sigma) {
    return TMath::Gaus(xx, 0.0, sigma, true);
}

// Integrand function 
double Integrand(double x_prime) {
    double sine_val = Sine (x_prime, gParams);
    double gauss_val = Gaussian(gXX - x_prime, gParams[4]);

    return sine_val * gauss_val;
}

// Initialize the integrator and set the function once
void InitializeIntegrator() {
    if (!gIntegrator) {
        gIntegrator = new ROOT::Math::IntegratorOneDim(ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
        if (gIntegrator) {
            gIntegrator->SetRelTolerance(1e-4);
            gIntegrator->SetFunction(Integrand);  // Set the function only once
        }
        else {
            std::cerr << "Error: Failed to initialize the integrator!" << std::endl;
        }
    }
}


// Continuous convolution of sine and Gaussian
double convolvedFunction(const double xx, const double *params) {
    if (!gIntegrator) {
        InitializeIntegrator();
    }

    // Set global variables for the integrand
    gParams = params;
    gXX = xx;

    const double x_min = xx - 5 * TMath::Abs(params[4]);
    const double x_max = xx + 5 * TMath::Abs(params[4]);

    double result = params[0] * gIntegrator->Integral(x_min, x_max);
  
    if (gIntegrator->Status() != 0) {
        std::cerr << "Warning: Integration failed for xx = " << xx << std::endl;
        return 0;
    }

    return result;
}


// Define a ROOT-compatible function using TF1 for the convolution
TF1* CreateConvolution(const double* params) {
    auto convolved = [&](double* x, double* par) {
        return convolvedFunction(x[0], par);
    };

    // Define convolution function in ROOT's TF1 class
    TF1 *fitFunc = new TF1("fitFunc", convolved, 0, 50, 5);
    fitFunc->SetParameters(params);
    fitFunc->SetNpx(1e5); 
    return fitFunc;
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
        const double fitValue = convolvedFunction(xx, par);

        double binError = gHist->GetBinError(ii);
        if (binError <= 0) {
            binError = 1.0;  // Assign a fallback value to avoid division by zero
        }
        chi2 += TMath::Power((fitValue - data), 2) / TMath::Power(binError, 2);
    }
    

    ff = chi2;
}

// Perform the fit using Minuit
int SingleFit(TF1* fitFunc, double* const params, double* const errors, double& chi2) {
    TMinuit Minuit(5);
    Minuit.SetFCN(fcn);

    Minuit.DefineParameter(0, "Normalisation", params[0], 0.1 * params[0], 0, 10 * params[0]);
    Minuit.DefineParameter(1, "Amplitude", params[1], 0.1 * params[1], 0, 10 * params[1]);  // Amplitude
    Minuit.DefineParameter(2, "Frequency", params[2], 0.1 * params[2], 0.1, 10 * params[2]); // Frequency
    Minuit.DefineParameter(3, "Phase", params[3], 0.1, -TMath::Pi(), TMath::Pi()); // Phase
    Minuit.DefineParameter(4, "Sigma", params[4], 0.1 * params[4], 0.1, 10 * params[4]); // Sigma

    Minuit.Migrad();  // Perform the minimization

    // Retrieve fitted parameters
    for (int ii = 0; ii < 5; ii++) {
        Minuit.GetParameter(ii, params[ii], errors[ii]);
    }

    // Update the fit function parameters
    fitFunc->SetParameters(params[0], params[1], params[2], params[3], params[4]);

    // Retrieve fitting status
    double amin, edm, errdef;
    int nvpar, nparx, icstat;
    Minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    chi2 = amin;

    return Minuit.GetStatus();  // Return fit status
}

// Plot the fit at each iteration and store it in the TList
void PlotSineFit(const int iteration, TF1* fitFunc, const double* params, const double* errors) {
    if (!gCanvas) {
        gCanvas = new TCanvas(Form("c1_iter_%d", iteration), "Sine Fit", 800, 600);
    } else {
        gCanvas->cd();  // Make sure to use the existing canvas
    }

    gCanvas->Clear();  // Clear the canvas for new plot

    // Set title and axis labels
    gHist->SetTitle("Sine Fitting");
    gStyle->SetOptTitle(1);
    gHist->GetXaxis()->SetTitle("x-axis");
    gHist->GetYaxis()->SetTitle("y-axis");
    gHist->GetYaxis()->SetRangeUser(0, 5e3);

    // Draw histogram and fit function
    style::ResetStyle(gHist);
    gHist->SetLineColor(kBlue);
    gHist->SetLineWidth(3);
    gHist->SetFillStyle(0);
    gHist->Draw();  // Draw after setting the title and axis label
    
    fitFunc->SetLineColor(kOrange);
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

    // Print out the plots
    gCanvas->Print(Form("outplot/SineFit_%d.png", iteration));

    delete lg;
}


// Perform iterative fitting with convergence check
void IterativeFit(double* const params, double* const errors, const int maxIterations, const double convergenceThreshold) {
    double previousChi2 = 1e10;
    double currentChi2 = 0.0;

    for (int ii = 0; ii < maxIterations; ii++) {

        TStopwatch timer;  // Start timing each iteration
	timer.Start();
      
        // Create the convolution function for each iteration
        TF1 *fitFunc = CreateConvolution(params);

        // Perform a single fit
        const int status = SingleFit(fitFunc, params, errors, currentChi2);

        // Output the chi-squared value for this iteration
        std::cout << "////////////////////////////Iteration " << ii + 1 << ": Chi2 = " << currentChi2 << "///////////////////////////" << std::endl;

        // Plot and store the results for this iteration
        PlotSineFit(ii + 1, fitFunc, params, errors);

	// Stop the timer and print the time consumed per loop
	timer.Stop();
	std::cout << "Iteration " << ii + 1 << " took: " << timer.RealTime() << "seconds (real time), " 
                  << timer.CpuTime() << " seconds (CPU time)." << std::endl;
	std::cout << "//////////////////////////////////////////////////////////////////////////////////" << std::endl;

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

        delete fitFunc; // Clean up the TF1 object
    }
}


// Extract sine initial guesses  using ROOT's internal fitting
void InitialGuess(double* params) { 
     auto sineFunc = [&](double* x, double* par) {
        return par[0] * Sine(x[0], par);  // Assuming Sine is already defined elsewhere
    };
    
    TF1 *sineFit = new TF1("sineFit", sineFunc, 0, 50, 4);
    
    // Initial guesses based on histogram data
    double maxVal = gHist->GetMaximum();
    sineFit->SetParameters(maxVal, 2, 3, 1);  // Guesses for amplitude, frequency, phase

    // // Adjust and relax parameter limits
    // sineFit->SetParLimits(0, 0.5 * maxVal, 2 * maxVal);  // Wider amplitude range
    // sineFit->SetParLimits(1, 0, 100);  // Frequency
    // sineFit->SetParLimits(2, -10, 10); // Phase

    gHist->Fit(sineFit, "M");  // Minuit for better precision

    // // Extract fitted parameters and their errors
    // for (int ii = 0; ii < 3; ii++) {
    //     params[ii+1] = sineFit->GetParameter(ii);
    // }
    params[0] = sineFit->GetParameter(0);
    params[1] = sineFit->GetParameter(1);
    params[2] = sineFit->GetParameter(2);
    params[3] = sineFit->GetParameter(3);
    params[4] = 1e-6;  // Initial guess for sigma

    // Log the fit results and errors
    std::cout << "Initial Fit Results:" << std::endl;
    std::cout << "Normalisation: " << params[0] << " ± " << sineFit->GetParError(3) << std::endl;
    std::cout << "Frequency: " << params[2] << " ± " << sineFit->GetParError(1) << std::endl;
    std::cout << "Phase: " << params[3] << " ± " << sineFit->GetParError(2) << std::endl;

    delete sineFit;
}

int main() {
    // Redirect ROOT output to the log file
    const Int_t status = gSystem->RedirectOutput("outplot/see.log", "w");  // "w" to overwrite the file
    if (status != 0) {
        std::cerr << "Error: Unable to redirect output to file." << std::endl;
        return 1;
    }
  
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

    gHist = new TH1D("gHist", "Generated Noisy Sine Signal", 500, 0, 50);

    // Parameters for sine function and Gaussian noise
    const double params[5] = {5, 2, 1, 1, 1}; // Normalisation, Amplitude, Frequency, Phase, Sigma

    // Perform smearing of the sine signal with Gaussian noise
    SmearSignal(gHist, params);

    double guess[5];
    double errors[5];

    InitialGuess(guess);

    IterativeFit(guess, errors, 20, 1e-6); // Iteration number, Convergence threshold

    // Restore output to the default streams
    gSystem->RedirectOutput(0);

    delete gHist;
    delete gCanvas;

    return 0;
}
