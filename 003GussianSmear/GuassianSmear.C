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
#include <TVirtualFFT.h>

#include "style.h"

<<<<<<< HEAD
// Global objects
TRandom3 gRand(123); 
TH1D *gHist = nullptr;  
TCanvas *gCanvas = nullptr;
TList *plotList = nullptr; 
=======
// Global pointers to the random number generator and histogram
TRandom3 gRand(123);
TH1D *gHist = nullptr;
TCanvas *gCanvas = nullptr;
TList *plotList = nullptr;
TMinuit *gMinuit;
>>>>>>> quarry

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
	    // Protection: check if yVal is negative
            if (yVal < 0) {
                std::cerr << "Error: Negative yVal encountered (yVal = " << yVal 
                          << ") at xVal = " << xVal << ". Check your parameters." << std::endl;
                return;  // Exit function on error
            }
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
  const double p0 = TMath::Abs(params[0]);
  const double p1 = TMath::Abs(params[1]); // Amplitude
  const double p2 = TMath::Abs(params[2]); // frequency
  const double p3 = TMath::Abs(params[3]); // phase
  const double uu = 1 / (gkxMax - gkxMin) *
    (1 + p1 / p2 * ((TMath::Cos(p2 * gkxMax + p3) - 
                     TMath::Cos(p2 * gkxMin + p3))));
  // const double uu = -(p2 - p1 * TMath::Cos(gkxMin * p2 + p3) + p1 * TMath::Cos(gkxMax * p2 + p3))/((gkxMin - gkxMax) * p2);

  //substitute p4 for params[3] in the sine function expression
  return p0 * (p1 * TMath::Sin(p2 * xx + p3) + uu);
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

<<<<<<< HEAD
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
=======

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

    double result = (gIntegrator->Integral(x_min, x_max));
  
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
    
>>>>>>> quarry

    ff = chi2;
}

// Perform the fit using Minuit
<<<<<<< HEAD
int SingleFit(TF1* const fitFunc, double* const params, double* const errors, double& chi2) {
    TMinuit Minuit(4);
    Minuit.SetFCN(fcn);

    Minuit.DefineParameter(0, "Amplitude", params[0], 0.1, 0.1, 1e4);
    Minuit.DefineParameter(1, "Frequency", params[1], 0.01, 0.5, 5.0);
    Minuit.DefineParameter(2, "Phase", params[2], 0.1, -5, 5);
    Minuit.DefineParameter(3, "Offset", params[3], 0.1, -10, 1e4);
=======
int SingleFit(TF1* fitFunc, double* const params, double* const errors, double& chi2) {
    TMinuit Minuit(5);
    Minuit.SetFCN(fcn);

    Minuit.DefineParameter(0, "Normalisation", params[0], 0.1, 0, 0);
    Minuit.DefineParameter(1, "Amplitude", params[1], 0.1, 0, 0);  
    Minuit.DefineParameter(2, "Frequency", params[2], 0.1, 0, 0); 
    Minuit.DefineParameter(3, "Phase", params[3], 0.5, 0, 0); 
    Minuit.DefineParameter(4, "Sigma", params[4], 1, 0, 0); 

    // Minuit.DefineParameter(0, "Normalisation", params[0], 0.1 , 0, 0);
    // Minuit.DefineParameter(1, "Amplitude", params[1], 0.1 , 0, 0);  //Amplitude
    // Minuit.DefineParameter(2, "Frequency", params[2], 0.1 , 0, 0); //Frequency
    // Minuit.DefineParameter(3, "Phase", params[3], 0.5, 0, 0); //Phase
    // Minuit.DefineParameter(4, "Sigma", params[4], 0.1, 0, 0); //Sigma
>>>>>>> quarry

    Minuit.Migrad();  // Perform the minimization

    // Retrieve fitted parameters
<<<<<<< HEAD
    for (int ii = 0; ii < 4; ii++) {
=======
    for (int ii = 0; ii < 5; ii++) {
>>>>>>> quarry
        Minuit.GetParameter(ii, params[ii], errors[ii]);
    }

    // Update the fit function parameters
<<<<<<< HEAD
    fitFunc->SetParameters(params[0], params[1], params[2], params[3]);
=======
    fitFunc->SetParameters(params[0], params[1], params[2], params[3], params[4]);
>>>>>>> quarry

    // Retrieve fitting status
    double amin, edm, errdef;
    int nvpar, nparx, icstat;
    Minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    chi2 = amin;

    return Minuit.GetStatus();  // Return fit status
}

<<<<<<< HEAD
// Plot the fit at each iteration and store it in the TList
void PlotFit(const int iteration, TF1 * const fitFunc, const double *params, const double *errors) {
    if (!gCanvas) {
        gCanvas = new TCanvas(Form("c1_iter_%d", iteration), "Smear Fit", 800, 600);
=======

// Plot the fit at each iteration and store it in the TList
void PlotSineFit(const int iteration, TF1* fitFunc, const double* params, const double* errors) {
    if (!gCanvas) {
        gCanvas = new TCanvas(Form("c1_iter_%d", iteration), "Sine Fit", 800, 600);
>>>>>>> quarry
    } else {
        gCanvas->cd();  // Make sure to use the existing canvas
    }

    gCanvas->Clear();  // Clear the canvas for new plot

<<<<<<< HEAD
    gHist->SetTitle("SmearSine Fitting");
=======
    // Set title and axis labels
    gHist->SetTitle("Sine Fitting");
    gStyle->SetOptTitle(1);
    gHist->GetXaxis()->SetTitle("x-axis");
    gHist->GetYaxis()->SetTitle("y-axis");
    gHist->GetYaxis()->SetRangeUser(0, 5e3);
>>>>>>> quarry

    // Draw histogram and fit function
    style::ResetStyle(gHist);
    gHist->SetLineColor(kBlue);
    gHist->SetLineWidth(3);
    gHist->SetFillStyle(0);
<<<<<<< HEAD
    gHist->Draw();
    
    fitFunc->SetLineColor(kRed);
=======
    gHist->Draw();  // Draw after setting the title and axis label
    
    fitFunc->SetLineColor(kOrange);
>>>>>>> quarry
    fitFunc->SetLineWidth(3);
    fitFunc->SetNpx(1e5);
    fitFunc->Draw("same");

    /////////////////////////////////////////
    // double params1[4] = {5, 0.01, 2, 0};
    // auto sineFunction = [&](double* x, double* par) {
    //   double sine = (2e3 / 0.02) * (Sine(x[0], par) + 1/50);
    //   return sine;
    // };  
        
    // Legend
    const double x0 = 0.2;
    const double y0 = 0.6;
    const double x1 = x0 + 0.34;
    TLegend *lg = new TLegend(x0, y0, x1, 0.92);
    style::ResetStyle(lg, 0.18, 0.68);
    lg->SetTextAlign(12);

<<<<<<< HEAD
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
=======
    lg->AddEntry(gHist, "Generated Sine", "l");
    lg->AddEntry(fitFunc, "Fitted Sine", "l");

    lg->Draw();

    // Print out the plots
    gCanvas->Print(Form("outplot/SineFit_%d.png", iteration));

    delete lg;
}


// Perform iterative fitting with convergence check
void IterativeFit(double* const params, double* const errors, const int maxIterations, const double convergenceThreshold) {
>>>>>>> quarry
    double previousChi2 = 1e10;
    double currentChi2 = 0.0;

    for (int ii = 0; ii < maxIterations; ii++) {
<<<<<<< HEAD
        // Perform a single fit
        const int status = SingleFit(fitFunc, params, errors, currentChi2);

	// Output the chi-squared value for this iteration
        std::cout << "////////////////////////////Iteration " << ii + 1 << ": Chi2 = " << currentChi2 << "////////////////////////////" << std::endl;

        // Plot and store the results for this iteration
        PlotFit(ii + 1, fitFunc, params, errors);
=======

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
>>>>>>> quarry

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
<<<<<<< HEAD
    }
}

int main() {

     // Redirect ROOT output to the log file
=======

        delete fitFunc; // Clean up the TF1 object
    }
}


//Extract data from the histogram directly with FFT
double FFT(TH1D *const hist) {
    if (!hist) {
        std::cerr << "Error: Histogram not properly initialized!" << std::endl;
        return 1;
    }

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
        return 1;
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
    std::vector<double> magnitudes(n_size / 2);
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
    const double freqStep = hist->GetXaxis()->GetBinWidth(1);  // Bin width in x-axis units
    const double dominantFreq = (maxBin * 1.0) / (n_size * freqStep); // Make sure division is correct

    // Output the results
    std::cout << "maxBin: " << maxBin << std::endl;
    // std::cout << "Offset: " << offset << " ± " << offsetErr << std::endl;
    std::cout << "Dominant Frequency: " << dominantFreq << std::endl;

    // Clean up
    delete[] in;
    delete[] real;
    delete[] imaginary;
    delete fft;

    return dominantFreq;
}

double Norsine(double *xx, double *par)
{
  const double p0 = TMath::Abs(par[0]);
  const double p1 = TMath::Abs(par[1]);
  const double p2 = TMath::Abs(par[2]);
  const double p3 = TMath::Abs(par[3]);

  // const  double uu = -1 / (gkxMax - gkxMin) *
  //   (1 + p1 / p2 * ((TMath::Cos(p2 * (gkxMax + p3)) - 
  //                    TMath::Cos(p2 * (gkxMin + p3)))));

  const double uu = -(p2 - p1 * TMath::Cos(gkxMin * p2 + p3) + p1 * TMath::Cos(gkxMax * p2 + p3))/((gkxMin - gkxMax) * p2);
  
  return p0*(p1*TMath::Sin(p2*xx[0]+p3) + uu);

}


// Extract sine initial guesses  using ROOT's internal fitting
void InitialGuess(double* params) { 
    //  auto sineFunc = [&](double* x, double* par) {
    //    return TMath::Abs(par[0]) * Sine(x, par);  
    // };

    // TF1 *sineFit = new TF1("sineFit", sineFunc, 0, 50, 3);

    // Define a TF1 object with the same sine formula
    TF1* sineFit = new TF1( "norsine", Norsine, gkxMin, gkxMax, 4);
    
    // Initial guesses
    const double norm = gHist->Integral(0,10000, "width");
    const double maxVal = gHist->GetMaximum() / norm;
    const double freq = TMath::TwoPi() * FFT(gHist);
    const double fitParams[] = {norm, maxVal, freq, 0, 1};
    sineFit->SetParameters(fitParams);

    // // Adjust and relax parameter limits
    // sineFit->SetParLimits(0, 0, 1e9);  // Wider amplitude range
    // sineFit->SetParLimits(1, 0, 1e3);  // Frequency
    // sineFit->SetParLimits(3, 0, 2 * TMath::Pi()); // Phase

    // sineFit->SetParLimits(0, 0, 0);  // Wider amplitude range
    // sineFit->SetParLimits(1, 0, 0);  // Frequency
    // sineFit->SetParLimits(3, 0, 0); // Phase

    gHist->Fit(sineFit, "V");  // Verbose output

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
    std::cout << "Normalisation: " << params[0] << " ± " << sineFit->GetParError(0) << std::endl;
    std::cout << "Amplitude: " << params[1] << " ± " << sineFit->GetParError(1) << std::endl;
    std::cout << "Frequency: " << params[2] << " ± " << sineFit->GetParError(2) << std::endl;
    std::cout << "Phase: " << params[3] << " ± " << sineFit->GetParError(3) << std::endl;

    delete sineFit;
}

int main() {
    // Redirect ROOT output to the log file
>>>>>>> quarry
    const Int_t status = gSystem->RedirectOutput("outplot/see.log", "w");  // "w" to overwrite the file
    if (status != 0) {
        std::cerr << "Error: Unable to redirect output to file." << std::endl;
        return 1;
    }
<<<<<<< HEAD

    // Canvas Setup
=======
  
>>>>>>> quarry
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
<<<<<<< HEAD
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
=======
    const double params[5] = {0.5, 2, 2, 1, 1}; // amp  freq  phase offset sigma 
    
    // const double params[5] = {0.5, 2, 5 * TMath::TwoPi() / 50, 1, 1}; // One period


    // Perform smearing of the sine signal with Gaussian noise
    SmearSignal(gHist, params);

        /////////////////////////////////////////
    // double params1[4] = {5, 0.01, 2, 0};
    // auto sineFunction = [&](double* x, double* par) {
    //   double sine = (2e3 / 0.02) * (Sine(x[0], par) + 1/50);
    //   return sine;
    // };

    double guess[5];
    double errors[5];

    InitialGuess(guess);

    IterativeFit(guess, errors, 20, 1e-9); // Iteration number, Convergence threshold

    // // Extract parameters from the histogram
    // FFT(gHist);

    // Restore output to the default streams
    gSystem->RedirectOutput(0);

>>>>>>> quarry
    delete gHist;
    delete gCanvas;

    return 0;
}
