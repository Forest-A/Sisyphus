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
#include <TStopwatch.h>

#include "style.h"

// Global pointers to the random number generator and histogram
TRandom3 gRand(123);
TH1D *gHist = nullptr;
TCanvas *gCanvas = nullptr;
TList *plotList = nullptr;  

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
    
    for (int ii = 0; ii < 1e4; ii++) {
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
    return params[0] * sin(params[1] * xx + params[2]) + params[3];
}

// Gaussian function
double Gaussian(double xx, double sigma) {
    return exp(-0.5 * xx * xx / (sigma * sigma));
}

// Continuous convolution of sine and Gaussian
double convolvedFunction(const double xx, const double *params) {
    ROOT::Math::IntegratorOneDim integrator;
    // ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Guass");

    auto integrand = [&](double x_prime) -> double {
        double sine_val = Sine(x_prime, params);          
        double gauss_val = Gaussian(xx - x_prime, params[4]);  
        return sine_val * gauss_val;
    };

    integrator.SetRelTolerance(1e-4);
    
    // Set lambda function to the integrator
    integrator.SetFunction(integrand);

    // Adjust integration limits based on sigma
    const double x_min = xx - 15 * params[4];  
    const double x_max = xx + 15 * params[4];  
    
    // Perform continuous integration
    double result = integrator.Integral(x_min, x_max);

    return result;
}

// Define a ROOT-compatible function using TF1 for the convolution
TF1* CreateConvolution(const double* params) {
    auto convolved = [&](double* x, double* par) {
        return convolvedFunction(x[0], par);
    };

    // Define convolution function in ROOT's TF1 class
    TF1 *fitFunc = new TF1("fitFunc", convolved, 0, 50, 5); // xmin = 0, xmax = 100
    fitFunc->SetParameters(params);
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

    Minuit.DefineParameter(0, "Amplitude", params[0], 1, 0.1, 1e7);  // Amplitude
    Minuit.DefineParameter(1, "Frequency", params[1], 0.01, 0.0, 1e4); // Frequency
    Minuit.DefineParameter(2, "Phase", params[2], 0.1, 0.0, 1e3);      // Phase
    Minuit.DefineParameter(3, "Offset", params[3], 0.1, 0, 1e3);     // Offset
    Minuit.DefineParameter(4, "Sigma", params[4], 0.1, 0.0, 1e3);       // Sigma

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
    gHist->GetYaxis()->SetRangeUser(0, 60);

    // Draw histogram and fit function
    style::ResetStyle(gHist);
    gHist->SetLineColor(kBlue);
    gHist->SetLineWidth(1);
    gHist->SetFillStyle(0);
    gHist->Draw();  // Draw after setting the title and axis labels

    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(2);
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
        std::cout << "////////////////////////////Iteration " << ii + 1 << ": Chi2 = " << currentChi2 << "////////////////////////////" << std::endl;

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
    const double params[5] = {5, 1, 0, 10, 1}; // Amplitude, Frequency, Phase, Offset, Sigma

    // Perform smearing of the sine signal with Gaussian noise
    SmearSignal(gHist, params);

    double guess[5] = {5.0, 1.0, 0.0, 10.0, 1.0}; 
    double errors[5];

    IterativeFit(guess, errors, 10, 1e-8); // Iteration number, Convergence threshold

    // Restore output to the default streams
    gSystem->RedirectOutput(0);
    
    delete gHist;
    delete gCanvas;

    return 0;
}


// #include <TFile.h>
// #include <TH1D.h>
// #include <TF1.h>
// #include <TCanvas.h>
// #include <TMath.h>
// #include <TRandom3.h>
// #include <TLegend.h>
// #include <Math/Integrator.h>
// #include <iostream>
// #include <cmath>

// // Sine function
// double Sine(double xx, const double *params) {
//     return params[0] * sin(params[1] * xx + params[2]) + params[3];
// }

// // Normalized Gaussian function
// double Gaussian(double xx, double sigma) {
//     return (1.0 / (sigma * sqrt(2 * TMath::Pi()))) * exp(-0.5 * xx * xx / (sigma * sigma));
// }

// // Continuous convolution of sine and Gaussian
// double convolvedFunction(const double xx, const double *params) {
//       ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");
//      ROOT::Math::IntegratorOneDim integrator;
//     integrator.SetRelTolerance(1e-6);  // Increase accuracy
    
    
//     // Define the integrand for convolution
//     auto integrand = [&](double x_prime) -> double {
//         double sine_val = Sine(xx - x_prime, params);          
//         double gauss_val = Gaussian( x_prime, params[4]); 
//         return sine_val * gauss_val;
//     };

//     // Set the lambda function to the integrator
//     integrator.SetFunction(integrand);

//     // Adjust integration limits based on sigma
//     const double x_min = xx + 10 * params[4];
//     const double x_max = xx - 10 * params[4];
    
//     // Perform continuous integration
//     double result = integrator.Integral(x_min, x_max);
//     return result;
// }

// // Define a ROOT-compatible function using TF1 for the convolution
// TF1* CreateConvolutionFunction(const double* params) {
//     // Lambda function to wrap the convolution function for TF1
//     auto convolved = [&](double* x, double* par) {
//         return convolvedFunction(x[0], par);
//     };

//     // Define convolution function in ROOT's TF1 class
//     TF1 *fitFunc = new TF1("fitFunc", convolved, 0, 20, 5); // xmin = 0, xmax = 20
//     fitFunc->SetParameters(params);
//     return fitFunc;
// }

// int main() {
//     // Parameters for the sine function and Gaussian noise
//   double params[5] = {5, 1, 0, 0, 2};  // Amplitude, Frequency, Phase, Offset, Sigma

//     // Create canvas for plotting
//     TCanvas *canvas = new TCanvas("canvas", "Convolved Sine and Gaussian", 800, 600);
    
//     // Create the convolution function and plot it
//     TF1* convolutionFunc = CreateConvolutionFunction(params);
    
//     convolutionFunc->SetLineColor(kRed);
//     convolutionFunc->SetTitle("Convolved Sine and Gaussian; x-axis; y-axis");
//     convolutionFunc->SetNpx(1e3);  // Increase number of points for smooth plot
//     convolutionFunc->Draw();

//     // Save the plot as a PNG file
//     canvas->SaveAs("outplot/convolution_plot.png");

//     // Display canvas
//     canvas->Update();
    
//     // Clean up
//     delete canvas;
//     delete convolutionFunc;

//     return 0;
// }


// #include <TFile.h>
// #include <TH1D.h>
// #include <TF1.h>
// #include <TCanvas.h>
// #include <TMath.h>
// #include <TRandom3.h>
// #include <TLegend.h>
// #include <TStyle.h>
// #include <TF1Convolution.h>

// int main() {
//     // Style and canvas setup
//     gStyle->SetOptStat(0);
//     TCanvas *canvas = new TCanvas("canvas", "Sine + Gaussian Convolution", 800, 600);

//     // Define the sine function
//     TF1 *sine = new TF1("sine", "[0] * TMath::Sin([1] * x + [2]) + [3]", 0, 100);
//     sine->SetParameters(5, 1, 0, 20);  // Amplitude, Frequency, Phase, Offset
//     sine->SetLineColor(kBlue);
//     sine->SetTitle("Sine and Gaussian Convolution; x-axis; y-axis");

//     // Define the Gaussian function (centered at 0)
//     TF1 *gaussian = new TF1("gaussian", "TMath::Gaus(x, 0, [0])", -10, 10);
//     gaussian->SetParameter(0, 4);  // Mean, Sigma
//     gaussian->SetLineColor(kGreen);

//     // Create a TF1Convolution object that convolves sine and Gaussian
//     TF1Convolution *conv = new TF1Convolution(sine, gaussian, 0, 100, true);
//     conv->SetNofPointsFFT(1000); // Set the number of points for FFT convolution

//     // Create a new TF1 to represent the convolution
//     TF1 *convolution = new TF1("convolution", *conv, 0, 100, conv->GetNpar());
//     convolution->SetNpx(10000);  // Increase the number of points for better resolution
//     convolution->SetLineColor(kRed);

//     // Set the same parameters for the convolution as in the original functions
//     convolution->SetParameters(sine->GetParameters());
//     convolution->SetParameter(4, gaussian->GetParameter(0));  // Sigma

//     // Draw the original sine, Gaussian, and convolution
//     // sine->Draw();
//     // gaussian->Draw("same");
//     convolution->Draw();

//     // Create a legend
//     TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
//     legend->AddEntry(sine, "Sine", "l");
//     legend->AddEntry(gaussian, "Gaussian", "l");
//     legend->AddEntry(convolution, "Convolution", "l");
//     legend->Draw();

//     // Save the canvas
//     canvas->SaveAs("sine_gaussian_convolution.png");

//     // Clean up
//     delete sine;
//     delete gaussian;
//     delete convolution;
//     delete conv;
//     delete canvas;

//     return 0;
// }


