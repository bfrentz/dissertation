/*******************************************
 * The analysis program from the 
 * 14N(p,g)15O implanted target tests
 * Usage with root:
 *
 *       Prep:
 *       Combine all histograms from separate files into one root file
 *       Adjust arrays and methods in this file for the correct number of histograms
 *       Test out fitting by running the methods individually to find your bounds for summing
 *
 *       Workflow:
 *       Change energy of interest in main() method
 *       (in root) .L aggregateAnalysis.cpp
 *       main()  
 *
 *       Loops over all histograms in file and fits for the specific peak to which you're interested, 
 *           first finding the centroid from the TSpectrum class. Then it outputs a figure of the fit,
 *           as well as the statistics about net area, centroid, and centroid error to a csv file.
 *
 * Created: 9 August 2019
 * Updated: 9 October 2019
 *
 * Bryce Frentz (bfrentz@nd.edu)
 ********************************************/


#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdio>
#include "TCanvas.h"
#include "TRandom2.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCut.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TKey.h"

using namespace std;


// Class for storing all fit results into one parameter
struct fitResults {
	double cts;
    double centroid;
	double unc;
    double unc2;
};

// Defines the piecewise line
// Fits background around peak by taking only user defined background
//     regions and excluding other areas between for BG fitting
// The surrounding area is fit with a cubic term (to account for
//     the compton background changing the height around the peak)
//     and the rejection under the peak is done with fline
Bool_t reject;
Double_t xExcludeMin;
Double_t xExcludeMax;
Double_t fLine(double *x, double *par)
{
    if (reject && x[0] > xExcludeMin && x[0] < xExcludeMax) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0];
}

Double_t fPoly(double *x, double *par)
{
    if (reject && x[0] > xExcludeMin && x[0] < xExcludeMax) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}


// Uses TSpectrum to find peaks with relative heights as specified and outputs peak number and centroid
double PeakFinder(int runNum, TH1D* h0, int peak)
{


    //Construct TSpectrum with max number of peaks 50
    TSpectrum *spec = new TSpectrum(250);
    
    //Find the peaks, nSigma is the approximate sigma of the peaks, threshold discards any peak with amplitude less than thresh*highest_peak
    Double_t nSigma = 1;
    Double_t thresh = 2e-3;
    Int_t nFound = spec->Search(h0, nSigma, "", thresh);
    
    //Print out the number of found peaks
    cout << "Number of Peaks Found: " << nFound << endl;
    
    
    //Get the X and Y positions of the peaks and sort them into arrays
    Double_t *xpos = spec->GetPositionX();
    Double_t *ypos = spec->GetPositionY();
    Double_t xswap = 0;
    Double_t yswap = 0;
    Int_t count = 0;
    
    while (count < nFound) {
        for (Int_t i = count; i < nFound; i++) {
            if (xpos[i] < xpos[count]) {
                xswap = xpos[count];
                xpos[count] = xpos[i];
                xpos[i] = xswap;
                
                yswap = ypos[count];
                ypos[count] = ypos[i];
                ypos[i] = yswap;
            }
        }
        count++;
    }
    
    //Print out the location of the different peaks
    for (Int_t i=1; i <= nFound; i++) {
        printf("%f\n", xpos[i-1]);
        if (xpos[i-1] > (peak-20) and xpos[i-1] < (peak+40)){
            return xpos[i-1];
        }
    }
    return 0;
}


//linear+gaussian+skew+step guess
// Gaussian on top of radware background
// Quick single fit to show how close parameters are
fitResults fitSkewGauss(TH1D* h1, double cent, double width, int xMin, int xMax)
{
    // Define functions and set initial parameters
    double pi = 3.14159265359;
    //TH1F* h1 = (TH1F*)inFile->Get(hist);

    // Create the Canvas
    //TCanvas *cFit = new TCanvas("cFit", "Skew Gaussian Fit", 700, 500);

    // Function definitions
    TF1* fLin = new TF1("fLin","[0]+[1]*x+[2]*TMath::Erfc((x-[3])/(sqrt(2.0)*[4]))*[5]/100", xMin, xMax);
    TF1* fLinGauss = new TF1("fLinGauss","[0]+[1]*x+[2]*((1-[5]/100)*exp(-((x-[3])/(sqrt(2.0)*[4]))**2)+TMath::Erfc((x-[3])/(sqrt(2.0)*[4])+[4]/(sqrt(2.0)*[6]))*[5]/100*exp((x-[3])/[6])+TMath::Erfc((x-[3])/(sqrt(2.0)*[4]))*[7]/100)");
    fLinGauss->SetParameters(0,0,0,cent,width,1,1,1);
    fLinGauss->FixParameter(3,cent);
    fLinGauss->FixParameter(4,width);
    fLinGauss->FixParameter(5,1);
    fLinGauss->FixParameter(6,1);
    fLinGauss->FixParameter(7,1);
    h1->Fit("fLinGauss","BLSQ","",xMin,xMax);   
    fLinGauss->SetParLimits(5,0,100.0);  //R
    fLinGauss->SetParLimits(6,0,1000.0); //beta
    fLinGauss->SetParLimits(7,0,100.0);  //Step
    fLinGauss->SetParLimits(3,cent-10,cent+10);
    fLinGauss->SetParLimits(4,0,width+10);


    // Iterative fit 10 times
    // Repeats the fit to reduce error each time and find best fit
    for(int i=1; i<=10; i++)
    {
        h1->Fit("fLinGauss","BLQM","",xMin,xMax);
    }
    h1->Fit("fLinGauss","BLSQ","",xMin,xMax);
    h1->GetXaxis()->SetRangeUser(xMin-10,xMax+10);

    // Store final fit as pointer for later reference
    TFitResultPtr r = h1->Fit("fLinGauss","LSQ","",xMin,xMax);

    // Set background fucntion parameters from the full fit function
    fLin->SetParameter(0, fLinGauss->GetParameter(0));
    fLin->SetParameter(1, fLinGauss->GetParameter(1));
    fLin->SetParameter(2, fLinGauss->GetParameter(2));
    fLin->SetParameter(3, fLinGauss->GetParameter(3));
    fLin->SetParameter(4, fLinGauss->GetParameter(4));
    fLin->SetParameter(5, fLinGauss->GetParameter(7));

    // Set the background function parameter errors from the full fit function errors
    fLin->SetParError(0, fLinGauss->GetParError(0));
    fLin->SetParError(1, fLinGauss->GetParError(1));
    fLin->SetParError(2, fLinGauss->GetParError(2));
    fLin->SetParError(3, fLinGauss->GetParError(3));
    fLin->SetParError(4, fLinGauss->GetParError(4));
    fLin->SetParError(5, fLinGauss->GetParError(7));
    Double_t c0_err = fLinGauss->GetParError(0);
    Double_t c1_err = fLinGauss->GetParError(1);


    // Bin contents analysis
    // Setting the centroid and sigma from the fit
    Double_t centroid = fLinGauss->GetParameter(3);
    Double_t centroidError = fLinGauss->GetParError(3);
    Double_t sigma = fLinGauss->GetParameter(4);
    Double_t sigmaError = fLinGauss->GetParError(4);
    Double_t sigTau = sqrt(TMath::Power(sigma,2) + TMath::Power(fLinGauss->GetParameter(6),2));

    // Bounds for integrating background: (centroid +/- 4sigma)
    Int_t xMinBin = h1->FindBin(centroid - (4*sigTau));
    Int_t xMaxBin = h1->FindBin(centroid + (4*sigma))+1;
    Double_t xMinBoundary = h1->GetBinLowEdge(xMinBin);
    Double_t xMaxBoundary = h1->GetBinLowEdge(xMaxBin);

    // Sum bins within boundary
    Double_t binWidth = h1->GetXaxis()->GetBinWidth(xMinBin);
    Double_t sumArea = 0;
    Double_t sumError = 0;
    Double_t sumCenter = 0;
    Double_t sumVariance = 0;
    Double_t sumSigma = 0;
    Double_t sumFWHM = 0;
    Double_t totalCentroidError = 0;
    Int_t binCounter = xMinBin;
    do {
            sumArea += h1->GetBinContent(binCounter);
            sumCenter += h1->GetBinContent(binCounter)*h1->GetBinCenter(binCounter);
            binCounter++;
    }   while(binCounter < xMaxBin);

    // Using bin contents find the center of mass (centroid) - for comparison only
    sumCenter = sumCenter*binWidth;
    sumArea = sumArea*binWidth;
    sumError = sqrt(sumArea);
    sumCenter = sumCenter/sumArea;

    // Variance and sigma for bin contents
    binCounter = xMinBin;
    do{
            sumVariance += h1->GetBinContent(binCounter)*TMath::Power((h1->GetBinCenter(binCounter)-sumCenter),2);
            binCounter++;
    }   while(binCounter < xMaxBin);

    sumVariance = sumVariance*binWidth;
    sumVariance = sumVariance/sumArea;
    sumSigma = sqrt(sumVariance);
    sumFWHM = sumSigma*2.35482;


    // Integrating fit results
    Double_t bkgdIntegral = fLin->Integral(xMinBoundary, xMaxBoundary);
    Double_t bkgdIntegralError = sqrt(bkgdIntegral);
    //Double_t bkgdIntegralError = TMath::Sqrt( TMath::Power((xMaxBoundary - xMinBoundary)*c0_err,2.0) + TMath::Power(0.5*(xMaxBoundary - xMinBoundary)*(xMaxBoundary - xMinBoundary)*c1_err,2.0) );
    Double_t fnIntegral = fLinGauss->Integral(xMinBoundary, xMaxBoundary);
    Double_t fnIntegralError = fLinGauss->IntegralError( xMinBoundary, xMaxBoundary, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray() );
    Double_t netArea = fnIntegral - bkgdIntegral;
    Double_t netAreaError = TMath::Sqrt( TMath::Power(fnIntegralError,2.0) + TMath::Power(bkgdIntegralError,2.0) );


    // Sum-Background 
    //Double_t netArea = sumArea - bkgdIntegral;
    //Double_t netError = TMath::Sqrt( TMath::Power(sumError,2.0) + TMath::Power(bkgdIntegralError,2.0) );

    // Print results
    //r->GetCovarianceMatrix().Print();
    //printf("\n\nFit area = %10.6f+/-%10.6f\nBkgd area = %10.6f+/-%10.6f\nNet area (fit) = %10.6f+/-%10.6f\n",fnIntegral,fnIntegralError,bkgdIntegral,bkgdIntegralError);

    // Total centroid error:
    // Contribution from background and weighting of bin contents center of mass bs
    // unc = sigma(centroid) * diff(Area)/Area
    Double_t diffFactor = TMath::Sqrt(fnIntegral + bkgdIntegral);
    //cout << "DEBUG: diffFactor (sqrt) = " << diffFactor << endl;
    diffFactor = diffFactor / netArea;
    totalCentroidError = sumSigma * diffFactor;

    //cout << "DEBUG: diffFactor (division) = " << diffFactor << endl;
    cout << endl;

    totalCentroidError = TMath::Sqrt(totalCentroidError*totalCentroidError + centroidError*centroidError);

    // Print results
    cout << endl;
    cout << endl;
    cout << "+-------------------------------------------------+" << endl;
    //cout << "b (bkgd):                 " << fLin->GetParameter(0) << " +/- " << fLin->GetParError(0) << endl;
    //cout << "m (bkgd):                 " << fLin->GetParameter(1) << " +/- " << fLin->GetParError(1) << endl;
    //cout << "scale (bkgd):             " << fLin->GetParameter(2) << " +/- " << fLin->GetParError(2) << endl;
    //cout << "R    (fit):               " << fLinGauss->GetParameter(5) << " +/- " << fLinGauss->GetParError(5) << endl;
    cout << "Tau  (fit):               " << fLinGauss->GetParameter(6) << " +/- " << fLinGauss->GetParError(6) << endl;
    //cout << "Step (fit):               " << fLinGauss->GetParameter(7) << " +/- " << fLinGauss->GetParError(7) << endl;
    cout << "MinBoundary (bins):       " << xMinBoundary << /*"\t" << xMinBin << */ endl;
    cout << "MaxBoundary (bins):       " << xMaxBoundary << /*"\t" << xMaxBin << */  endl;
    cout << "+-------------------------------------------------+" << endl;
    cout << "Centroid  (fit):          " << centroid << " +/- " << centroidError << endl;
    cout << "Centroid (bins):          " << sumCenter << endl;
    cout << "+-------------------------------------------------+" << endl;
    cout << "Sigma  (fit):             " << sigma << " +/- " << sigmaError << endl;
    cout << "Sigma (bins):             " << sumSigma <<  endl;
    cout << "+-------------------------------------------------+" << endl;
    cout << "Tau  (fit):               " << fLinGauss->GetParameter(6) << " +/- " << fLinGauss->GetParError(6) << endl;
    //cout << "FWHM (bins):              " << sumFWHM << endl;
    cout << "+-------------------------------------------------+" << endl;
    cout << "Area (bkgd):              " << bkgdIntegral << " +/- " << bkgdIntegralError << endl;
    cout << "Area (bins):              " << sumArea << " +/- " << sumError << endl;
    cout << "Area  (net):              " << netArea << " +/- " << netAreaError << endl;
    cout << "+-------------------------------------------------+" << endl;
    cout << endl;
    cout << endl;

    // Draw results
    h1->Draw();
    fLin->SetLineColor(9);
    fLin->Draw("same");

    // Store and return peak info
    auto result = fitResults({netArea, centroid, centroidError, totalCentroidError*3});
    return result;
}

//
// In honor of Laura:
// ...cheese and puffy crust...
// Cheese and puffy crust
// Cheese and Puffy Crust!
// CHEESE AND PUFFY CRUST!
//
//



// Actual summing function
// This function doesn't fit the peak, only sums the bin contents
// Returns results (area, uncertainty) as a struct (so you can get both values from the one run)
fitResults sumPoly(TH1D* h1, double xMin1, double xMax1, double xMin2, double xMax2, double xMin3, double xMax3){

	// Define functions and set initial parameters
	double pi = 3.14159265359;
	//TH1F* h1 = (TH1F*)inFile->Get(hist);

	// Create the Canvas
    //TCanvas *cPoly = new TCanvas("cPoly", "Peak Summing - Cubic", 700, 500);
	h1->GetXaxis()->SetRangeUser(xMin1-100,xMax3+100);


    // Bounds for summing background/peak
   	Int_t xMinBin1 = h1->FindBin(xMin1);
   	Int_t xMaxBin1 = h1->FindBin(xMax1);
    Double_t xMin1Boundary = h1->GetBinLowEdge(xMinBin1);
    Double_t xMax1Boundary = h1->GetBinLowEdge(xMaxBin1+1);
   	Int_t xMinBin2 = h1->FindBin(xMin2);
   	Int_t xMaxBin2 = h1->FindBin(xMax2);
    Double_t xMin2Boundary = h1->GetBinLowEdge(xMinBin2);
    Double_t xMax2Boundary = h1->GetBinLowEdge(xMaxBin2+1);
   	Int_t xMinBin3 = h1->FindBin(xMin3);
   	Int_t xMaxBin3 = h1->FindBin(xMax3);
    Double_t xMin3Boundary = h1->GetBinLowEdge(xMinBin3);
    Double_t xMax3Boundary = h1->GetBinLowEdge(xMaxBin3+1);

    // Sum bins within boundary
    Double_t binWidth = h1->GetXaxis()->GetBinWidth(xMinBin1);
    Double_t sumPeak = 0;
    Double_t sumError = 0;
    Double_t sumCenter = 0;
    Double_t sumVariance = 0;
    Double_t sumSigma = 0;
    Double_t sumFWHM = 0;
    Double_t centroidError = 0;
    Double_t totalCentroidError = 0;
    Int_t binCounter1 = xMinBin1;   // BG1
    Int_t binCounter2 = xMinBin2;   // Peak
    Int_t binCounter3 = xMinBin3;   // BG2

    // Peak - Counter 2
    do {
            sumPeak += h1->GetBinContent(binCounter2);
            sumCenter += h1->GetBinContent(binCounter2)*h1->GetBinCenter(binCounter2);
            binCounter2++;
    }   while(binCounter2 <= xMaxBin2);

    // Using bin contents find the center of mass (centroid)
    sumCenter = sumCenter*binWidth;
    sumPeak = sumPeak*binWidth;
    sumError = sqrt(sumPeak);
    sumCenter = sumCenter/sumPeak;

    // Variance and sigma for peak contents
    binCounter2 = xMinBin2;
    do{
            sumVariance += h1->GetBinContent(binCounter2)*TMath::Power((h1->GetBinCenter(binCounter2)-sumCenter),2);
            centroidError += (TMath::Power(h1->GetBinError(binCounter2), 2) * TMath::Power((h1->GetBinCenter(binCounter2) - sumCenter), 2));
            binCounter2++;
    }   while(binCounter2 <= xMaxBin2);

    sumVariance = sumVariance*binWidth;
    sumVariance = sumVariance/sumPeak;
    sumSigma = sqrt(sumVariance);
    centroidError = centroidError*binWidth;
    centroidError = sqrt(centroidError);
    centroidError = centroidError/sumPeak;
    sumFWHM = sumSigma*2.35482;

    // Background fitting
    xExcludeMin = xMax1Boundary;
    xExcludeMax = xMin3Boundary;
    TF1 *fBG = new TF1("fBG", fLine, xMin1Boundary, xMax3Boundary, 4);
    fBG->SetParameters(1000, -0.005, 0.000001, 0.00000001);
    //fBG->SetParLimits(0, 0, 10000);
    //fBG->SetParLimits(1, 0, 100);
    //fBG->SetParLimits(2, 0, 100);
    //fBG->SetParLimits(3, 0, 100);
    reject = kTRUE;
    for (int i = 0; i < 10; i++) {
        h1->Fit(fBG, "LQM0", "", xMin1Boundary, xMax3Boundary);       // Fit 10 times, improving on the fit everytime
    }
    h1->Fit(fBG, "0Q");
    reject = kFALSE;

    Double_t a0 = fBG->GetParameter(0);
    Double_t a1 = fBG->GetParameter(1);
    Double_t a2 = fBG->GetParameter(2);
    Double_t a3 = fBG->GetParameter(3);

    // Integrating fit results
    Double_t bkgdIntegral = fBG->Integral(xMin2Boundary, xMax2Boundary);
    Double_t bkgdError = sqrt(bkgdIntegral);

    if(bkgdIntegral < 0){
        bkgdIntegral = 0;
        bkgdError = 0;
    }


    // Sum-Background 
    Double_t netArea;
    if((sumPeak - bkgdIntegral) < 0){
        netArea = sumPeak;
    }
    else {
        netArea = sumPeak - bkgdIntegral;
    }
    Double_t netError = TMath::Sqrt( TMath::Power(sumError,2.0) + TMath::Power(bkgdError,2.0) );

    // Total centroid error:
    // Contribution from background and weighting of bin contents center of mass bs
    // unc = sigma(centroid) * diff(Area)/Area
    Double_t diffFactor = TMath::Sqrt(sumPeak + bkgdIntegral);
    //cout << "DEBUG: diffFactor (sqrt) = " << diffFactor << endl;
    diffFactor = diffFactor / netArea;
    totalCentroidError = sumSigma * diffFactor;

    //cout << "DEBUG: diffFactor (division) = " << diffFactor << endl;
    cout << endl;

    totalCentroidError = TMath::Sqrt(totalCentroidError*totalCentroidError + centroidError*centroidError);

	// Print results
	// One day I'll write out to textfile or csv.
	cout << endl;
	cout << endl;
    cout << "+-------------------------------------------------+" << endl;
    cout << "Centroid:                 " << sumCenter << endl;
    cout << "Sigma:                    " << sumSigma <<  endl;
    cout << "FWHM (bins):              " << sumFWHM << endl;
    cout << "+-------------------------------------------------+" << endl;
    cout << "a0:                       " << a0 << endl;
    cout << "a1:                       " << a1 <<  endl;
    cout << "a2:                       " << a2 << endl;
    cout << "a3:                       " << a3 << endl;
    cout << "+-------------------------------------------------+" << endl;
    cout << "Area (bkgd):              " << bkgdIntegral << " +/- " << bkgdError << endl;
    cout << "Area (peak):              " << sumPeak << " +/- " << sumError << endl;
    cout << "Area  (net):              " << netArea << " +/- " << netError << endl;
    cout << "+-------------------------------------------------+" << endl;
    cout << endl;
    cout << endl;
    

    // Draw results
    h1->Draw();
    fBG->SetLineWidth(4);
    TH1 *hleft = (TH1*)h1->Clone("hleft");
    hleft->GetXaxis()->SetRange(xMinBin1,xMaxBin1);
    hleft->SetFillColor(40);
    hleft->Draw("hist same");
    TH1 *hright = (TH1*)h1->Clone("hright");
    hright->GetXaxis()->SetRange(xMinBin3,xMaxBin3);
    hright->SetFillColor(40);
    hright->Draw("hist same");
    TH1 *hsig = (TH1*)h1->Clone("hsig");
    hsig->GetXaxis()->SetRange(xMinBin2,xMaxBin2);
    hsig->SetFillColor(42);
    hsig->Draw("hist same");
    fBG->SetLineColor(9);
    fBG->Draw("same");

    auto result = fitResults({netArea, sumCenter, totalCentroidError, totalCentroidError*3});
    return result;

}




// Run method
// inputs:  runNum - the run number
//			f - file containing histograms
// 			counts - array for peak area
//			uncertainties - array for peak area uncertainties		   
// Function takes the input information and fits the histogram from runNum using the sumPoly method. 
// 		It takes both the fit area and uncertainty and subsequently adds them to the corresponding position 
// 		in the arrays counts and uncertainties, respectively. 
void Run(int runNum, TH1D* h0, double cent, double counts[31], double centroids[31], double uncertainties[31], double uncertainties2[31]){
    // CHANGE THE ARRAY SIZES FOR YOUR NUMBER OF HISTOGRAMS/FILES


	// Debug:
    //cout << "DEBUG: Will attempt to fit peak at " << cent << " keV for run number: " << runNum << "." << endl;
	//h0->Draw();

	// Fit peak and get results
    // Use the specified centroid location to define regions
    fitResults peak = fitResults();
    if (cent > 3000 and cent < 3100){
        // 3040 keV
        peak = fitSkewGauss(h0, cent, 1, int(cent)-40, int(cent)+40);
        //peak = sumPoly(h0, cent - 100, cent - 40, cent - 25, cent + 25, cent + 40, cent + 110);
    }
    else if (cent > 5100 and cent < 5220){
        // 5180 keV
        peak = fitSkewGauss(h0, cent, 1, int(cent)-40, int(cent)+40);
        //peak = sumPoly(h0, cent - 150, cent - 100, cent - 15, cent + 15, cent + 70, cent + 120);
    }
    else if (cent > 5220 and cent < 5300){
        // 5240 keV
        peak = fitSkewGauss(h0, cent, 1, int(cent)-40, int(cent)+40);
        //peak = sumPoly(h0, cent - 200, cent - 150, cent - 8, cent + 8, cent + 20, cent + 90);
    }
    else if (cent > 6090 and cent < 6150){
        // 6129 keV - Fluorine contamination
        peak = fitSkewGauss(h0, cent, 1, int(cent)-40, int(cent)+40);
        //peak = sumPoly(h0, cent - 90, cent - 25, cent - 10, cent + 10, cent + 70, cent + 150);
    }
    else if (cent > 6150 and cent < 6210){
        // 6180 keV
        peak = fitSkewGauss(h0, cent, 1, int(cent)-20, int(cent)+40);
        //peak = sumPoly(h0, cent - 150, cent - 100, cent - 10, cent +  8, cent + 40, cent + 90);
    }
    else if (cent > 6340 and cent < 6230){
        // 6180 keV
        peak = fitSkewGauss(h0, cent, 1, int(cent)-40, int(cent)+40);
        //peak = sumPoly(h0, cent - 70, cent - 25, cent - 10, cent + 10, cent + 40, cent + 90);
    }
    else if (cent > 6750 and cent < 6830){
        // 6793 keV
        peak = fitSkewGauss(h0, cent, 1, int(cent)-40, int(cent)+40);
        cout << "Fitting peak at " << cent << " keV." << endl;
        //peak = sumPoly(h0, cent - 90, cent - 25, cent - 10, cent + 10, cent + 20, cent + 50);
    }
    else {
        cout << "No peak found." << endl;
        peak = fitSkewGauss(h0, cent, 1, int(cent)-40, int(cent)+40);
        //peak = sumPoly(h0, cent - 90, cent - 25, cent - 10, cent + 10, cent + 20, cent + 50);
    }

    // Save peak info into arrays (as passed into the function) by modifying 
    //     the array location corresponding to the run number
    //cout << "DEBUG: Using a struct for the peak info. Counts, unc: " << peak.counts << "," << peak.unc << endl;
    counts[runNum] = peak.cts;
    centroids[runNum] = peak.centroid;
    uncertainties[runNum] = peak.unc;
    uncertainties2[runNum] = peak.unc2;

}




// ---------------
// THIS IS FUCKING UGLY
// But i didn't know a better way other than hard coding in the histogram 
// and it's corresponding title.
// I Saved all of the histograms from their root files into a single file
// And then read that file in here

// getHistograms method
// inputs: f - the file containing histograms
//         h[30] - the array in which to fill the histograms
//               - This needs to change depending on the number of files/histograms
// Function fills the histogram array with each histogram in the file
void getHistograms(TFile *f, TH1D *h[31]){
    // CHANGE THE h ARRAY SIZE TO MATCH YOUR NUMBER OF HISTOGRAMS

    TKey *key = f->FindKey("hist_Ta30_0deg;2");
    if (key == 0){
        cout << "!!Histogram 0 does not exist!!" << endl;
        throw 1;
    }
    h[0] =  (TH1D*)f->Get("hist_Ta30_0deg;2");
    key = f->FindKey("hist_Ta30_0deg;2");
    if (key == 0){
        cout << "!!Histogram 1 does not exist!!" << endl;
        throw 1;
    }
    h[1] =  (TH1D*)f->Get("hist_Ta30_0deg;2");
    key = f->FindKey("hist_Ta30_0deg;2");
    if (key == 0){
        cout << "!!Histogram 2 does not exist!!" << endl;
        throw 1;
    }
    h[2] =  (TH1D*)f->Get("hist_Ta30_0deg;2");
    key = f->FindKey("hist_Ta30_0deg;2");
    if (key == 0){
        cout << "!!Histogram 3 does not exist!!" << endl;
        throw 1;
    }
    h[3] =  (TH1D*)f->Get("hist_Ta30_0deg;2");
    key = f->FindKey("hist_Ta30_0deg;2");
    if (key == 0){
        cout << "!!Histogram 4 does not exist!!" << endl;
        throw 1;
    }
    h[4] =  (TH1D*)f->Get("hist_Ta30_0deg;2");
    key = f->FindKey("hist_Ta30_0deg;2");
    if (key == 0){
        cout << "!!Histogram 5 does not exist!!" << endl;
        throw 1;
    }
    h[5] =  (TH1D*)f->Get("hist_Ta30_0deg;2");
    key = f->FindKey("hist_Ta30_0deg;2");
    if (key == 0){
        cout << "!!Histogram 6 does not exist!!" << endl;
        throw 1;
    }
    h[6] =  (TH1D*)f->Get("hist_Ta30_0deg;2");
    key = f->FindKey("hist_Ta30_0deg;2");
    if (key == 0){
        cout << "!!Histogram 7 does not exist!!" << endl;
        throw 1;
    }
    h[7] =  (TH1D*)f->Get("hist_Ta30_0deg;2");
    key = f->FindKey("hist_Ta30_0deg;2");
    if (key == 0){
        cout << "!!Histogram 8 does not exist!!" << endl;
        throw 1;
    }
    h[8] =  (TH1D*)f->Get("hist_Ta30_0deg;2");
    key = f->FindKey("hist_Ta30_45deg;2");
    if (key == 0){
        cout << "!!Histogram 9 does not exist!!" << endl;
        throw 1;
    }
    h[9] =  (TH1D*)f->Get("hist_Ta30_45deg;2");
    key = f->FindKey("hist_Ta30_60deg;2");
    if (key == 0){
        cout << "!!Histogram 10 does not exist!!" << endl;
        throw 1;
    }
    h[10] =  (TH1D*)f->Get("hist_Ta30_60deg;2");
    key = f->FindKey("hist_Ta30_75deg;2");
    if (key == 0){
        cout << "!!Histogram 11 does not exist!!" << endl;
        throw 1;
    }
    h[11] =  (TH1D*)f->Get("hist_Ta30_75deg;2");
    key = f->FindKey("hist_Ta30_90deg;2");
    if (key == 0){
        cout << "!!Histogram 12 does not exist!!" << endl;
        throw 1;
    }
    h[12] =  (TH1D*)f->Get("hist_Ta30_90deg;2");
    key = f->FindKey("hist_Ta30_111deg;2");
    if (key == 0){
        cout << "!!Histogram 13 does not exist!!" << endl;
        throw 1;
    }
    h[13] =  (TH1D*)f->Get("hist_Ta30_111deg;2");
    key = f->FindKey("hist_Ta30_135deg;2");
    if (key == 0){
        cout << "!!Histogram 14 does not exist!!" << endl;
        throw 1;
    }
    h[14] =  (TH1D*)f->Get("hist_Ta30_135deg;2");
    key = f->FindKey("hist_Ta30_75deg;3");
    if (key == 0){
        cout << "!!Histogram 15 does not exist!!" << endl;
        throw 1;
    }
    h[15] =  (TH1D*)f->Get("hist_Ta30_75deg;3"); // Actually 60 degrees, mislabeled
    key = f->FindKey("hist_ZrN_0deg;1");
    if (key == 0){
        cout << "!!Histogram 16 does not exist!!" << endl;
        throw 1;
    }
    h[16] =  (TH1D*)f->Get("hist_ZrN_0deg;1");
    key = f->FindKey("hist_ZrN_90deg;1");
    if (key == 0){
        cout << "!!Histogram 17 does not exist!!" << endl;
        throw 1;
    }
    h[17] =  (TH1D*)f->Get("hist_ZrN_90deg;1");
    key = f->FindKey("hist_Mo30_0deg;1");
    if (key == 0){
        cout << "!!Histogram 18 does not exist!!" << endl;
        throw 1;
    }
    h[18] =  (TH1D*)f->Get("hist_Mo30_0deg;1");
    key = f->FindKey("hist_Mo30_45deg;1");
    if (key == 0){
        cout << "!!Histogram 19 does not exist!!" << endl;
        throw 1;
    }
    h[19] =  (TH1D*)f->Get("hist_Mo30_45deg;1");
    key = f->FindKey("hist_Mo30_60deg;1");
    if (key == 0){
        cout << "!!Histogram 20 does not exist!!" << endl;
        throw 1;
    }
    h[20] =  (TH1D*)f->Get("hist_Mo30_60deg;1");
    key = f->FindKey("hist_Mo30_75deg;1");
    if (key == 0){
        cout << "!!Histogram 21 does not exist!!" << endl;
        throw 1;
    }
    h[21] =  (TH1D*)f->Get("hist_Mo30_75deg;1");
    key = f->FindKey("hist_Mo30_90deg;1");
    if (key == 0){
        cout << "!!Histogram 22 does not exist!!" << endl;
        throw 1;
    }
    h[22] =  (TH1D*)f->Get("hist_Mo30_90deg;1");
    key = f->FindKey("hist_Mo30_111deg;1");
    if (key == 0){
        cout << "!!Histogram 23 does not exist!!" << endl;
        throw 1;
    }
    h[23] =  (TH1D*)f->Get("hist_Mo30_111deg;1");
    key = f->FindKey("hist_Mo30_135deg;1");
    if (key == 0){
        cout << "!!Histogram 24 does not exist!!" << endl;
        throw 1;
    }
    h[24] =  (TH1D*)f->Get("hist_Mo30_135deg;1");
    key = f->FindKey("hist_W30_135deg;1");
    if (key == 0){
        cout << "!!Histogram 25 does not exist!!" << endl;
        throw 1;
    }
    h[25] =  (TH1D*)f->Get("hist_W30_135deg;1");
    key = f->FindKey("hist_W30_90deg;1");
    if (key == 0){
        cout << "!!Histogram 26 does not exist!!" << endl;
        throw 1;
    }
    h[26] =  (TH1D*)f->Get("hist_W30_90deg;1");
    key = f->FindKey("hist_W30_0deg;1");
    if (key == 0){
        cout << "!!Histogram 27 does not exist!!" << endl;
        throw 1;
    }
    h[27] =  (TH1D*)f->Get("hist_W30_0deg;1");
    key = f->FindKey("hist_W30_45deg;1");
    if (key == 0){
        cout << "!!Histogram 28 does not exist!!" << endl;
        throw 1;
    }
    h[28] =  (TH1D*)f->Get("hist_W30_45deg;1");
    key = f->FindKey("hist_Ta10_0deg;1");
    if (key == 0){
        cout << "!!Histogram 29 does not exist!!" << endl;
        throw 1;
    }
    h[29] =  (TH1D*)f->Get("hist_Ta10_0deg;1");
    key = f->FindKey("hist_Ta10_90deg;1");
    if (key == 0){
        cout << "!!Histogram 30 does not exist!!" << endl;
        throw 1;
    }
    h[30] =  (TH1D*)f->Get("hist_Ta10_90deg;1");
    
    

    cout << endl;
    cout << "Histograms retrieved. " << endl;
    cout << endl;
}



// fitSpectra method	   
// Function is the main method for this program. It opens the histogram file and calls the run function 
// 		to fit each histogram in the file. After each fit it saves an image of the canvas to file so it 
//		can be verified later. It finally takes the output from these fits (area and uncertainty) and 
//		writes them to a text file so that they can be imported into excel. 
void fitSpectra(int energy){

	// Arrays for fit data
    // CHANGE THE SIZE FOR YOUR SITUATION
    double counts[31] = {};
    double centroids[31] = {};
    double uncertainties[31] = {};
    double uncertainties2[31] = {};
    double cents[31] = {};

	// Open file of the summed spectra
    // YOUR FILENAME HERE CONTAINING ALL HISTOGRAMS
	TFile *f = new TFile("./summedSpectra.root");

    // Get all the histograms
    TH1D *h[31];
    getHistograms(f, h);
    
    
    for(int j = 0; j < 31; j++){
        cout << "DEBUG: Histograms - " << j << " - " << h[j]->GetTitle() << endl;
    }

    cout << endl;
    cout << "Histograms retrieved. " << endl;
    cout << endl;

    // Debug
    //h[0]->Draw();

	// Create canvases
    // Works on two canvases so that output images aren't messy
    TCanvas *canvasFinder = new TCanvas("finder", "finder", 1400, 1000);
	TCanvas *canvas = new TCanvas("canvas", "14N(p,g) Lifetime", 1400, 1000);

    char canvasName[40];
    char canvasNameOpen[40];
    char canvasNameClose[40];
    sprintf(canvasNameOpen, "./output/FITsummedSpectra_%04d.pdf[", energy);
    sprintf(canvasName, "./output/FITsummedSpectra_%04d.pdf", energy);
    sprintf(canvasNameClose, "./output/FITsummedSpectra_%04d.pdf]", energy);

    // Opens output pdf file
	canvas->Print(canvasNameOpen);

    // Debug
    /*
    for(int i = 0; i < 31; i++){
        h[i]->Draw();
        canvas->Print(canvasName);
        cout << "DEBUG: Print canvas h[" << i << "]" << endl;
    }
    */

	// Find peak
    //double cent = 0;
    canvasFinder->cd();
    for (int i = 0; i < 31; i++){
        cout << endl;
        cout << "Run Number " << i << endl;
        cents[i] = PeakFinder(i, h[i], energy); 
        if (cents[i] == 0) {
            cents[i] = energy;
        }
    }

    // Rerun and save data/image
    canvas->cd();
	for(int i = 0; i < 31; i++){
        cout << endl;
        cout << "Run Number:\t" << i << endl;
        //cent = PeakFinder(i, f, energy);
        //cout << "Centroid:\t" << cents[i] << endl;
		Run(i, h[i], cents[i], counts, centroids, uncertainties, uncertainties2);
		canvas->Print(canvasName);
		// Debug:
		//cout << "DEBUG: Counts = " << counts[0] << "\t Uncertainty = " << uncertainties[0] << endl;
	}

	// Close canvas
	canvas->Print(canvasNameClose);

	// Close histogram file
	f->Close();


	// Output data to text file
    char fileName[40];
    sprintf(fileName, "./output/fitData_%04d.csv", energy);
    ofstream fOut(fileName, ios::app);
    for(int i = 0; i < 31; i++){
    	fOut << counts[i] << "," << centroids[i] << "," << uncertainties[i] << "," << uncertainties2[i] << endl;
    }

    fOut.close();

    cout << "File finished. " << endl;
}

int main(){

    int energy = 6793;
    fitSpectra(energy);

    return 0;
}
