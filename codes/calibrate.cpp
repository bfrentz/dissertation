/*******************************************
 * The calibration program from the 
 * 14N(p,g)15O 6.79 MeV state lifetime
 * measurement (data replay and calibration)
 * Usage with root:
 *       .L calibrate.cpp
 *       Run(runNumber) - for individal run
 *       addData() - for all runs (need to change loop)
 *
 * The way I wrote this was to first "RECRETE" the calibratedSpectra.root file
 * so I interactively do Run(1), then I change this file so that the run() function
 * now updates the TFile you create. After that point, simply load into root and use
 * the runAll() function, after changing your run number bounds, to create and add the
 * remaining histograms to the new storage file.
 *
 * Created 8 Aug 2019
 *
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

using namespace std;


// Run function
// Inputs: integer runNum
// Outputs: 
// This method takes a run number and takes the data from the given run, calibrates the spectra, 
// 		and outputs the resulting histogram to a new file (calibratedSpectra.root). 
void Run(int runNum){

	// Create a random # generator in range [0,1]
	TRandom2 myRandom(0);


    //Optimize cache and disk reading because this isn't default
    //in this version of ROOT.
    TChain chain("evtTree");
    chain.SetCacheSize(1E8);
    //chain.AddBranchToCache("*");

	// For adding files in chain
    ifstream fRuns;
    char dataf[255];


    // Get data
    sprintf(dataf,"./lifetime_%03d.root", runNum);
    TFile *f = new TFile(dataf);
    cout << "Data File: " << dataf << endl;
	chain.Add(dataf);
    

	// Create an output file for the histograms
    TFile file("./calibratedSpectra.root","UPDATE");


    // Create histogram
    char title[255];
    sprintf(title, "Georgina Spectrum run%03d", runNum);
    TH1D *h0 = new TH1D("h0", title, 8192, -0.5, 9999.5);

    char hname[7];
    sprintf(hname, "run%03d", runNum);

    // Centroid of 1460.83: 8718.83
    // Centroid of 6129.59: 36609.6
    // Calibration slope and offset:
    Double_t a = 0.167393;
    Double_t b = 1.381436;
    // Event energy
    Double_t ecal;

    // Get event tree
    TTreeReader reader(&chain);
    TTreeReaderArray<UShort_t> event_e(reader,"event.e");

    ULong64_t nEvent = chain.GetEntries();
    cout << nEvent << " events to be analyzed." << endl;


    // Start the event loop
    while(reader.Next()) {
    	if(event_e[0] > 0) {
    		// Calibrate from channel to energy and fill output run histogram
    		ecal = a * (event_e[0] + myRandom.Rndm()-0.5) + b;
    		h0->Fill(ecal);
    	}
    	else ecal = 0;
    }

    h0->SetName(hname);
    h0->Write(hname, TObject::kWriteDelete);

    //file.Write();
    file.Close();

    cout << hname << " has been added." << endl;
    cout << endl;

    f->Close();

}


// Main method
// Inputs:
// Outputs:
// This method loops over all files, invoking the Run method in order to add all the data to one file
void addData(){

	for(int i = 192; i < 241; i++){
		Run(i);
	}

	cout << endl;
	cout << "Finished adding all runs." << endl;
	cout << "Have a nice day." << endl;
	cout << endl;

}
