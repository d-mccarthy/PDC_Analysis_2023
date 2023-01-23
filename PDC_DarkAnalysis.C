#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <typeinfo>
#include "TClass.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH2D.h"
#include "THStack.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

using namespace std;


void PDC_DarkAnalysis(string inDirectoryName){

    const char* fileName = inDirectoryName.c_str();

    filesystem::


    // Create a histogram for the values we read.
    TH1D *histTimeBetween = new TH1D("TimeBetween", "Time Difference at -100C and 2.7 VOV", 500, 0.00, 1);

    // Open the file containing the tree.
    TFile *myFile = TFile::Open(fileName);
    TTreeReader myReader("TPulse", myFile);
    TTreeReaderArray<Double_t> myMTS(myReader, "mts");
    TTreeReaderArray<Double_t> myTIME(myReader, "time");
    TTreeReaderArray<Double_t> myTrigTIME(myReader,"trigtime");
    TTreeReaderArray<Double_t> myWIDTH(myReader,"ToT");

    
    bool first = true;

    //counter to deal with trigger time steps (see DAQ manual for details)
    int trigCounter = 0;
    double trigLast = 0;

    //logic variables for creating histograms
    double t0;
    double timeBefore = 0;
    double timeNow;
    double timeBetween;

    double windowResolution = 2.0/1000000000.0; //time resolution within window is 2ns
    double triggerResolution = 16.0/1000000000.0; //time resolution within full run is 16ns

    // Loop over all entries of the TTree.
    while (myReader.Next()) {
        if (first){
            //cout<<(myTrigTIME[0]+trigCounter*2147483648)*(triggerResolution) + myTIME[0] * (windowResolution)<<endl;
            first = false;
        }
        cout.precision(15);
        
        for (int j = 0; j<myTIME.GetSize(); j++){
            if(myTrigTIME[j] < trigLast){
                trigCounter++;
            }
            
            timeNow = (myTrigTIME[j]+trigCounter*2147483648)*(triggerResolution) + myTIME[j] * (windowResolution);
            timeBetween = timeNow - timeBefore;
            
            histTimeBetween->Fill(timeBetween);

            //logic to deal with two pulses back to back
            //fill the hist with time between of one pulse width (50 ns)
            //increase the current time
            if(myWIDTH[j]>40.0){
                histTimeBetween->Fill(25*(2.0/1000000000.0));
                timeNow += 25*(2.0/1000000000.0);
            }

            timeBefore = timeNow;
            trigLast = myTrigTIME[j];
        }
        
    }
    
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    histTimeBetween->Draw();
}