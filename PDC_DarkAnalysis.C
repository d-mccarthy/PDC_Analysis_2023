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


vector<TH1F*> histCollector(float VoV[7], const char *filenames[7], int temp)
{
    std::vector< TH1F* > hist;

    double windowResolution = 2.0/1000000000.0; //time resolution within window is 2ns
    double triggerResolution = 16.0/1000000000.0; //time resolution within full run is 16ns

    for(int i=0; i<7; i++)
    {
        //create histogram
        hist.push_back(new TH1F(Form("hist%d",i),Form("VoV = %f",VoV[i]), 500, 0, 2));
        //open file and TTreeReader
       
        TFile *myFile = TFile::Open(filenames[i]);
        TTreeReader myReader("TPulse", myFile);
        TTreeReaderArray<Double_t> myMTS(myReader, "mts");
        TTreeReaderArray<Double_t> myTIME(myReader, "time");
        TTreeReaderArray<Double_t> myTrigTIME(myReader,"trigtime");
        TTreeReaderArray<Double_t> myWIDTH(myReader,"ToT");

        //counter to deal with trigger time steps (see DAQ manual for details)
        int trigCounter = 0;
        double trigLast = 0;

        //logic variables for creating histograms
        double t0;
        double timeBefore = 0;
        double timeNow;
        double timeBetween;

        while (myReader.Next()) 
        {
        
            for (int j = 0; j<myTIME.GetSize(); j++)
            {
                if(myTrigTIME[j] < trigLast)
                {
                    trigCounter++;
                }
                
                timeNow = (myTrigTIME[j]+trigCounter*2147483648)*(triggerResolution) + myTIME[j] * (windowResolution);
                timeBetween = timeNow - timeBefore;
                
                hist[i]->Fill(timeBetween);

                //logic to deal with two pulses back to back
                //fill the hist with time between of one pulse width (50 ns)
                //increase the current time
                if(myWIDTH[j]>40.0)
                {
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    timeNow += 25*(2.0/1000000000.0);
                }

                timeBefore = timeNow;
                trigLast = myTrigTIME[j];
            }
            
        }   

        delete myFile;

    }
    return hist;
}

void PDC_DarkAnalysis(){

    float overVoltages[7] = {0.1,0.2,0.7,1.2,1.7,2.7};
    const char *files[7] = {"root_output_files/output00414.root", "root_output_files/output00415.root","root_output_files/output00416.root","root_output_files/output00417.root","root_output_files/output00418.root","root_output_files/output00419.root","root_output_files/output00420.root"};
    vector<TH1F*> histograms = histCollector(overVoltages, files, 160);

    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    
    TFile *out = new TFile("histOutput160", "RECREATE");

    for (int count = 0; count < 7; count ++){
        histograms[count]->Write(Form("Hist%d",count));
    }
    
    
}