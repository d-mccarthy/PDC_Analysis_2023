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
#include "TGraphErrors.h"
#include "TH2D.h"
#include "THStack.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

using namespace std;


vector<TH1F*> histCollector(float VoV[], const char *filenames[], int temp, int size)
{
    std::vector< TH1F* > hist;

    double windowResolution = 2.0/1000000000.0; //time resolution within window is 2ns
    double triggerResolution = 16.0/1000000000.0; //time resolution within full run is 16ns

    for(int i=0; i<size; i++)
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

    std::vector< TF1* > fit;


    int size = 6;
    int temp = 40;

    float fitRate[6];
    float fitError[6];

    int dataNumbers[6] = {463,464,465,466,467,468};
    //40 {463,464,465,466,467,468};
    //60 {450,451,452,453,454,455};
    //80 {443,444,445,446,447,448};
    //100 {436,437,438,439,440,441};
    //120 {429,430,431,432,433,434};
    //140 {422,423,424,425,426,427};
    //160 {415,416,417,418,419,420};




    float overVoltages[6] = {0.2,0.3,0.7,1.2,1.7,2.7};
    float overVolErrors[6] = {0.05,0.05,0.05,0.05,0.05,0.05};
    const char *files[6] = {Form("root_output_files/output00%d.root",dataNumbers[0]), Form("root_output_files/output00%d.root",dataNumbers[1]),Form("root_output_files/output00%d.root",dataNumbers[2]),Form("root_output_files/output00%d.root",dataNumbers[3]),Form("root_output_files/output00%d.root",dataNumbers[4]),Form("root_output_files/output00%d.root",dataNumbers[5])};
    
    vector<TH1F*> histograms = histCollector(overVoltages, files, temp, size);
    
    TFile *out = new TFile(Form("histOutput%d.root",temp), "RECREATE");

    for (int count = 0; count < 6; count ++){
        histograms[count]->Fit("expo","","",0.004,2);
        fit.push_back(histograms[count]->GetFunction("expo"));

        //write histo and fit
        histograms[count]->Write(Form("Hist%d",count));
        fit[count]->Write(Form("Fit%d",count));

        fitRate[count] = (fit[count]->GetParameter(1))*(-1000000000)/(1.296);
        fitError[count] = fit[count]->GetParError(1)*(1000000000)/(1.296);
    }
    auto c1 = new TCanvas("c1","VoV vs Slope Fit",200,10,700,500);
    c1->SetFillColor(0);
    c1->SetGridx();
    c1->SetGridy();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);

    auto gr = new TGraphErrors(size, overVoltages,fitRate,overVolErrors,fitError);
    gr->SetMarkerStyle(22);
    gr->GetXaxis()->SetTitle("VoV [V]");
    gr->GetYaxis()->SetTitle("DCR [Hz/m^{2}]");
    gr->Draw();
    gr->Write(Form("Graph%d",temp));
    
}