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
    double binArray[401];

    //rebin histograms so all bins have statistics
    for(int j =0; j<401; j++){
        binArray[j]= ((pow(2.7162,j*0.02+2) - 1)/10000);
    }


    for(int i=0; i<size; i++)
    {
        //create histogram
        hist.push_back(new TH1F(Form("hist%d",i),Form("VoV = %f",VoV[i]), 400, binArray));
        //open file and TTreeReader
        if(gSystem->AccessPathName(filenames[i]))
        {
            std::cout << "file does not exist" << std::endl;
        } 
        else 
        {
            std::cout << "file exists" << std::endl;
        };
       
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

                //logic to deal with two-seven pulses back to back (rarely see more than 2 except at room temp)
                //fill the hist with time between of one pulse width (50 ns)
                //increase the current time
                if(myWIDTH[j]>40.0)
                {
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    timeNow += 25*(2.0/1000000000.0);
                }
                else if(myWIDTH[j]>55)
                {
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    timeNow += 50*(2.0/1000000000.0);
                }
                else if(myWIDTH[j]>80)
                {
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    timeNow += 75*(2.0/1000000000.0);
                }
                else if(myWIDTH[j]>100)
                {
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    timeNow += 100*(2.0/1000000000.0);
                }
                else if(myWIDTH[j]>130)
                {
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    timeNow += 125*(2.0/1000000000.0);
                }
                else if(myWIDTH[j]>150)
                {
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    timeNow += 150*(2.0/1000000000.0);
                }
                else if(myWIDTH[j]>180)
                {
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    hist[i]->Fill(25*(2.0/1000000000.0));
                    timeNow += 175*(2.0/1000000000.0);
                };
                

                timeBefore = timeNow;
                trigLast = myTrigTIME[j];
                
            }
            
        }   
        //get error per each bin and then scale each bin by its width to normalize
        hist[i]->Sumw2();
        hist[i]->Scale(1./hist[i]->Integral(),"width");

        hist[i]->GetYaxis()->SetTitle("Probability");
        hist[i]->GetXaxis()->SetTitle("Time After Primary Pulse [s]");

        delete myFile;

    }
    return hist;
}

void PDC_DarkAnalysis(){

    //store the fits in a vector to retrieve later for error propagation
    std::vector< TF1* > fit;

    const int size = 7; // will need to change number of filenames in files array to match (should do as a vector, but root was seg faulting for me)

    int temp = 20;

    float fitRate[size];
    float fitError[size];



    int dataNumbers[size] = {528,529,530,531,533,534,535};

    //40 {463,464,465,466,467,468};
    //60 {450,451,452,453,454,455};
    //80 {443,444,445,446,447,448};
    //100 {436,437,438,439,440,441};
    //120 {429,430,431,432,433,434};
    //140 {422,423,424,425,426,427};
    //160 {415,416,417,418,419,420};

    //room temp{528,529,530,531,533,534,535};

    float overVoltages[size] = {0.5,1.5,3,4.5,6,7,8};
    //old data {0.1,0.2,0.7,1.2,1.7,2.7};
    //room temp{0.5,1.5,3,4.5,6,7,8};

    float fitEndRange[size] = {2,1.6,1,.5,.5};

    const char *files[size] = {Form("root_output_files/output00%d.root",dataNumbers[0]), Form("root_output_files/output00%d.root",dataNumbers[1]),Form("root_output_files/output00%d.root",dataNumbers[2]),Form("root_output_files/output00%d.root",dataNumbers[3]),Form("root_output_files/output00%d.root",dataNumbers[4]),Form("root_output_files/output00%d.root",dataNumbers[5]),Form("root_output_files/output00%d.root",dataNumbers[6])};

    float overVolErrors[size] = {0.05,0.05,0.05,0.05,0.05,0.05,0.05};

    // initialize files
    
    //Form("root_output_files/output00%d.root",dataNumbers[4]),Form("root_output_files/output00%d.root",dataNumbers[5])
    vector<TH1F*> histograms = histCollector(overVoltages, files, temp, size);
    // output file
    TFile *out = new TFile(Form("histOutput%d.root",temp), "RECREATE");

    for (int count = 0; count < size; count ++){

        //fit range below is a workaround for the strange gaussian behavior of my time difference plots... don't understand the underlying distribution (PROBLEM!)

        histograms[count]->Fit("expo","WL","",0.0001,fitEndRange[count]); // L specifies log likelihood (which deals with the non-gaussian bin statistics in low count bins). We only fit after the first several bins to ignore afterpulsing. 

        fit.push_back(histograms[count]->GetFunction("expo")); // save the fit paramters to a vector

        //write histo and fit
        histograms[count]->Write(Form("Hist%d",count));
        fit[count]->Write(Form("Fit%d",count));

        fitRate[count] = (fit[count]->GetParameter(1))*(-1000000000)/(1.296);
        fitError[count] = fit[count]->GetParError(1)*(1000000000)/(1.296);
    }
    // make the graphs
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