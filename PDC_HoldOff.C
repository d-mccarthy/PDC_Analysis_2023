#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <typeinfo>
#include "TRandom.h"
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

double windowResolution = 2.0/1E9; //time resolution within window is 2ns
double triggerResolution = 16.0/1E9; //time resolution within full run is 16ns

double binArray[121];
float pulseWidth[6] = {80.0,37.0,25.0,19.0,15.0,13.0}; //in ADC counts -- 50 ns pulses are 25 ADC
float pulseWidthNS[6] = {160.0,74.0,50.0,38.0,30.0,26.0};


void filler(int repeats, TH1F* toFill, float width){
    for(int inx = 0; inx<repeats;inx++){
        toFill->Fill(width*(2.0/1E9));
    }
    
}

// Function to give the afterpulsing distribution with the probability distribution function of the total DNR as input.
vector<TH1F*> unShadower(vector<TH1F*> inputHist){
    //define the binning of the histogram. This was done by eye to give reasonable statistics in each bin.
    for(int j =0; j<121; j++){

        binArray[j]= ((pow(2.7162,j*0.19))/1E9);    
        
    }
    // initialize a vector of histograms to give as output
    std::vector< TH1F* > unShadowed;
    // initialize variables: following from 
    double lambda;
    double beta;
    double sigPi;
    double sigLambda;
    double sigBeta;
    double Pi;
    double PiWidth;
    int numBins;

    for(int hInx = 0; hInx<inputHist.size(); hInx++){
        unShadowed.push_back(new TH1F(Form("unShadowedHist%d",hInx),Form("PulseWidth = %f [ns]",2*pulseWidth[hInx]),120,binArray));
        beta=0.;
        sigBeta=0.;
        numBins = inputHist[hInx]->GetNbinsX();


        for(int ti = 0; ti < numBins-5 ; ti++){

            Pi = inputHist[hInx]->GetBinContent(ti);

            sigPi = (inputHist[hInx]->GetBinError(ti));
           
            lambda = -log(1-(Pi/exp(-beta)));

            sigLambda = sqrt(sigPi*sigPi+Pi*Pi*sigBeta*sigBeta)/(exp(-beta)-Pi);
            cout<< lambda << " with error: " << sigLambda << endl;
            sigBeta = sqrt(sigPi*sigPi+exp(-2*beta)*sigBeta*sigBeta)/(exp(-beta)-Pi);
            
            beta += lambda;

            if(!isnan(lambda)){
                unShadowed[hInx]->SetBinContent(ti,lambda);
                unShadowed[hInx]->SetBinError(ti,sigLambda);
            }
            
        
            
        }
        float fitEnd[6] = {.1,.1,.1,.1,.1,.1};
        unShadowed[hInx]->Sumw2();
        unShadowed[hInx]->Scale(1./unShadowed[hInx]->Integral());  
        unShadowed[hInx]->Scale(1.,"width");
        unShadowed[hInx]->Fit("pol0","WL","",1E-2,fitEnd[5-hInx]);
    }

    return unShadowed;
        
}

vector<TH1F*> histCollector(float VoV[], const char *filenames[], int temp, int size)
{
    std::vector< TH1F* > hist;
    
    const int length = size;
    
    //rebin histograms so all bins have statistics
    for(int j =0; j<121; j++){

        binArray[j]= ((pow(2.7162,j*0.19))/1E9);  
        
    }

    for(int i=0; i<size; i++)
    {
        //create histogram
        hist.push_back(new TH1F(Form("hist%d",i),Form("PulseWidth = %f [ns]",2*pulseWidth[i]), 120, binArray));
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
                
                // if (i == 2 && timeBetween < 50/1E9){
                //     std::cout<<"Time is too short at event = "<< myEvent[j] <<" time is : "<< timeBetween <<std::endl;
                // }

                hist[i]->Fill(timeBetween); 

                //logic to deal with two->seven pulses back to back (rarely see more than 2 except at room temp)
                //fill the hist with time between of one pulse width (50 ns)
                //increase the current time

                if(myWIDTH[j]>(pulseWidth[i]+5))
                {
                    filler(1,hist[i],pulseWidth[i]);
                    timeNow += pulseWidth[i]*(2.0/1000000000.0);
                }
                else if(myWIDTH[j]>(2*pulseWidth[i]+5))
                {
                    filler(2,hist[i],pulseWidth[i]);
                    timeNow += 2*pulseWidth[i]*(2.0/1000000000.0);
                }
                else if(myWIDTH[j]>(3*pulseWidth[i]+5))
                {
                    filler(3,hist[i],pulseWidth[i]);
                    timeNow += 3*pulseWidth[i]*(2.0/1000000000.0);
                }
                else if(myWIDTH[j]>(4*pulseWidth[i]+5))
                {
                    filler(4,hist[i],pulseWidth[i]);
                    timeNow += 4*pulseWidth[i]*(2.0/1000000000.0);
                }
                else if(myWIDTH[j]>(5*pulseWidth[i]+5))
                {
                    filler(5,hist[i],pulseWidth[i]);
                    timeNow += 5*pulseWidth[i]*(2.0/1000000000.0);
                }
                else if(myWIDTH[j]>(6*pulseWidth[i]+5))
                {
                    filler(6,hist[i],pulseWidth[i]);
                    timeNow += 6*pulseWidth[i]*(2.0/1000000000.0);
                }
                else if(myWIDTH[j]>(7*pulseWidth[i]+5))
                {
                    filler(7,hist[i],pulseWidth[i]);
                    timeNow += 7*pulseWidth[i]*(2.0/1000000000.0);
                };
                
                timeBefore = timeNow;
                trigLast = myTrigTIME[j];
                
            }
            
        }   
        //get error per each bin and then scale each bin by its width to normalize
        hist[i]->Sumw2();
        
        hist[i]->Scale(1./hist[i]->Integral());     

        hist[i]->GetYaxis()->SetTitle("Probability/s");
        hist[i]->GetXaxis()->SetTitle("Time After Primary Pulse [s]");

        delete myFile;

    }
    return hist;
}

// make histograms for afterpulsing (subtract the fit, etc)


void PDC_HoldOff(){

    //store the fits in a vector to retrieve later for error propagation
    std::vector< TF1* > fit;
    std::vector< TF1* > fitUS;

    const int size = 6; // will need to change number of filenames in files array to match (should do as a vector, but root was seg faulting for me)

    int temp = 999;

    float fitRate[size];
    float fitError[size];
    float fitRateUS[size];
    float fitErrorUS[size];

    float temperature[1] = {999.0};

    int dataNumbers[size] = {602,603,599,604,605,606};


    //LOW STAT DATA
    //40 {463,464,465,466,467,468};
    //60 {450,451,452,453,454,455};
    //80 {443,444,445,446,447,448};
    //100 {436,437,438,439,440,441};
    //120 {429,430,431,432,433,434};
    //140 {422,423,424,425,426,427};
    //160 {415,416,417,418,419,420};

    //HIGH STAT DATA
    //300{528,529,530,531,533,534,535};
    //270{580,581,582,583,584,585,586,587};
    //260{564,565,566,567,568,569,570,571};
    //250{548,549,550,551,552,553,554,555};
    //240{556,557,558,559,560,561,562,563};
    //230{572,573,574,575,576,577,578,579};
    //-60 {539,540,541,542,543};
    //-80 {537,544,545,546,547};

    //hold off {589,590,591,592,593,594};
    //hold off -50 {602,603,599,604,605,606};


    float overVoltages[size] = {0.2,0.5,0.75,1,1.25,1.5};
    //old data {0.1,0.2,0.7,1.2,1.7,2.7};
    //new data 
    //room temp{0.5,1.5,3,4.5,6,7,8};
    //hold off{0.1,0.2,0.5,0.75,1,1.25,1.5,1.75}

    float fitEndRange[size] = {1,1,1,1,1,1};
    float fitStartRange[size] = {1E-3,1E-3,1E-3,1E-3,1E-3,1E-3};

    const char *files[size] = {Form("root_output_files/output00%d.root",dataNumbers[0]),
    Form("root_output_files/output00%d.root",dataNumbers[1]),
    Form("root_output_files/output00%d.root",dataNumbers[2]),
    Form("root_output_files/output00%d.root",dataNumbers[3]),
    Form("root_output_files/output00%d.root",dataNumbers[4]),
    Form("root_output_files/output00%d.root",dataNumbers[5])
    };

    float overVolErrors[size] = {1.0,1.0,1.0,1.0,1.0,1.0};
    double histScales[size] = {1,2,4,20,40,60};


    // initialize files
    
    //Form("root_output_files/output00%d.root",dataNumbers[4]),Form("root_output_files/output00%d.root",dataNumbers[5])
    vector<TH1F*> histograms = histCollector(overVoltages, files, temp, size);
    vector<TH1F*> histogramsUnshadowed = unShadower(histograms);
    // output file
    TFile *out = new TFile(Form("histOutput%d.root",temp), "RECREATE");

    for (int count = 0; count < size; count ++){

   
        histograms[count]->Scale(histScales[count]);
        histograms[count]->Scale(1.,"width");
        histograms[count]->Fit("expo","WL","",fitStartRange[count],fitEndRange[count]); // L specifies log likelihood (which deals with the non-gaussian bin statistics in low count bins). We only fit after the first several bins to ignore afterpulsing. 

        fit.push_back(histograms[count]->GetFunction("expo")); // save the fit paramters to a vector
        fitUS.push_back(histogramsUnshadowed[count]->GetFunction("pol0")); // save the fit paramters to a vector

        //write histo and fit
        histograms[count]->Write(Form("Hist%d",count));
        fit[count]->Write(Form("Fit%d",count));
//  scaled
        fitRate[count] = (fit[count]->GetParameter(1))*(-1000000000)/(1.296);
        fitError[count] = fit[count]->GetParError(1)*(1000000000)/(1.296);
        fitRateUS[count] = fitUS[count]->GetParameter(0)*(1000000000)/(1.296);
        fitErrorUS[count] = fitUS[count]->GetParError(0)*(1000000000)/(1.296);

    }

// TGraph 
    auto c1 = new TCanvas("c1","HoV vs Slope Fit",200,10,700,500);
    c1->SetFillColor(0);
    c1->SetGridx();
    c1->SetGridy();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);

    auto gr = new TGraphErrors(size, pulseWidthNS,fitRate,overVolErrors,fitError);
    gr->SetMarkerStyle(22);
    gr->GetXaxis()->SetTitle("Pulse width [ns]");
    //scaled
    gr->GetYaxis()->SetTitle("DCR [Hz/]");
    //unscaled
    //gr->GetYaxis()->SetTitle("DCR [Hz]");
    gr->Draw();
    gr->Write(Form("Graph%d",temp));

//Time Difference Plot
    auto c2 = new TCanvas("c2","Time Difference Distributions",200,10,700,500);
    c2->SetFillColor(0);
    c2->SetGridx();
    c2->SetGridy();
    c2->GetFrame()->SetFillColor(21);
    c2->GetFrame()->SetBorderSize(12);
    
    gPad->SetLogy();
    gPad->SetLogx();

    histograms[0]->Draw("PLC PMC");

    histograms[size-1]->Draw("PLC PMC");
    for (int num = 0; num < size-1; num++){
        histograms[size-2-num]->Draw("same PLC PMC");
    }
    gPad->BuildLegend();
// Afterpulsing plot
    auto c3 = new TCanvas("c3","Unshadowed Distributions",200,10,700,500);
    c3->SetFillColor(0);
    c3->SetGridx();
    c3->SetGridy();
    c3->GetFrame()->SetFillColor(21);
    c3->GetFrame()->SetBorderSize(12);

    gPad->SetLogx();  
    gPad->SetLogy();
    std::cout << "size of unshadowed: " <<histogramsUnshadowed.size()<<std::endl;

    histogramsUnshadowed[size-1]->Draw(" PLC PMC");
    for (int num = 1; num < size; num++){
        histogramsUnshadowed[size-num-1]->Draw(" same PLC PMC");
    }
    gPad->BuildLegend();
    
    ofstream myFile;
    myFile.open(Form("%d_VoV_DNR.txt",temp));
    myFile << "VoV" <<", ";
    myFile << "VoV Error" <<", ";
    myFile << "DNR" <<", ";
    myFile << "DNR Error" <<" "<<std::endl;
    for (int n = 0; n<size; n++){
        myFile << overVoltages[n] <<", ";
        myFile << overVolErrors[n] <<", ";
        myFile << fitRate[n] <<", ";
        myFile << fitError[n] <<" "<<std::endl;

    }
    myFile.close();

}