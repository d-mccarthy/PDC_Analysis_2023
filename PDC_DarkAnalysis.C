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

double pulseWidth[5] = {84.0,84.0}; //in ADC counts -- 50 ns pulses are 25 ADC
double overV[2] = {1,2};
double binArray[121];

vector<TH1F*> unShadower(vector<TH1F*> inputHist){

    for(int j =0; j<121; j++){

        binArray[j]= ((pow(2.7162,j*0.19))/1E9);    
        
    }

    std::vector< TH1F* > unShadowed; 
    double lambda;
    double beta;
    double sigPi;
    double sigLambda;
    double sigBeta;
    double Pi;
    double PiWidth;
    int numBins;

    for(int hInx = 0; hInx<inputHist.size(); hInx++){
        unShadowed.push_back(new TH1F(Form("unShadowedHist%d",hInx),Form("VoV = %f [V]",overV[hInx]),120,binArray));
        beta=0.;
        sigBeta=0.;
        numBins = inputHist[hInx]->GetNbinsX();

        for(int ti = 0; ti < numBins-5 ; ti++){

            Pi = inputHist[hInx]->GetBinContent(ti);

            sigPi = (inputHist[hInx]->GetBinError(ti));
            lambda = -log(1-(Pi/exp(-beta)));

            std::cout<< Pi << "  " << lambda << "  " << beta <<std::endl;

            sigLambda = sqrt(sigPi*sigPi+Pi*Pi*sigBeta*sigBeta)/(exp(-beta)-Pi);
            sigBeta = sqrt(sigPi*sigPi+exp(-2*beta)*sigBeta*sigBeta)/(exp(-beta)-Pi);
            
            beta += lambda;

            if (Pi != 0){
                std::cout << "lambda, ti, width = " << lambda << ", " << ti << ", "<<  inputHist[hInx]->GetBinWidth(ti) <<std::endl;
            }
            unShadowed[hInx]->SetBinContent(ti,lambda);
            //unShadowed[hInx]->SetBinError(ti,sigLambda);
        }
        for (int vi = 0; vi < unShadowed[hInx]->GetNbinsX(); vi++) {
            if (unShadowed[hInx]->GetBinCenter(vi) <= 0) {
                cout << "bad binning!" <<endl;
            }
        }
        unShadowed[hInx]->Scale(1.,"width");
    }

    return unShadowed;
        
}

void filler(int repeats, TH1F* toFill, float width){
    for(int inx = 0; inx<repeats;inx++){
        toFill->Fill(width*(2.0/1E9));
    }
    
}

vector<TH1F*> histCollector(float VoV[], const char *filenames[], int temp, int size)
{
    std::vector< TH1F* > hist;

    double windowResolution = 2.0/1000000000.0; //time resolution within window is 2ns
    double triggerResolution = 16.0/1000000000.0; //time resolution within full run is 16ns
    double binArray[121];

    const int length = size;
    

    //rebin histograms so all bins have statistics
    for(int j =0; j<121; j++){

        binArray[j]= ((pow(2.7162,j*0.19))/1E9);    
        
    }

    
    
    for(int i=0; i<size; i++)
    {
        //create histogram
        hist.push_back(new TH1F(Form("hist%d",i),Form("VoV = %f",VoV[i]), 120, binArray));
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
            //for only looking at timing within event windows
            // if (myTIME.GetSize()<2){
            //     continue;
            // }

            for (int j = 0; j<myTIME.GetSize(); j++)
            {
                // if(j==0){
                //     timeBefore = (myTrigTIME[j])*(triggerResolution);
                //     continue;
                // }
                if(myTrigTIME[j] < trigLast)
                {
                    trigCounter++;
                }
                
                timeNow = ((myTrigTIME[j])+trigCounter*2147483648)*(triggerResolution) + (myTIME[j]) * (windowResolution);

                //timeNow = (myTrigTIME[j]+trigCounter*2147483648)*(triggerResolution);
                //timeNow = (myTIME[j]*windowResolution);
                timeBetween = timeNow - timeBefore;

                hist[i]->Fill(timeBetween);

                if(myWIDTH[j]>(4*pulseWidth[i]+5))
                {
                    filler(4,hist[i],pulseWidth[i]);
                    timeNow += 4*pulseWidth[i]*(2.0/1000000000.0);
                    std::cout << myWIDTH[j] << std::endl;
                }
                else if(myWIDTH[j]>(3*pulseWidth[i]+5))
                {
                    filler(3,hist[i],pulseWidth[i]);
                    timeNow += 3*pulseWidth[i]*(2.0/1000000000.0);
                    std::cout << myWIDTH[j] << std::endl;
                }
                else if(myWIDTH[j]>(2*pulseWidth[i]+5))
                {
                    filler(2,hist[i],pulseWidth[i]);
                    timeNow += 2*pulseWidth[i]*(2.0/1000000000.0);
                    std::cout << myWIDTH[j] << std::endl;
                }
                else if(myWIDTH[j]>(pulseWidth[i]+12))
                {
                    filler(1,hist[i],pulseWidth[i]);
                    timeNow += pulseWidth[i]*(2.0/1000000000.0);
                    std::cout << myWIDTH[j] << std::endl;
                }

                timeBefore = timeNow;
                trigLast = myTrigTIME[j];

                
            }
            
        }   
        //get error per each bin and then scale each bin by its width to normalize
        hist[i]->Sumw2();
        //hist[i]->Scale(1.,"width");
        hist[i]->Scale(1./hist[i]->Integral());  

        hist[i]->GetYaxis()->SetTitle("Probability/s");
        hist[i]->GetXaxis()->SetTitle("Time After Primary Pulse [s]");

        delete myFile;

    }
    return hist;
}

void PDC_DarkAnalysis(){

    //store the fits in a vector to retrieve later for error propagation
    std::vector< TF1* > fit;

    const int size = 2; // will need to change number of filenames in files array to match (should do as a vector, but root was seg faulting for me)

    int temp = 230;

    float fitRate[size];
    float fitError[size];

    float temperature[1] = {230.0};

    int dataNumbers[size] = {698,699};
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


    float overVoltages[size] = {1,2};
    //old data {0.1,0.2,0.7,1.2,1.7,2.7};
    //new data 
    //room temp{0.5,1.5,3,4.5,6,7,8};
    //hold off{0.1,0.2,0.5,0.75,1,1.25,1.5,1.75}

    float fitEndRange[size] = {4.9,4.9 };

    const char *files[size] = {Form("root_output_files/output00%d.root",dataNumbers[0]),
    Form("root_output_files/output00%d.root",dataNumbers[1])
    };
    float overVolErrors[size] = {0.1,0.1};
    // initialize files
    
    //Form("root_output_files/output00%d.root",dataNumbers[4]),Form("root_output_files/output00%d.root",dataNumbers[5])
    vector<TH1F*> histograms = histCollector(overVoltages, files, temp, size);
    vector<TH1F*> histogramsUnshadowed = unShadower(histograms);
    // output file
    TFile *out = new TFile(Form("histOutput%d.root",temp), "RECREATE");

    for (int count = 0; count < size; count ++){

        //Scale to plot expo fit nicely
        
        histograms[count]->Fit("expo","WL","",0.01,fitEndRange[count]); // L specifies log likelihood (which deals with the non-gaussian bin statistics in low count bins). We only fit after the first several bins to ignore afterpulsing. 

        fit.push_back(histograms[count]->GetFunction("expo")); // save the fit paramters to a vector


        //write histo and fit
        histograms[count]->Write(Form("Hist%d",count));
        //fit[count]->Write(Form("Fit%d",count));

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


    auto c2 = new TCanvas("c2","Time Difference Distributions",200,10,700,500);
    c1->SetFillColor(0);
    c1->SetGridx();
    c1->SetGridy();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);

    gPad->SetLogy();
    gPad->SetLogx();    

    gStyle->SetPalette(kLightTemperature);

    histograms[0]->Scale(1.,"width");
    histograms[1]->Scale(1.,"width");
    histograms[0]->Draw("hist PLC PMC");
    for (int num = 1; num < size; num++){
        histograms[num]->Draw("same hist PLC PMC");
    }
    gPad->BuildLegend();


   // Afterpulsing plot
    auto c3 = new TCanvas("c3","Unshadowed Distributions",200,10,700,500);
    c3->SetFillColor(0);
    c3->SetGridx();
    c3->SetGridy();
    c3->GetFrame()->SetFillColor(21);
    c3->GetFrame()->SetBorderSize(12);

    gPad->SetLogy();
    gPad->SetLogx();  
    

    histogramsUnshadowed[0]->Draw("hist PLC PMC");
    histogramsUnshadowed[1]->Draw("hist same PLC PMC");
    // for (int num = 0; num < size; num++){
    //     cout << num << endl;
    //     histogramsUnshadowed[size-2-num]->Draw("hist same PLC PMC");
    // }
    
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