#include <iostream>
#include <fstream>
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
// Global Variables
const int fileNumber = 5; // will need to change number of filenames in files array to match (should do as a vector, but root was seg faulting for me)
double pulseWidth[5] = {83.0,83.0,83.0,83.0,83.0}; //in ADC counts -- 50 ns pulses are 25 ADC
double binArray[91]; //defines total number of bins


//define a class to determine the temperature dependant data (i.e. fit end ranges, VOVs measured)
class TemperatureDependentData {
public:
    int temperature;
    std::vector<double> overVoltage;
    std::vector<float> fitRangeStart;
    std::vector<float> fitRangeEnd;

    TemperatureDependentData(int temp, const std::vector<double>& OV, const std::vector<float>& fRS, const std::vector<float>& fES){
        temperature = temp;
        overVoltage = OV;
        fitRangeStart = fRS;
        fitRangeEnd = fES;

    };
};

// Function to give the afterpulsing distribution with the probability distribution function of the total DNR as input.
// vector<TH1F*> unShadower(vector<TH1F*> inputHist, vector<double> VoV, vector<float> fitStart, vector<float> fitEnd){

//     //define the binning of the histogram. This was done by eye to give reasonable statistics in each bin.
//     for(int j =0; j<91; j++){

//         binArray[j]= ((pow(2.7162,j*0.25)*1.1)/1E7);    
        
//     }
//     // initialize a vector of histograms to give as output
//     std::vector< TH1F* > unShadowed;
//     // initialize variables: following from 
//     double lambda;
//     double beta;
//     double sigPi;
//     double sigLambda;
//     double sigBeta;
//     double Pi;
//     double PiWidth;
//     int numBins;

//     for(int hInx = 0; hInx<inputHist.size(); hInx++){
//         unShadowed.push_back(new TH1F(Form("unShadowedHist%d",hInx),Form("VoV = %f [V]",VoV[hInx]),90,binArray));
//         beta=0.;
//         sigBeta=0.;
//         numBins = inputHist[hInx]->GetNbinsX();


//         for(int ti = 0; ti < numBins-5 ; ti++){

//             Pi = inputHist[hInx]->GetBinContent(ti);

//             sigPi = (inputHist[hInx]->GetBinError(ti));
           
//             lambda = -log(1-(Pi/exp(-beta)));

//             sigLambda = sqrt(sigPi*sigPi+Pi*Pi*sigBeta*sigBeta)/(exp(-beta)-Pi);
//             sigBeta = sqrt(sigPi*sigPi+exp(-2*beta)*sigBeta*sigBeta)/(exp(-beta)-Pi);
            
//             beta += lambda;

//             if(!isnan(lambda) && sigLambda < lambda){
//                 unShadowed[hInx]->SetBinContent(ti,lambda);
//                 unShadowed[hInx]->SetBinError(ti,sigLambda);
//             }
            
//         }
//         unShadowed[hInx]->Sumw2();
//         unShadowed[hInx]->Scale(1./unShadowed[hInx]->Integral());  
//         unShadowed[hInx]->Scale(1.,"width");
//         unShadowed[hInx]->Fit("pol0","WL","",fitStart[hInx],fitEnd[hInx]);

//         unShadowed[hInx]->GetYaxis()->SetTitle("Pulse Rate [Hz/m^{2}]");
//         unShadowed[hInx]->GetXaxis()->SetTitle("Time After Primary Pulse [s]");
//     }

//     return unShadowed;
        
// }
//for use in dealing with back to back pulses
// void filler(int repeats, TH1F* toFill, float width){
//     for(int inx = 0; inx<repeats;inx++){
//         toFill->Fill(width*(2.0/1E9));
//     }
    
// }

vector<TH1F*> histCollector(vector<double> VoV, const char *filenames[], int temp, int inSize)
{
    vector<TH1F*> hist;

    double windowResolution = 2.0/1E9; //time resolution within window is 2ns
    double triggerResolution = 16.0/1E9; //time resolution within full run is 16ns

    double binArray[91];

    const int length = inSize;
    

    //rebin histograms so all bins have statistics
    for(int j =0; j<91; j++){

        binArray[j]= ((pow(2.7162,j*0.25)*1.1)/1E7);
        
    }
    
    for(int i=0; i<inSize; i++)
    {
        //create histogram
        hist.push_back(new TH1F(Form("hist%d",i),Form("VoV [V] = %f, Temp[K] = %d",VoV[i],temp), 90, binArray));

        //open file and TTreeReader
        if(gSystem->AccessPathName(filenames[i]))
        {
            std::cout << "file does not exist" << std::endl;
            gApplication->Terminate();
        } 
        else 
        {
            std::cout << "file exists" << std::endl;
        };


        //open file
        TFile *myFile = TFile::Open(filenames[i]);
        TTreeReader myReader("TPulse", myFile);

        //open timing data from file
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

                
                timeNow = ((myTrigTIME[j])+trigCounter*2147483648)*(triggerResolution) + (myTIME[j]) * (windowResolution);
                //timeNow = ((myTrigTIME[j])+trigCounter*2147483648)*(triggerResolution);
                
                timeBetween = timeNow - timeBefore;
                
                //hist[i]->Fill(myTIME[j]);
                hist[i]->Fill(timeBetween);
                
                //This is for resolving back-to-back pulses in the case where there are no overlapping pulses (1 SPAD)
                /*********************
                if(myWIDTH[j]>(7*pulseWidth[i]+5))
                {
                    filler(7,hist[i],pulseWidth[i]);
                    timeNow += 7*pulseWidth[i]*(2.0/1E9);
                }
                else if(myWIDTH[j]>(6*pulseWidth[i]+5))
                {
                    filler(6,hist[i],pulseWidth[i]);
                    timeNow += 6*pulseWidth[i]*(2.0/1E9);
                }
                else if(myWIDTH[j]>(5*pulseWidth[i]+5))
                {
                    filler(5,hist[i],pulseWidth[i]);
                    timeNow += 5*pulseWidth[i]*(2.0/1E9);
                }
                else if(myWIDTH[j]>(4*pulseWidth[i]+5))
                {
                    filler(4,hist[i],pulseWidth[i]);
                    timeNow += 4*pulseWidth[i]*(2.0/1E9);
                }
                else if(myWIDTH[j]>(3*pulseWidth[i]+5))
                {
                    filler(3,hist[i],pulseWidth[i]);
                    timeNow += 3*pulseWidth[i]*(2.0/1E9);

                }
                else if(myWIDTH[j]>(2*pulseWidth[i]+5))
                {
                    filler(2,hist[i],pulseWidth[i]);
                    timeNow += 2*pulseWidth[i]*(2.0/1E9);

                }
                else if(myWIDTH[j]>(pulseWidth[i]+10))
                {
                    filler(1,hist[i],pulseWidth[i]);
                    timeNow += pulseWidth[i]*(2.0/1E9);

                }
                **********************************************/
                timeBefore = timeNow;
                trigLast = myTrigTIME[j];
                
            }
            
        }   

        //get error per each bin and then scale each bin by its width to normalize
        hist[i]->Sumw2();
        //hist[i]->Scale(1.,"width");
        hist[i]->Scale(1./hist[i]->Integral());  

        hist[i]->GetYaxis()->SetTitle("Pulse Rate [Hz/m^{2}]");
        hist[i]->GetXaxis()->SetTitle("Time After Primary Pulse [s]");

        hist[i]->SetMarkerColor(10*(i+2));
        delete myFile;
    }
    return hist;
}
//input should be (csv run numbers, anode)
void PDC_DarkAnalysis(
    int inRun1 = 2298,
    int inRun2 = 2299,
    int inRun3 = 2300,
    int inRun4 = 2301,
    int inRun5 = 2302,
    int inSPAD = 7, 
    int inTemp = 293)
    
    // 440nm IR-VIS PDE runs
    //'R02298','R02299','R02300','R02301','R02302'
{

    //store the fits in a vector to retrieve later for error propagation
    std::vector< TF1* > fit;
    std::vector< TF1* > fitUS;
    std::vector<TemperatureDependentData> tempData;

    
    //create a look up to find the correct fit range for each temperature
    //fit ranges given in seconds, temperatures in K. Temp, VoV,FitStart,FitEnd
    tempData.push_back(TemperatureDependentData(293,{1,2,3,4,5},{5e-7,5e-7,1E-4,1E-4,1E-4},{2E-2,1.5E-2,1E-2,9E-3,8E-3}));
    tempData.push_back(TemperatureDependentData(273,{1,2,3,4,5},{1E-3,1E-3,9E-4,8E-4,8E-4},{9E-2,7E-2,5E-2,3E-2,2E-2}));
    tempData.push_back(TemperatureDependentData(223,{4,5,6,7,8},{3E-3,4E-3,2E-3,1.3E-3,1E-3},{2E-1,1.8E-1,9E-2,7E-2,4E-2}));
    tempData.push_back(TemperatureDependentData(193,{4,5,6,7,8},{3E-3,2.8E-3,2.4E-3,2.2E-3,2E-3},{1.4E-1,1E-1,9E-2,7E-2,5E-2}));
    tempData.push_back(TemperatureDependentData(163,{4,5,6,7,8},{3E-3,2.8E-3,2.4E-3,2.2E-3,2E-3},{1.4E-1,1E-1,9E-2,7E-2,5E-2}));
    tempData.push_back(TemperatureDependentData(133,{4,5,6,7,8},{3E-3,2.8E-3,2.4E-3,2.2E-3,2E-3},{1.4E-1,1E-1,9E-2,7E-2,5E-2}));
    tempData.push_back(TemperatureDependentData(93,{4,5,6,7,8},{1E-3,1E-3,1E-3,1E-3,1E-3},{3,2.5,2,1.5,1}));
    tempData.push_back(TemperatureDependentData(276,{1,2,3,4,5},{1E-3,1E-3,9E-4,8E-4,8E-4},{9E-2,7E-2,5E-2,3E-2,2E-2}));
    tempData.push_back(TemperatureDependentData(233,{4,5,6,7,8},{3E-3,4E-3,2E-3,1.3E-3,1E-3},{2E-1,1.8E-1,9E-2,7E-2,4E-2}));
    tempData.push_back(TemperatureDependentData(208,{4,5,6,7,8},{3E-3,2.8E-3,2.4E-3,2.2E-3,2E-3},{1.6,0.7,0.5,0.3,0.1}));
    tempData.push_back(TemperatureDependentData(163,{4,5,6,7,8},{3E-3,2.8E-3,2.4E-3,2.2E-3,2E-3},{1.4E-1,1E-1,9E-2,7E-2,5E-2}));
    tempData.push_back(TemperatureDependentData(133,{4,5,6,7,8},{3E-3,2.8E-3,2.4E-3,2.2E-3,2E-3},{1.4E-1,1E-1,9E-2,7E-2,5E-2}));
    tempData.push_back(TemperatureDependentData(117,{4,5,6,7,8},{1E-3,1E-3,1E-3,1E-3,1E-3},{3,2.5,2,1.5,1}));
    
    int anode = inSPAD;
    int temp = inTemp;
    double area = (1E9)/(1.296);

    vector<float> fitStartRange;
    vector<float> fitEndRange;
    vector<double> overV;

    bool found = false;

    for (const TemperatureDependentData& item :tempData){
        if (temp == item.temperature){
            found = true;
            overV = item.overVoltage;
            fitStartRange = item.fitRangeStart;
            fitEndRange = item.fitRangeEnd;
            break;
        }
        else
        {
            continue;
        }
    }
    if (!found){
        cout<< "temperature not found, resorting to defaults" <<endl;
        overV = tempData[0].overVoltage;
        fitStartRange = tempData[0].fitRangeStart;
        fitEndRange = tempData[0].fitRangeEnd;
    }
    

    double fitRate[fileNumber];
    double fitError[fileNumber];
    double fitRateUS[fileNumber];
    double fitErrorUS[fileNumber];
    double fitChiSqr[fileNumber];
    double fitChiSqrUS[fileNumber];
    double reducedChiSquaredErrors[fileNumber];
    double DoF = 88.0;

    int dataNumbers[fileNumber] = {inRun1,inRun2,inRun3,inRun4,inRun5};

    // initialize files
    const char *files[fileNumber] = 
    {
    Form("root_output_files/output0%d.root",dataNumbers[0]),
    Form("root_output_files/output0%d.root",dataNumbers[1]),
    Form("root_output_files/output0%d.root",dataNumbers[2]),
    Form("root_output_files/output0%d.root",dataNumbers[3]),
    Form("root_output_files/output0%d.root",dataNumbers[4])
    };

    double overVolErrors[fileNumber] = {0.1,0.1,0.1,0.1,0.1};
   
    
    
    vector<TH1F*> histograms = histCollector(overV, files, temp, fileNumber);
    //vector<TH1F*> histogramsUnshadowed = unShadower(histograms, overV,fitStartRange,fitEndRange);
    //vector<TH1F*> histogramsUnshadowedSubtracted;

    // output file

    TFile *out = new TFile(Form("PDCOutput/histOutput_%d.root",temp), "UPDATE");

    for (int count = 0; count < fileNumber; count ++){

        //Scale to plot expo fit nicely
        histograms[count]->Scale(1.,"width");
        //histograms[count]->Fit("expo","WL","",fitStartRange[count],fitEndRange[count]*7); // L specifies log likelihood (which deals with the non-gaussian bin statistics in low count bins). We only fit on long time scales to ignore afterpulsing. 
        
        //USING FIXED START END END POINTS FOR FIT
        histograms[count]->Fit("expo","WL","",5e-6,1e-5); // L specifies log likelihood (which deals with the non-gaussian bin statistics in low count bins). We only fit on long time scales to ignore afterpulsing. 

        fit.push_back(histograms[count]->GetFunction("expo")); // save the fit paramters to a vector
        //fitUS.push_back(histogramsUnshadowed[count]->GetFunction("pol0"));

        //write histo and fit
        histograms[count]->Write(Form("Hist%d_anode%d",count,anode));
        fit[count]->Write(Form("Fit%d_anode%d",count,anode));

        fitRate[count] = (fit[count]->GetParameter(1))*-area;
        fitError[count] = fit[count]->GetParError(1)*area;
        fitRateUS[count] = fitUS[count]->GetParameter(0);
        fitErrorUS[count] = fitUS[count]->GetParError(0);
        fitChiSqr[count] = fit[count]->GetChisquare();
        fitChiSqrUS[count] = fitUS[count]->GetChisquare();

    }
    // scale the errors by the reduced chi squared
    for (size_t i = 0; i < fileNumber; i++) {
            reducedChiSquaredErrors[i] = fitError[i] * fitChiSqr[i]/88;
        }

    // make the graphs

    auto c1 = new TCanvas("c1","VoV vs Slope Fit",200,10,700,500);
    c1->SetFillColor(0);
    c1->SetGridx();
    c1->SetGridy();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);

    auto gr = new TGraphErrors(fileNumber, overV.data(),fitRate,overVolErrors,reducedChiSquaredErrors);
    gr->SetMarkerStyle(22);
    gr->GetXaxis()->SetTitle("VoV [V]");
    gr->GetYaxis()->SetTitle("DCR [Hz/m^{2}]");
    gr->Draw();
    gr->Write(Form("Graph%d_anode%d",temp,anode));


    auto c2 = new TCanvas("c2","Time Difference Distributions",200,10,700,1100);
    c2->SetFillColor(0);
    c2->Divide(2,3);
    c2->GetFrame()->SetFillColor(21);
    c2->GetFrame()->SetBorderSize(12); 

    THStack *tempStack = new THStack("tempStack",Form("Time Difference Distribution @ %d K", temp));
    for (int i = 0; i < 5; i++)
    {
        histograms[i]->SetLineColor(10*(i+2));
        histograms[i]->SetStats(0);
        histograms[i]->GetXaxis()->SetTitleOffset(1.2);
        
    }
    
    gStyle->SetOptFit(1);
    tempStack->Add(histograms[0]);
    tempStack->Add(histograms[1]);
    tempStack->Add(histograms[2]);
    tempStack->Add(histograms[3]);
    tempStack->Add(histograms[4]);

    tempStack->Draw("PADS hist");
    gPad->Update();

    for (int i = 1; i <= 5; i++)
    {
    c2->cd(i);
    gPad->SetLogy(i);
    gPad->SetLogx(i);  
    }  
    
    TLegend* legend = new TLegend(0.49, 0.8, 0.9, 0.9); // Adjust the position and size as needed
    legend->AddEntry(histograms[0], "VoV = 4V", "f");
    legend->AddEntry(histograms[1], "VoV = 5V", "f");
    legend->AddEntry(histograms[2], "VoV = 6V", "f");
    legend->AddEntry(histograms[3], "VoV = 7V", "f");
    legend->AddEntry(histograms[4], "VoV = 8V", "f");
    legend->SetTextSize(0.02);
    //legend->Draw();


   //Afterpulsing plot
    // auto c3 = new TCanvas("c3","Unshadowed Distributions",200,10,700,500);
    // c3->SetFillColor(0);
    // c3->SetGridx();
    // c3->SetGridy();
    // c3->GetFrame()->SetFillColor(21);
    // c3->GetFrame()->SetBorderSize(12);

    // gPad->SetLogy();
    // gPad->SetLogx();  


    float currentVal;
    float currentErr;

    // for (int i = 0; i < fileNumber; i++)
    // {
    //     histogramsUnshadowedSubtracted.push_back(new TH1F(Form("unShadowedHistSubtracted%d",i),Form("VoV = %f [V]",overV[i]),90,binArray));
    //     for (int ti = 0; ti < histogramsUnshadowed[0]->GetNbinsX();ti++)
    //     {
    //         currentVal = histogramsUnshadowed[i]->GetBinContent(ti);
    //         currentErr = histogramsUnshadowed[i]->GetBinError(ti);
    //         histogramsUnshadowedSubtracted[i]->SetBinContent(ti, currentVal - fitRateUS[i]);
    //         histogramsUnshadowedSubtracted[i]->SetBinError(ti,pow(pow(currentErr,2)+pow(fitErrorUS[i],2),.5));
    //     }
    // }
    

    // histogramsUnshadowed[0]->Draw("PLC PMC");
    // histogramsUnshadowed[1]->Draw("same PLC PMC");
    // histogramsUnshadowed[2]->Draw("same PLC PMC");
    // histogramsUnshadowed[3]->Draw("same PLC PMC");
    // histogramsUnshadowed[4]->Draw("same PLC PMC");

    // histogramsUnshadowedSubtracted[0]->Draw("hist PLC PMC");
    // histogramsUnshadowedSubtracted[1]->Draw("same hist PLC PMC");
    // histogramsUnshadowedSubtracted[2]->Draw("same hist PLC PMC");
    // histogramsUnshadowedSubtracted[3]->Draw("same hist PLC PMC");
    // histogramsUnshadowedSubtracted[4]->Draw("same hist PLC PMC");

    gPad->BuildLegend();

    
    
    double integrals[fileNumber];
    double intError[fileNumber];
    double VOver[fileNumber] = {1,2,3,4,5};
    double VOverError[fileNumber] = {0.1,0.1,0.1,0.1,0.1};

    for (int i = 0; i < fileNumber; i++)
    {   
        //integrals[i] = histogramsUnshadowedSubtracted[i]->IntegralAndError(0,histogramsUnshadowedSubtracted[i]->FindBin(1E-5),intError[i],"width");
        intError[i] = intError[i];
    }
    cout<<"afterpulsing and error at 4 VoV: "<<endl;
    cout<<integrals[3]<<endl;
    cout<<intError[3]<<endl;

    auto c4 = new TCanvas("c4","VoV vs Afterpulsing",200,10,700,500);
    c4->SetFillColor(0);
    c4->SetGridx();
    c4->SetGridy();
    c4->GetFrame()->SetFillColor(21);
    c4->GetFrame()->SetBorderSize(12);
    c4->SetLeftMargin(0.12);

    auto gr2 = new TGraphErrors(fileNumber, VOver,integrals,VOverError,intError);
    gr2->SetMarkerStyle(22);
    gr2->GetXaxis()->SetTitle("VoV [V]");
    gr2->GetYaxis()->SetTitle("Number of Afterpulses/Pulse");
    gr2->Draw();
   
    cout << "Fit DCR at 4VoV = " << fitRate[0] << "[Hz/m^{2}]" << endl;
    cout << "Fit DCR Error at 4VoV = " << fitError[0] << endl;


}