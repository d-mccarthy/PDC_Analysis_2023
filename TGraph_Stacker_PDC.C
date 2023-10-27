#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <typeinfo>
#include "TSystem.h"
#include "TClass.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TFile.h"
#include "TMultiGraph.h"

using namespace std;

void TGraph_Stacker_PDC(){

    TMultiGraph *stackOfTemps = new TMultiGraph();

    vector<TGraphErrors*> errorGraphs;

    const int size = 7;
    //naming is bad right now, silly mistake -- 60 is -60, 80 is -80
    const char* files[size] = {"histOutput90.root","histOutput_133_7.root","histOutput_163_7.root","histOutput200.root","histOutput_223_32.root","histOutput_273_7.root","histOutput300.root"};
    const char* graphs[size] = {"Graph90","Graph133_anode7","Graph163_anode7","Graph200","Graph223_anode32","Graph273_anode7","Graph300"};
    const char* temps[size] = {"93","133","163","193","223","273","293"};
    const float temperatures[size] = {93,133,163,193,223,273,293};

    auto tempGr = new TGraphErrors(7);
    tempGr->SetTitle(Form("DCR vs Temp @4VoV"));
    tempGr->GetXaxis()->SetTitle("Temp [K]");
    tempGr->GetYaxis()->SetTitle("DCR [Hz/m^{2}]");

    for(int i = 0; i < size; i++)
    {
        TFile *myFile = TFile::Open(files[i]);
        if(gSystem->AccessPathName(files[i]))
        {
            std::cout << "file does not exist" << std::endl;
        } 
        else 
        {
            std::cout << "file exists" << std::endl;
        };

        errorGraphs.push_back( (TGraphErrors*)myFile->Get(graphs[i]) );
        errorGraphs[i]->SetMarkerStyle(104);
        errorGraphs[i]->SetTitle(Form("Temp: %s C",temps[i]));
        
        if (temperatures[i]<250 || temperatures[i]==293){
            std::cout << errorGraphs[i]->GetPointX(0)<<std::endl;
            std::cout << errorGraphs[i]->GetPointY(0)<< " " << errorGraphs[i]->GetErrorY(0) << std::endl;
            std::cout << temperatures[i]<< " " << 5 << std::endl;
            tempGr->SetPoint(i,temperatures[i],errorGraphs[i]->GetPointY(0));
            tempGr->SetPointError(i,5,errorGraphs[i]->GetErrorY(0));
        }
        else{
            std::cout << errorGraphs[i]->GetPointX(3)<<std::endl;
            std::cout << errorGraphs[i]->GetPointY(3)<< " " << errorGraphs[i]->GetErrorY(3) << std::endl;
            std::cout << temperatures[i]<< " " << 5 << std::endl;
            tempGr->SetPoint(i,temperatures[i],errorGraphs[i]->GetPointY(3));
            tempGr->SetPointError(i,5,errorGraphs[i]->GetErrorY(3));
        }
        

        delete myFile;
    }

    auto c1 = new TCanvas("c1","VoV vs Slope Fit",200,10,700,500);
    c1->SetFillColor(0);
    c1->SetGridx();
    c1->SetGridy();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);
    
    // 
    for (int j = 0; j<size;j++){
        errorGraphs[j]->SetMarkerColor(j+1);
        stackOfTemps->Add(errorGraphs[j]);
    
    }

    gPad->SetLogy();

    stackOfTemps->SetTitle(Form("DCR vs VoV for 7 Temperatures"));
    stackOfTemps->GetXaxis()->SetTitle("VoV [V]");
    stackOfTemps->GetYaxis()->SetTitle("DCR [Hz/m^{2}]");

    stackOfTemps->Draw("ALP");
    c1->BuildLegend();

    auto c2 = new TCanvas("c2","Temp vs Slope Fit",200,10,700,500);
    c2->SetFillColor(0);
    c2->SetGridx();
    c2->SetGridy();
    c2->GetFrame()->SetFillColor(21);
    c2->GetFrame()->SetBorderSize(12);
    
    gPad->SetLogy();
    tempGr->SetMarkerStyle(24);
    tempGr->Draw("ALP");


}