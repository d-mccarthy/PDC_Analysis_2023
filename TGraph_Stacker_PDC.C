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

    const int size = 9;
    //naming is bad right now, silly mistake -- 60 is -60, 80 is -80
    const char* files[size] = {"histOutput170.root","histOutput190.root","histOutput210.root","histOutput230.root","histOutput240.root","histOutput250.root","histOutput260.root","histOutput270.root","histOutput300.root"};
    const char* graphs[size] = {"Graph170","Graph190","Graph210","Graph230","Graph240","Graph250","Graph260","Graph270","Graph300"};
    const char* temps[size] = {"-100","-80","-60","-40","-30","-20","-10","0","20"};
    const float temperatures[size] = {175.0,195.0,215.0,235.0,245.0,255.0,265.0,275.0,300.0};

    auto tempGr = new TGraphErrors(5);
    tempGr->SetTitle(Form("DNR vs Temp @4VoV"));
    tempGr->GetXaxis()->SetTitle("Temp [K]");
    tempGr->GetYaxis()->SetTitle("DCR [Hz/m^{2}]");

    for(int i = 0; i < size; i++)
    {
        TFile *myFile = TFile::Open(files[i]);
        // if(gSystem->AccessPathName(files[i]))
        // {
        //     std::cout << "file does not exist" << std::endl;
        // } 
        // else 
        // {
        //     std::cout << "file exists" << std::endl;
        // };

        errorGraphs.push_back( (TGraphErrors*)myFile->Get(graphs[i]) );
        errorGraphs[i]->SetMarkerStyle(104);
        errorGraphs[i]->SetTitle(Form("Temp: %s C",temps[i]));
        
        if (i<3){
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

    tempGr->Draw("ALP");


}