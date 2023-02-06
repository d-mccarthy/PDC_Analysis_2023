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

    const char* files[7] = {"histOutput40.root","histOutput60.root","histOutput80.root","histOutput100.root","histOutput120.root","histOutput140.root","histOutput160.root"};
    const char* graphs[7] = {"Graph40","Graph60","Graph80","Graph100","Graph120","Graph140","Graph160"};

    for(int i = 0; i < 7; i++)
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
        errorGraphs[i]->SetMarkerStyle(104+i);
        errorGraphs[i]->SetTitle(Form("Temp: %d C",-40-(20*i)));
        delete myFile;
    }

    auto c1 = new TCanvas("c1","VoV vs Slope Fit",200,10,700,500);
    c1->SetFillColor(0);
    c1->SetGridx();
    c1->SetGridy();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);
    
    // 

    errorGraphs[0]->SetMarkerColor(1);
    errorGraphs[1]->SetMarkerColor(2);
    errorGraphs[2]->SetMarkerColor(3);
    errorGraphs[3]->SetMarkerColor(4);
    errorGraphs[4]->SetMarkerColor(5);
    errorGraphs[5]->SetMarkerColor(6);
    errorGraphs[6]->SetMarkerColor(7);

    stackOfTemps->Add(errorGraphs[0]);
    stackOfTemps->Add(errorGraphs[1]);
    stackOfTemps->Add(errorGraphs[2]);
    stackOfTemps->Add(errorGraphs[3]);
    stackOfTemps->Add(errorGraphs[4]);
    stackOfTemps->Add(errorGraphs[5]);
    stackOfTemps->Add(errorGraphs[6]);

    stackOfTemps->GetXaxis()->SetTitle("VoV [V]");
    stackOfTemps->GetYaxis()->SetTitle("DCR [Hz/m^{2}]");

    stackOfTemps->Draw("ALP");
    c1->BuildLegend();


}