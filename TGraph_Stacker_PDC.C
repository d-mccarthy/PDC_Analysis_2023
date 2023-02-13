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

    const int size = 2;

    const char* files[size] = {"histOutput20.root","histOutput60.root"};
    const char* graphs[size] = {"Graph20","Graph60"};

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
    for (int j = 0; j<size;j++){
        errorGraphs[j]->SetMarkerColor(j+1);

        stackOfTemps->Add(errorGraphs[j]);
    
    }

    stackOfTemps->GetXaxis()->SetTitle("VoV [V]");
    stackOfTemps->GetYaxis()->SetTitle("DCR [Hz/m^{2}]");

    stackOfTemps->Draw("ALP");
    c1->BuildLegend();


}