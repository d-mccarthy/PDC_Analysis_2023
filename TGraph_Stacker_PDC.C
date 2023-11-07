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
    TMultiGraph *stackForColoring = new TMultiGraph();
    std::vector< TF1* > fits;
    vector<TGraphErrors*> errorGraphs;

    const int size = 6;

    const char* temps[size] = {"117","176","208","233","276","293"};
    const int temperatures[size] = {117,176,208,233,276,293};
    const int anode[3] = {7,20,32};
    
    std::vector<const char* > files; 
    std::vector<std::vector<const char *> > graphs(6, std::vector<const char *>(3));
    
    

    for (int i = 0; i < size; i++){
        
        files.push_back(Form("PDCOutput/histOutput_%d.root",temperatures[i]));
        
        for (int j = 0; j < 3; j++){
            graphs[i][j] = (Form("Graph%d_anode%d",temperatures[i],anode[j]));
        }
    }

    auto tempGr = new TGraphErrors(size*3);
    auto tempGr2 = new TGraphErrors(size*3);
    tempGr->SetTitle(Form("DCR vs Temp @4VoV"));
    tempGr->GetXaxis()->SetTitle("Temp [K]");
    tempGr->GetYaxis()->SetTitle("DCR [Hz/m^{2}]");
    tempGr2->SetTitle(Form("DCR vs Temp @4VoV"));
    tempGr2->GetXaxis()->SetTitle("Temp [K]");
    tempGr2->GetYaxis()->SetTitle("DCR [Hz/m^{2}]");

    for(int i = 0; i < size; i++){
        TFile *myFile = TFile::Open(files[i]);
            if(gSystem->AccessPathName(files[i]))
            {
                std::cout << "file does not exist" << std::endl;
            } 
            else 
            {
                std::cout << "file exists" << std::endl;
            };
            
        for (int j = 0; j < 3; j++){
            errorGraphs.push_back( (TGraphErrors*)myFile->Get(graphs[i][j]) );
            
            errorGraphs[3*i+j]->SetMarkerStyle(104);
            cout << 1 <<endl;
            if (temperatures[i]<235){
                std::cout << errorGraphs[3*i+j]->GetPointX(0)<<std::endl;
                std::cout << errorGraphs[3*i+j]->GetPointY(0)<< " " << errorGraphs[3*i+j]->GetErrorY(0) << std::endl;
                std::cout << temperatures[i]<< " " << 5 << std::endl;

                if (errorGraphs[3*i+j]->GetErrorY(0) < errorGraphs[3*i+j]->GetPointY(0)*0.5){
                    tempGr2->SetPoint(3*i+j,temperatures[i],errorGraphs[3*i+j]->GetPointY(0));
                    tempGr2->SetPointError(3*i+j,5,errorGraphs[3*i+j]->GetErrorY(0));
                    tempGr2->SetMarkerColor(kRed);
                }
                else{
                    tempGr->SetPoint(3*i+j,temperatures[i],errorGraphs[3*i+j]->GetPointY(0));
                    tempGr->SetPointError(3*i+j,5,errorGraphs[3*i+j]->GetErrorY(0));
                }
                
                
            }
            else{
                std::cout << errorGraphs[3*i+j]->GetPointX(3)<<std::endl;
                std::cout << errorGraphs[3*i+j]->GetPointY(3)<< " " << errorGraphs[3*i+j]->GetErrorY(3) << std::endl;
                std::cout << temperatures[i]<< " " << 5 << std::endl;

                if (errorGraphs[3*i+j]->GetErrorY(3) < errorGraphs[3*i+j]->GetPointY(3)*0.5){
                    tempGr2->SetPoint(3*i+j,temperatures[i],errorGraphs[3*i+j]->GetPointY(3));
                    tempGr2->SetPointError(3*i+j,5,errorGraphs[3*i+j]->GetErrorY(3));
                    tempGr2->SetMarkerColor(kRed);
                }
                else{
                    tempGr->SetPoint(3*i+j,temperatures[i],errorGraphs[3*i+j]->GetPointY(3));
                    tempGr->SetPointError(3*i+j,5,errorGraphs[3*i+j]->GetErrorY(3));
                }
                
            }
            cout << 2 <<endl;
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
    for (int j = 0; j<size*3;j++){
        
        if (errorGraphs[j]->GetErrorY(0) < errorGraphs[j]->GetPointY(0)*0.5 && errorGraphs[j]->GetErrorY(3) < errorGraphs[j]->GetPointY(3)*0.5){
            errorGraphs[j]->SetMarkerColor(2);
            errorGraphs[j]->SetLineColor(2);
            stackOfTemps->Add(errorGraphs[j]);
        }
        else{
            errorGraphs[j]->SetMarkerColor(1);
        }
        if (j<3){
            errorGraphs[j]->SetMarkerStyle(53);
        }
        else if (3<=j && j<6){
            errorGraphs[j]->SetMarkerStyle(54);
        }
        else if (6<=j && j<9){
            errorGraphs[j]->SetMarkerStyle(55);
        }
        else if (9<=j && j<12){
            errorGraphs[j]->SetMarkerStyle(56);
        }
        else if (9<=j && j<12){
            errorGraphs[j]->SetMarkerStyle(57);
        }
        else if (12<=j && j<15){
            errorGraphs[j]->SetMarkerStyle(58);
        }
        else{
            errorGraphs[j]->SetMarkerStyle(59);
        }
        
    
    }

    gPad->SetLogy();

    stackOfTemps->SetTitle(Form("DCR vs VoV for 6 Temperatures"));
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
    tempGr->SetLineStyle(0);
    tempGr2->SetMarkerStyle(24);
    tempGr2->SetLineStyle(0);
    tempGr2->SetLineColor(2);
    stackForColoring->Add(tempGr);
    stackForColoring->Add(tempGr2);

    stackForColoring->Draw("AP");


}