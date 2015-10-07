#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <algorithm>
//#include <TTreeReader.h>
//#include <TTreeReaderValue.h>

void loadRootFiles(const std::string fileName[], int nSamples, TTree *tree[] , int Length[]){
    TFile *f;
    TDirectory *dir;
    std::string dirName;

    for(int i = 0; i < nSamples; i++){
        dirName = fileName[i] + ".root";
        f = TFile::Open(dirName.c_str());
        dirName = dirName + ":/NtupleAnalyzer";
        dir = (TDirectory*)f->Get(dirName.c_str());
        dir->ls();
        tree[i] =  (TTree*)dir->Get("ntuple");
        Length[i] = tree[i]->GetEntries();
        std::cout << Length[i] << std::endl;
    }
}

void makeHistograms(const std::string sampleNames[], const int nSamples,  const std::string doubleVariableNames[],const int nDoubleVariables, const int nBins, const int minValue[], const int maxValue[], TH1F **TH1F_histograms){

    for (int i = 0; i < nSamples; i++){
        TH1F_histograms[i] = new TH1F [nDoubleVariables];
    }

    std::string histName;
    for(int i = 0; i < nSamples; i++){
        for(int j = 0; j < nDoubleVariables; j++){
            histName = "TH1F_" + sampleNames[i] + "_" + doubleVariableNames[j];
            std::cout << histName << std::endl;
            TH1F_histograms[i][j].SetName(histName.c_str());
            TH1F_histograms[i][j].SetTitle(histName.c_str());
            TH1F_histograms[i][j].SetBinsLength(nBins);
            TH1F_histograms[i][j].SetMinimum(minValue[j]);
            TH1F_histograms[i][j].SetMaximum(maxValue[j]);
        }
    } 

    
}

void getBranches(const std::string doubleVariableNames[], const int nDoubleVariables, const int doubleArraySize[], double  ***doubleVariables, TTree *sampleTrees[], const int nSamples){

    for(int i = 0; i < nSamples; i++){
        doubleVariables[i] = new double*[nDoubleVariables];
        for(int j = 0; j < nDoubleVariables; j++){
            doubleVariables[i][j] = new double [doubleArraySize[j]];
        }
        
    }



}


void plotHistograms(){

}

void drawHistograms(){

}

void correlationStudy(){
    
    ///constants
    const std::string sampleNames[] = {"MonoTop", "WJet", "TTbar"};
    const int nSamples = sizeof(sampleNames)/sizeof(sampleNames[0]);

    const std::string doubleVariableNames[] =   {"MetPt",   "PFMuonPt"};
    const int minValue[] =                      {10,        10};
    const int maxValue[] =                      {700,       700};
    const int doubleArraySize[] =               {40,        13};

    const int nBins = 100;
    const int nDoubleVariables = sizeof(doubleVariableNames)/sizeof(doubleVariableNames[0]);
    

    ///variables
    int sampleLength[nSamples];
    double ***doubleVariables = new double**[nSamples];

    //root objects
    TTree *sampleTrees[nSamples];

    TH1F **TH1F_histograms = new TH1F*[nSamples];


    loadRootFiles(sampleNames, nSamples, sampleTrees, sampleLength);

    makeHistograms(sampleNames, nSamples, doubleVariableNames, nDoubleVariables, nBins, minValue, maxValue,  TH1F_histograms);

    getBranches(doubleVariableNames, nDoubleVariables, doubleArraySize, doubleVariables, sampleTrees, nSamples);    

    plotHistograms();

    drawHistograms();

}
