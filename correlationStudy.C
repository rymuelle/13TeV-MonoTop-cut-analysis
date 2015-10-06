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

        f = TFile::Open(fileName[i].c_str());
        dirName = fileName[i] + ":/NtupleAnalyzer";
        dir = (TDirectory*)f->Get(dirName.c_str());
        dir->ls();
        tree[i] =  (TTree*)dir->Get("ntuple");
        Length[i] = tree[i]->GetEntries();
        std::cout << Length[i] << std::endl;
    }
}

void makeHistograms(){

}

void plotHistograms(){

}

void drawHistograms(){

}

void correlationStudy(){
    
    ///constants
    const std::string sampleNames[] = {"MonoTop.root", "WJet.root", "TTbar.root"};
    const int nSamples = sizeof(sampleNames)/sizeof(sampleNames[0]);

    ///variables
    int sampleLength[nSamples];
   
    TTree *sampleTrees[nSamples];

    loadRootFiles(sampleNames, nSamples, sampleTrees, sampleLength);

    makeHistograms();

    plotHistograms();

    drawHistograms();

}
