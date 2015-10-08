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

void makeTH1FHistograms(const std::string sampleNames[], const int nSamples,  const std::string doubleVariableNames[],const int nDoubleVariables, const int nBins, const int minValue[], const int maxValue[], TH1F **TH1F_histograms){

    for (int i = 0; i < nSamples; i++){
        TH1F_histograms[i] = new TH1F [nDoubleVariables];
    }

    std::string histName;
    for(int i = 0; i < nSamples; i++){
        for(int j = 0; j < nDoubleVariables; j++){
            histName = "TH1F_" + sampleNames[i] + "_" + doubleVariableNames[j];
            std::cout << histName << std::endl;
            TH1F_histograms[i][j].SetName(histName.c_str()); //change to pointers with time
            TH1F_histograms[i][j].SetTitle(histName.c_str());
            TH1F_histograms[i][j].SetBinsLength(nBins);
            TH1F_histograms[i][j].SetMinimum(minValue[j]);
            TH1F_histograms[i][j].SetMaximum(maxValue[j]);
            std::cout<< maxValue[j] << std::endl;   
        }
    } 

    
}

void getBranches(const std::string doubleVariableNames[], const int nDoubleVariables, const int doubleArraySize[], double  ***doubleVariables, TTree *sampleTrees[], const int nSamples){


    for(int i = 0; i < nSamples; i++){
        doubleVariables[i] = new double*[nDoubleVariables];
        for(int j = 0; j < nDoubleVariables; j++){
            //doubleVariables[i][j] = new double [doubleArraySize[j]]; should work but doesn't? Look into with time. 
            doubleVariables[i][j] = new double [100];
            sampleTrees[i]->SetBranchAddress(doubleVariableNames[j].c_str(), &*doubleVariables[i][j]);
        }
        
    }


}


void plotTH1FHistograms(TTree *sampleTrees[], const int nSamples, double ***doubleVariables, const int nDoubleVariables, TH1F **TH1F_histograms, int  sampleLength[], const int maxEvents, const int doubleArrayIndex[]){
    TH1F *test = new TH1F("test", "test", 100, 10, 700);
    TCanvas *cst = new TCanvas("cst", "TH1F histograms", 10, 10, 700, 700);

    int numberOfIterations = 0;
    for(int i = 0; i < nSamples; i++){
        numberOfIterations = max(numberOfIterations,sampleLength[i]);
    }
    numberOfIterations = min(numberOfIterations, maxEvents);

    for(int i = 0; i < nSamples; i++){
        for(int j = 0; j < numberOfIterations; j++){
        sampleTrees[i]->GetEntry(j);
            for(int k = 0; k < nDoubleVariables; k++){
                TH1F_histograms[i][k].Fill(doubleVariables[i][k][doubleArrayIndex[k]]);
                test.Fill(doubleVariables[i][k][doubleArrayIndex[k]]);
            }
        }
    }
    test.Draw();
    cst->SaveAs("test.png");
}

void drawTH1FHistograms(TH1F **TH1F_histograms,const std::string sampleNames[], const int nSamples, const std::string doubleVariableNames[], const int nDoubleVariables){
            std::string fileName;
            TCanvas *cst = new TCanvas("cst", "TH1F histograms", 10, 10, 700, 700);
            cst->SetLogy();

            for(int i = 0; i < nSamples; i++){
                for(int j = 0; j < nDoubleVariables; j++){
                    fileName = "TH1F_" + sampleNames[i] + "_" + doubleVariableNames[j] + ".png";
                    TH1F_histograms[i][j].Draw();
                    cst->SaveAs(fileName.c_str());  
                }
            }
}

void correlationStudy(){
    
    ///constants
    const std::string sampleNames[] = {"MonoTop", "WJet", "TTbar"};
    const int nSamples = sizeof(sampleNames)/sizeof(sampleNames[0]);

    const std::string doubleVariableNames[] =   {"MetPt",   "PFMuonPt"};
    const int minValue[] =                      {10,        10};
    const int maxValue[] =                      {700,       700};
    const int doubleArraySize[] =               {40,        13};
    const int doubleArrayIndex[] =              {10,        0};

    const int nBins = 100;
    const int nDoubleVariables = sizeof(doubleVariableNames)/sizeof(doubleVariableNames[0]);
    const int maxEvents = 5000; 

    ///variables
    int sampleLength[nSamples];
    double ***doubleVariables = new double**[nSamples];

    //root objects
    TTree *sampleTrees[nSamples];

    TH1F **TH1F_histograms = new TH1F*[nSamples];

    

    loadRootFiles(sampleNames, nSamples, sampleTrees, sampleLength);

    makeTH1FHistograms(sampleNames, nSamples, doubleVariableNames, nDoubleVariables, nBins, minValue, maxValue,  TH1F_histograms);

    getBranches(doubleVariableNames, nDoubleVariables, doubleArraySize, doubleVariables, sampleTrees, nSamples);    

    plotTH1FHistograms(sampleTrees, nSamples, doubleVariables,  nDoubleVariables, TH1F_histograms, sampleLength, maxEvents, doubleArrayIndex);

    drawTH1FHistograms(TH1F_histograms, sampleNames, nSamples, doubleVariableNames, nDoubleVariables);

}
