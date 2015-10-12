#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <algorithm>
#include <vector>
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

void makeTH1FHistograms(const std::string sampleNames[], const int nSamples,  const std::string doubleVariableNames[],const int nDoubleVariables, const int nBins, const int minDoubleValue[], const int maxDoubleValue[], std::vector < std::vector< TH1F* > > TH1F_Vector){
//void makeTH1FHistograms(const std::string sampleNames[], const int nSamples,  const std::string doubleVariableNames[],const int nDoubleVariables, const int nBins, const int minDoubleValue[], const int maxDoubleValue[], TH1F *TH1F_histograms[][2]){
//
////    for (int i = 0; i < nSamples; i++){
////        TH1F_histograms[i] = new (TH1F*) [nDoubleVariables];
////    }
//
    std::string histName;
    for(int i = 0; i < nSamples; i++){
        for(int j = 0; j < nDoubleVariables; j++){
            histName = "TH1F_" + sampleNames[i] + "_" + doubleVariableNames[j];
            std::cout << histName << std::endl;
            TH1F_Vector[i].push_back( new TH1F(histName.c_str(),histName.c_str(),nBins, minDoubleValue[j], maxDoubleValue[j]));                
        }
    } 

    for(int i = 0; i < nSamples; i++){
            for(int k = 0; k < nDoubleVariables; k++){
                TH1F_Vector[i][k]->Fill(rand() % 100); //doubleVariables[i][k][doubleArrayIndex[k]]);
            }
        }

//
//    
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


void plotTH1FHistograms(TTree *sampleTrees[], const int nSamples, double ***doubleVariables, const int nDoubleVariables, std::vector < std::vector< TH1F* > > TH1F_Vector, int  sampleLength[], const int maxEvents, const int doubleArrayIndex[]){

    int numberOfIterations = 0;
    for(int i = 0; i < nSamples; i++){
        numberOfIterations = max(numberOfIterations,sampleLength[i]);
    }
    numberOfIterations = min(numberOfIterations, maxEvents);

    for(int i = 0; i < nSamples; i++){
        for(int j = 0; j < min(numberOfIterations,sampleLength[i]); j++){
        sampleTrees[i]->GetEntry(j);
            for(int k = 0; k < nDoubleVariables; k++){
                //TH1F_Vector[i][k]->Fill(rand() % 100); //doubleVariables[i][k][doubleArrayIndex[k]]);
                TH1F_Vector[i][k]->Fill(doubleVariables[i][k][doubleArrayIndex[k]]);
                //std::cout << doubleVariables[i][k][doubleArrayIndex[k]] << std::endl;
            }
        }
    }
}

void plotTH2FHistograms(TTree *sampleTrees[], const int nSamples, double ***doubleVariables, const int nDoubleVariables,std::vector < std::vector < std::vector< TH2F* > > > TH2F_Vector, int  sampleLength[], const int maxEvents, const int doubleArrayIndex[]){

    int numberOfIterations = 0;
    for(int i = 0; i < nSamples; i++){
        numberOfIterations = max(numberOfIterations,sampleLength[i]);
    }
    numberOfIterations = min(numberOfIterations, maxEvents);

    for(int i = 0; i < nSamples; i++){
        for(int j = 0; j < min(numberOfIterations, sampleLength[i]); j++){
        sampleTrees[i]->GetEntry(j);
            for(int k = 0; k < nDoubleVariables; k++){
                for( int l = 0; l < nDoubleVariables; l++){
                    if(k < l){
                        TH2F_Vector[i][k][l]->Fill(doubleVariables[i][k][doubleArrayIndex[k]], doubleVariables[i][l][doubleArrayIndex[l]]);
                    }
                }
            }
        }
    }
}




void sqrtTH1F(TH1F *histogram){
    double errorOfBin;
    for(int i = 0; i < histogram->GetSize(); i++){
        errorOfBin = histogram->GetBinError(i);
        histogram->SetBinContent(i, errorOfBin);    
    }
}

void sqrtTH2F(TH2F *histogram){
    double errorOfBin;
    for(int i = 0; i < histogram->GetXaxis()->GetNbins(); i++){
        for(int j = 0; j < histogram->GetYaxis()->GetNbins(); j++){

        errorOfBin = histogram->GetBinError(i,j);
        histogram->SetBinContent(i,j, errorOfBin);    

        }
    }
}

void integrateTH1F(TH1F *histogram){
        double integral;
        for(int i = 0; i < histogram->GetXaxis()->GetNbins(); i++){
            integral = histogram->Integral(i,histogram->GetXaxis()->GetNbins());
            histogram->SetBinContent(i, integral);
        }
}

void drawTH1FHistograms(std::vector < std::vector< TH1F* > > TH1F_Vector, const std::string sampleNames[], const int nSamples, const std::string doubleVariableNames[], const int nDoubleVariables){
            std::string fileName;
            TCanvas *cst = new TCanvas("cst", "TH1F histograms", 10, 10, 700, 700);
            cst->SetLogy();
            TH1F *SandB;
            TH1F *Signal;
            Double_t Norm[] = {.000007207, 61.5267, .83176};


            for(int i = 0; i < nSamples; i++){
                for(int j = 0; j < nDoubleVariables; j++){
                    fileName = "TH1F_" + sampleNames[i] + "_" + doubleVariableNames[j] + ".png";
                    TH1F_Vector[i][j]->Scale(Norm[i]);
                    TH1F_Vector[i][j]->Draw();
                    cst->SaveAs(fileName.c_str());  
                }
            }

            for(int j = 0; j < nDoubleVariables; j++){
                for(int i = 1; i < nSamples; i++){
                    fileName = "TH1F_" + sampleNames[i] + "_" + doubleVariableNames[j] +"_significance.png";
                    Signal = (TH1F*)TH1F_Vector[0][j]->Clone();       
                    Signal->SetTitle(fileName.c_str());
                    SandB = (TH1F*)TH1F_Vector[0][j]->Clone();       
                    SandB->SetTitle(fileName.c_str());
                    SandB->Add(TH1F_Vector[i][j]);
                    integrateTH1F(Signal);
                    integrateTH1F(SandB);
                    sqrtTH1F(SandB);
                    Signal->Divide(SandB);
                    Signal->Draw();

                    //TH1F_Vector[i][j]->Add(TH1F_Vector[0][j]);
                    //sqrtTH1F(TH1F_Vector[i][j]);
                    //TH1F_Vector[0][j]->Divide(TH1F_Vector[i][j]);
                    //TH1F_Vector[0][j]->SetTitle(fileName.c_str());
                    //TH1F_Vector[0][j]->Draw();
                    cst->SaveAs(fileName.c_str());
                    //TH1F_Vector[0][j]->Multiply(TH1F_Vector[i][j]);
                }
            }
        

}

void drawTH2FHistograms(std::vector < std::vector < std::vector< TH2F* > > > TH2F_Vector, const std::string sampleNames[], const int nSamples, const std::string doubleVariableNames[], const int nDoubleVariables){
            std::string fileName;
            TCanvas *cst2F = new TCanvas("cst", "TH2F histograms", 10, 10, 700, 700);
            Double_t Norm[] = {.000007207, 61.5267, .83176};

            for(int i = 0; i < nSamples; i++){
                for(int j = 0; j < nDoubleVariables; j++){
                    for(int k = 0; k< nDoubleVariables; k++){
                        if(j < k){
                            fileName = "TH2F_" + sampleNames[i] + "_" + doubleVariableNames[j] + "_" + doubleVariableNames[k] + ".png";
                            TH2F_Vector[i][j][k]->Scale(Norm[i]);
                            TH2F_Vector[i][j][k]->Draw("colz");
                            cst2F->SaveAs(fileName.c_str());  
                        }
                    }
                }
            }

            for(int j = 0; j < nDoubleVariables; j++){
                for(int k = 0; k < nDoubleVariables; k++){
                    for(int i = 1; i < nSamples; i++){
                        if(j < k){
                            TH2F_Vector[i][j][k]->Add(TH2F_Vector[0][j][k]);
                            sqrtTH2F(TH2F_Vector[i][j][k]);
                            TH2F_Vector[0][j][k]->Divide(TH2F_Vector[i][j][k]);
                            fileName = "TH2F_" + sampleNames[i] + "_" + doubleVariableNames[j] + "_" + doubleVariableNames[k] + "_significance.png";
                            TH2F_Vector[0][j][k]->SetTitle(fileName.c_str());
                            TH2F_Vector[0][j][k]->Draw("colz");
                            cst2F->SaveAs(fileName.c_str());
                            TH2F_Vector[0][j][k]->Multiply(TH2F_Vector[i][j][k]);
                        }
                    }
                }
            }



}

void correlationStudy(){
    
    ///constants
    const std::string sampleNames[] = {"MonoTop", "WJet", "TTbar"};
    const int nSamples = sizeof(sampleNames)/sizeof(sampleNames[0]);

    const std::string intVariableNames[] =      {"NPFMuon"};
    const int minIntValue[] =                   {1}; 
    const int maxIntValue[] =                   {13};
    const int intArraySize[] =                  {1};
    const int intArrayIndex[] =                 {0};

    const std::string doubleVariableNames[] =   {"MetPt",   "PFMuonPt", "MetPhi",   "PFMuonPhi",    "MetSumEt",     "PFMuonEta",    "PFMuonChargedHadronIso"};
    const int minDoubleValue[] =                {10,        10,         -4,         -4,             10,             -3,             0};
    const int maxDoubleValue[] =                {700,       700,        4,          4,              4000,           3,              10};
    const int doubleArraySize[] =               {40,        13,         40,         13,             40,             13,             13};
    const int doubleArrayIndex[] =              {10,        0,          10,         0,              10,             0,              0};

    const int nBins = 100;
    const int nDoubleVariables = sizeof(doubleVariableNames)/sizeof(doubleVariableNames[0]);
    const int maxEvents = 10000; 

    ///variables
    int sampleLength[nSamples];
    double ***doubleVariables = new double**[nSamples];
    std::string histName;
    std::string hist2DName;

    //root objects
    TTree *sampleTrees[nSamples];

    //TH1F **TH1F_histograms = new TH1F*[nSamples];
    TH1F *TH1F_histograms[nSamples][nDoubleVariables];
    std::vector < std::vector< TH1F* > > TH1F_Vector;
    TH1F_Vector.resize(nSamples);

    std::vector < std::vector < std::vector < TH2F* > > > TH2F_Vector;
    TH2F_Vector.resize(nSamples);
    
 //   TH1F **ptrTo_TH1F_histograms = TH1F_histograms;

    //TH1F_histograms[0][0] = new TH1F("test","test",100,10, 100);
    

    loadRootFiles(sampleNames, nSamples, sampleTrees, sampleLength);

    for(int i = 0; i < nSamples; i++){
        TH2F_Vector[i].resize(nDoubleVariables);
        for(int j = 0; j < nDoubleVariables; j++){
            histName = "TH1F_" + sampleNames[i] + "_" + doubleVariableNames[j];
            std::cout << histName << std::endl;
            TH1F_Vector[i].push_back( new TH1F(histName.c_str(),histName.c_str(),nBins, minDoubleValue[j], maxDoubleValue[j]));                
            for(int k = 0; k < nDoubleVariables; k++){
                hist2DName = "TH2F_" + sampleNames[i] + "_" + doubleVariableNames[j] + "_" + doubleVariableNames[k];
                std::cout << hist2DName << std::endl;
                TH2F_Vector[i][j].push_back(new TH2F(hist2DName.c_str(), hist2DName.c_str(), nBins, minDoubleValue[j], maxDoubleValue[j], nBins, minDoubleValue[k], maxDoubleValue[k]));
            }            
        }
    } 
    

    //makeTH1FHistograms(sampleNames, nSamples, doubleVariableNames, nDoubleVariables, nBins, minDoubleValue, maxDoubleValue,  TH1F_Vector);

    getBranches(doubleVariableNames, nDoubleVariables, doubleArraySize, doubleVariables, sampleTrees, nSamples);    

    plotTH1FHistograms(sampleTrees, nSamples, doubleVariables,  nDoubleVariables, TH1F_Vector, sampleLength, maxEvents, doubleArrayIndex);

    plotTH2FHistograms(sampleTrees, nSamples, doubleVariables,  nDoubleVariables, TH2F_Vector, sampleLength, maxEvents, doubleArrayIndex);

    drawTH1FHistograms(TH1F_Vector, sampleNames, nSamples, doubleVariableNames, nDoubleVariables);

    drawTH2FHistograms(TH2F_Vector, sampleNames, nSamples, doubleVariableNames, nDoubleVariables);

}
