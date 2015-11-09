#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <algorithm>
#include <vector>
//#include <TTreeReader.h>
//#include <TTreeReaderValue.h>

    struct cuts{
        std::string cutVariable;
        int count;
    };

double transverseMass(Double_t PFMuonPhiTree[13],Double_t MetPhiTree[40], Double_t PFMuonPtTree[13], Double_t MetPtTree[40]){
            double transverseMassVariable;

            transverseMassVariable = cos(PFMuonPhiTree[0]-MetPhiTree[10]);     
            transverseMassVariable = 1 - transverseMassVariable;
            transverseMassVariable = 2*PFMuonPtTree[0]*MetPtTree[10]*transverseMassVariable;
            transverseMassVariable = sqrt(transverseMassVariable);

            return transverseMassVariable;
}








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

void makeTH1FHistograms(const std::string sampleNames[], const int nSamples,  const std::string doubleVariableNames[],const int nDoubleVariables, const int nBins, const double minDoubleValue[], const double maxDoubleValue[], std::vector < std::vector< TH1F* > > TH1F_Vector){
//void makeTH1FHistograms(const std::string sampleNames[], const int nSamples,  const std::string doubleVariableNames[],const int nDoubleVariables, const int nBins, const double minDoubleValue[], const double maxDoubleValue[], TH1F *TH1F_histograms[][2]){
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

void getIntBranches(const std::string intVariableNames[], const int nIntVariables, const int intArraySize[], int  ***intVariables, TTree *sampleTrees[], const int nSamples){


    for(int i = 0; i < nSamples; i++){
        intVariables[i] = new int*[nIntVariables];
        for(int j = 0; j < nIntVariables; j++){
            //doubleVariables[i][j] = new double [doubleArraySize[j]]; should work but doesn't? Look into with time. 
            intVariables[i][j] = new int [100];
            sampleTrees[i]->SetBranchAddress(intVariableNames[j].c_str(), &*intVariables[i][j]);
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
bool postCut(){


}

double getDouble(int i, int variable_index, const int doubleArrayIndex[], double ***doubleVariables){
    
    return doubleVariables[i][variable_index][doubleArrayIndex[variable_index]];
}

bool precut(int i,const int nDoubleVariables, double ***doubleVariables,const int doubleArrayIndex[], int ***intVariables, const int nIntVariables, const int intArrayIndex[], const double minDoubleCutValue[], const double maxDoubleCutValue[], const int minIntCutValue[], const int maxIntCutValue[], std::vector < std::vector< TH1F* > > TH1F_cutRatios, const std::string doubleVariableNames[], const std::string intVariableNames[], std::vector < std::vector< cuts* > > cutFlow){
    bool noCut = true;
    TH1F_cutRatios[i][0]->Fill(-1);
    int PFMuonPt_index;
    int PFMuonChargedHadronIso_index;
    int PFMuonPhotonIso_index;
    int PFMuonNeutralHadronIso_index;
    int PFMuonEta_index;
    int NPFMuon_index;
    int MetPhi_index, PFAK4JetPhi_index, PFMuonPhi_index, MetSumEt_index, PFAK4JetPt_index, MetPt_index;

 //std::cout << "____________------------------" << std::endl;
    cutFlow[i][0]->cutVariable = "All events: "; 
    cutFlow[i][0]->count = cutFlow[i][0]->count + 1;


    for(int k = 0; k < nDoubleVariables; k++){
        if((doubleVariables[i][k][doubleArrayIndex[k]] > minDoubleCutValue[k] and doubleVariables[i][k][doubleArrayIndex[k]] < maxDoubleCutValue[k]) and noCut == true){
            //std::cout << "double " << k << " " << doubleVariables[i][k][doubleArrayIndex[k]] <<  " not in range " << minDoubleCutValue[k] << " " << maxDoubleCutValue[k] << std::endl;
       //     std::cout << doubleVariableNames[k] << " " << noCut << std::endl;
            noCut = true;
            cutFlow[i][k + 1]->cutVariable = doubleVariableNames[k] + " " + std::to_string(minDoubleCutValue[k]) + " to " + std::to_string(maxDoubleCutValue[k]) + ": ";
            //std::cout << cutFlow[i][k]->cutVariable << std::endl;
            cutFlow[i][k + 1]->count = cutFlow[i][k + 1]->count + 1;
            TH1F_cutRatios[i][0]->Fill(k + nIntVariables);
        } else{
            noCut = false;
     //       std::cout << "false" <<std::endl;
        } 

        if(doubleVariableNames[k] == "PFMuonPt"){
            PFMuonPt_index = k; 
        }
        if(doubleVariableNames[k] == "PFMuonChargedHadronIso"){
            PFMuonChargedHadronIso_index = k; 
        }
        if(doubleVariableNames[k] == "PFMuonPhotonIso"){
            PFMuonPhotonIso_index = k; 
        }
        if(doubleVariableNames[k] == "PFMuonNeutralHadronIso"){
            PFMuonNeutralHadronIso_index = k; 
        }
        if(doubleVariableNames[k] == "PFMuonEta"){
            PFMuonEta_index = k; 
        }
        if( doubleVariableNames[k] == "MetPhi"){
            MetPhi_index = k;
        }
        if( doubleVariableNames[k] == "PFAK4JetPhi"){
            PFAK4JetPhi_index = k;
        } 
        if( doubleVariableNames[k] == "PFMuonPhi"){
            PFMuonPhi_index = k;
        } 
        if( doubleVariableNames[k] == "MetSumEt"){
            MetSumEt_index = k;
        } 
        if( doubleVariableNames[k] == "PFAK4JetPt"){
            PFAK4JetPt_index = k;
        } 
        if( doubleVariableNames[k] == "MetPt"){
            MetPt_index = k;
        } 
    }

    for(int k = 0; k < nIntVariables; k++){
           // std::cout << i <<  intVariableNames[k] << " " << intVariables[i][k][intArrayIndex[k]] << std::endl;
        if((intVariables[i][k][intArrayIndex[k]] >= minIntCutValue[k] and intVariables[i][k][intArrayIndex[k]] <= maxIntCutValue[k]) and noCut == true){
            //std::cout << "int "  << k << " " <<  intVariables[i][k][intArrayIndex[k]] <<  " not in range " << minIntCutValue[k] << " " << maxIntCutValue[k] << std::endl;
            noCut = true; 
            cutFlow[i][k + 1 + nDoubleVariables ]->cutVariable = intVariableNames[k] + " " + std::to_string(minIntCutValue[k]) + " to " + std::to_string(maxIntCutValue[k]);
            //std::cout << cutFlow[i][k]->cutVariable << std::endl;
            cutFlow[i][k + 1  + nDoubleVariables]->count = cutFlow[i][k + nDoubleVariables + 1]->count + 1;
            TH1F_cutRatios[i][0]->Fill(k);
        } else {
            noCut = false;
           // std::cout << "false" <<std::endl;
        } 
        if(intVariableNames[k] == "NPFMuon"){
            NPFMuon_index = k;
        }

        





    }

        if((intVariables[i][NPFMuon_index][intArrayIndex[NPFMuon_index]] > 1) and noCut == true){
          //std::cout << intVariables[i][k][intArrayIndex[k]] << std::endl;
            for(int l = 1; l <= intVariables[i][NPFMuon_index][intArrayIndex[NPFMuon_index]]; l++){
                if(not (doubleVariables[i][PFMuonPt_index][l] > 10 and abs(doubleVariables[i][PFMuonEta_index][l]) < 2.4 and doubleVariables[i][PFMuonChargedHadronIso_index][l] < .2 and doubleVariables[i][PFMuonPhotonIso_index][l] < .2 and doubleVariables[i][PFMuonNeutralHadronIso_index][l] < .2)){
   //                 std::cout << i << " true  muon pt " << doubleVariables[i][PFMuonPt_index][l] << " " << abs(doubleVariables[i][PFMuonEta_index][l]) << " hadron " << doubleVariables[i][PFMuonChargedHadronIso_index][l] << " photon " << doubleVariables[i][PFMuonPhotonIso_index][l] << " neutral " << doubleVariables[i][PFMuonNeutralHadronIso_index][l] << std::endl;
                    noCut = true;
                    if(l == intVariables[i][NPFMuon_index][intArrayIndex[NPFMuon_index]]-1){
                    }
        //            std::cout << "Other isolated lepton" << std::endl;
                    TH1F_cutRatios[i][0]->Fill(nIntVariables+nDoubleVariables);
                } else{
     //               std::cout << i << " false muon pt " << doubleVariables[i][PFMuonPt_index][l] << " " << abs(doubleVariables[i][PFMuonEta_index][l]) << " hadron " << doubleVariables[i][PFMuonChargedHadronIso_index][l] << " photon "  << doubleVariables[i][PFMuonPhotonIso_index][l] << " neutral " << doubleVariables[i][PFMuonNeutralHadronIso_index][l] << std::endl;
                    noCut = false;
                    cutFlow[i][nDoubleVariables + nIntVariables + 1]->cutVariable = "No other isolated leptons: "; 
                    cutFlow[i][nDoubleVariables + nIntVariables + 1]->count = cutFlow[i][nDoubleVariables + nIntVariables + 1]->count + 1;
                    //std::cout << "false" <<std::endl;
                }
            } 
        }

    ///special cuts:
        //if( (abs((abs(doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - doubleVariables[i][PFMuonPhi_index][doubleArrayIndex[PFMuonPhi_index]])-3.14159265359)+(abs(doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - doubleVariables[i][PFAK4JetPhi_index][doubleArrayIndex[PFAK4JetPhi_index]])-3.14159265359)) < 1) and noCut == true){
        //if( (abs((abs(doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - doubleVariables[i][PFMuonPhi_index][doubleArrayIndex[PFMuonPhi_index]])-3.14159265359)+(abs(doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - doubleVariables[i][PFAK4JetPhi_index][doubleArrayIndex[PFAK4JetPhi_index]])-3.14159265359)) < 1) and noCut == true){
        if(1==1){//abs(abs(doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - doubleVariables[i][PFMuonPhi_index][doubleArrayIndex[PFMuonPhi_index]]) -3.14159265359) < 1 and noCut == true){
        //if((abs(doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - doubleVariables[i][PFMuonPhi_index][doubleArrayIndex[PFMuonPhi_index]])-3.14159265359) < 1 and noCut == true){
                    //std::cout<< abs((abs(doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - doubleVariables[i][PFMuonPhi_index][doubleArrayIndex[PFMuonPhi_index]])-3.14159265359)+(abs(doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - doubleVariables[i][PFAK4JetPhi_index][doubleArrayIndex[PFAK4JetPhi_index]])-3.14159265359))  << " makes sense?" << std::endl;
                    cutFlow[i][nDoubleVariables + nIntVariables + 2]->cutVariable = "MetPhi, MuonPhi, and Jet Phi cut : "; 
                    cutFlow[i][nDoubleVariables + nIntVariables + 2]->count = cutFlow[i][nDoubleVariables + nIntVariables + 2]->count + 1;
                }else{
                noCut = false;
                }

         if(1==1){//doubleVariables[i][PFAK4JetPt_index][doubleArrayIndex[PFAK4JetPt_index]] - doubleVariables[i][MetSumEt_index][doubleArrayIndex[MetSumEt_index]]/2 > -130 and  doubleVariables[i][PFAK4JetPt_index][doubleArrayIndex[PFAK4JetPt_index]]/doubleVariables[i][MetSumEt_index][doubleArrayIndex[MetSumEt_index]] > .21 and doubleVariables[i][MetSumEt_index][doubleArrayIndex[MetSumEt_index]] < 900 and noCut == true){
                    //std::cout << PFAK4JetPt_index << " " << i << " "  << MetSumEt_index << std::endl;
                    cutFlow[i][nDoubleVariables + nIntVariables + 3]->cutVariable = "PKA4JetPt/MetsumEt : "; 
                    cutFlow[i][nDoubleVariables + nIntVariables + 3]->count = cutFlow[i][nDoubleVariables + nIntVariables + 3]->count + 1;

         }   else{

                noCut = false;
         }    
        std::cout << getDouble(i, MetPt_index, doubleArrayIndex, doubleVariables)/getDouble(i, PFAK4JetPt_index, doubleArrayIndex, doubleVariables) << " " << i << std::endl; 
         if(1==1){//getDouble(i, MetPt_index, doubleArrayIndex, doubleVariables)/getDouble(i, PFAK4JetPt_index, doubleArrayIndex, doubleVariables) > .7 and getDouble(i, MetPt_index, doubleArrayIndex, doubleVariables) > 100 and  noCut == true){
        std::cout << getDouble(i, MetPt_index, doubleArrayIndex, doubleVariables)/getDouble(i, PFAK4JetPt_index, doubleArrayIndex, doubleVariables) << " " << i << std::endl; 
                    cutFlow[i][nDoubleVariables + nIntVariables + 4]->cutVariable = "MetPt v PFAk4JetPt : "; 
                    cutFlow[i][nDoubleVariables + nIntVariables + 4]->count = cutFlow[i][nDoubleVariables + nIntVariables + 4]->count + 1;

         }   else{

                noCut = false;
        }
         if(1==1){//getDouble(i, MetSumEt_index, doubleArrayIndex, doubleVariables)-40.0/17.0*getDouble(i, MetPt_index, doubleArrayIndex, doubleVariables) <3700.0/17.0 and  noCut == true){
           std::cout << getDouble(i, MetPt_index, doubleArrayIndex, doubleVariables)/getDouble(i, PFAK4JetPt_index, doubleArrayIndex, doubleVariables) << " " << i << std::endl; 
                    cutFlow[i][nDoubleVariables + nIntVariables + 5]->cutVariable = "METSumEt MetPt cut : ";
                    cutFlow[i][nDoubleVariables + nIntVariables + 5]->count = cutFlow[i][nDoubleVariables + nIntVariables + 5]->count + 1;

         }   else{

                noCut = false;
        }


    return noCut;
}

void sqrtTH1F(TH1F *histogram){
    double errorOfBin;
    for(int i = 0; i < histogram->GetSize(); i++){
        errorOfBin = histogram->GetBinError(i);
        histogram->SetBinContent(i, errorOfBin);    
    }
}

void integrateCenterTH1F(int nSamples, std::vector<TH1F*> histogram, TCanvas *cst){
        double integral;
        int range; 
        //std::string name;
        std::vector<TH1F*> tempHistogram;
        for(int i = 0; i < nSamples; i++){
            std::cout << nSamples << std::endl;
            range =histogram[i]->GetXaxis()->GetNbins()/2;
            tempHistogram.push_back( (TH1F*)histogram[i]->Clone());
            for(int j = 0; j < range; j++){
                integral = histogram[i]->Integral(j,histogram[i]->GetXaxis()->GetNbins()-j);
                tempHistogram[i]->SetBinContent(j, integral);
                tempHistogram[i]->SetBinContent(histogram[i]->GetXaxis()->GetNbins()-j, integral);
            }
            

            if (i > 0){
                //tempHistogram[i]->Add(tempHistogram[0]);
                sqrtTH1F(tempHistogram[i]);
            }
        }
        for(int i = 0; i < nSamples; i++){
            if(i>0){
            tempHistogram[0]->Divide(tempHistogram[i]);
            }

            tempHistogram[0]->Draw();
            string name(tempHistogram[i]->GetName(),0,100);
            name = name + "_sig.png";
            cst->SaveAs(name.c_str());

            if(i>0){
            tempHistogram[0]->Multiply(tempHistogram[i]);
            }
            
        }

}




void plotTH1FHistograms(TTree *sampleTrees[], const int nSamples, double ***doubleVariables, const int nDoubleVariables, std::vector < std::vector< TH1F* > > TH1F_Vector, int  sampleLength[], const int maxEvents, const int doubleArrayIndex[], int ***intVariables, const int nIntVariables, const int intArrayIndex[], const double minDoubleCutValue[], const double maxDoubleCutValue[], const int minIntCutValue[], const int maxIntCutValue[], std::vector < std::vector< TH1F* > > TH1F_cutRatios, const std::string doubleVariableNames[], const std::string intVariableNames[], std::vector < std::vector< cuts* > > cutFlow){

    std::vector <TH1F*> MetPhiMinusPFAK4JetPhi;
    std::vector <TH1F*> TH1F_transverseMass;
    std::vector <TH2F*> DeltaMuonPhiVDeltaJetPhi;
    std::string name;
    std::vector <std::vector < TH2F*>> TH2F_transverseMass_v;
    TH2F_transverseMass_v.resize(nSamples);
    double MetPhiVariable, PFAK4JetPhiVariable, x_variable, y_variable;
    int MetPhi_index, PFAK4JetPhi_index, PFMuonPhi_index, MetPt_index, PFMuonPt_index;

    for(i = 0; i < nSamples; i++){
        name = "MetPhiMinusPFAK4JetPhi_" + std::to_string(i);   
        MetPhiMinusPFAK4JetPhi.push_back(new TH1F(name.c_str(), name.c_str(), 100, -10, 10));
        name = "TH1F_transverseMass_" + std::to_string(i);   
        TH1F_transverseMass.push_back(new TH1F(name.c_str(), name.c_str(), 100, -10, 700));
        name = "DeltaMuonPhiVDeltaJetPhi_" + std::to_string(i);
        DeltaMuonPhiVDeltaJetPhi.push_back(new TH2F(name.c_str(), name.c_str(), 100, -7,7, 100, -7,7)); 
        for(int k = 0; k < nDoubleVariables; k++){
             name = "TH2F_transverseMass_v_" +doubleVariableNames[k] + "_"+ std::to_string(i);
             TH2F_transverseMass_v[i].push_back(new TH2F(name.c_str(), name.c_str(), 100, -7,700, 100, -7,700));
        }

    }

    int numberOfIterations = 0;
    for(int i = 0; i < nSamples; i++){
        numberOfIterations = max(numberOfIterations,sampleLength[i]);
    }
    numberOfIterations = min(numberOfIterations, maxEvents);
    for(int i = 0; i < nSamples; i++){
        for(int j = 0; j < min(numberOfIterations,sampleLength[i]); j++){
        sampleTrees[i]->GetEntry(j);
        //std::cout << i << " " << intVariableNames[0]<< " "  << intVariables[i][0][intArrayIndex[0]] << " " << intVariableNames[1] << " " << intVariables[i][1][intArrayIndex[1]] << std::endl; 
        //precut(i, nDoubleVariables, doubleVariables, doubleArrayIndex, intVariables, nIntVariables, intArrayIndex, minDoubleCutValue,  maxDoubleCutValue, minIntCutValue,  maxIntCutValue, TH1F_cutRatios, doubleVariableNames, intVariableNames, cutFlow);
            if(precut(i, nDoubleVariables, doubleVariables, doubleArrayIndex, intVariables, nIntVariables, intArrayIndex, minDoubleCutValue,  maxDoubleCutValue, minIntCutValue,  maxIntCutValue, TH1F_cutRatios, doubleVariableNames, intVariableNames, cutFlow)){

                for(int k = 0; k < nDoubleVariables; k++){
                    if( doubleVariableNames[k] == "MetPhi"){
                        MetPhiVariable = doubleVariables[i][k][doubleArrayIndex[k]];
                        MetPhi_index = k;
                    }
                    if( doubleVariableNames[k] == "PFAK4JetPhi"){
                        PFAK4JetPhiVariable = doubleVariables[i][k][doubleArrayIndex[k]];
                        PFAK4JetPhi_index = k;
                    } 
                    if( doubleVariableNames[k] == "PFMuonPhi"){
                        PFMuonPhi_index = k;
                    } 
                    if( doubleVariableNames[k] == "MetPt"){
                        MetPt_index = k;
                    } 
                    if( doubleVariableNames[k] == "PFMuonPt"){
                        PFMuonPt_index = k;
                    }    
                } 
                for(int k = 0; k < nDoubleVariables; k++){

                    //TH1F_Vector[i][k]->Fill(rand() % 100); //doubleVariables[i][k][doubleArrayIndex[k]]);
                    TH1F_Vector[i][k]->Fill(doubleVariables[i][k][doubleArrayIndex[k]]);
                    TH2F_transverseMass_v[i][k]->Fill(transverseMass(doubleVariables[i][PFMuonPhi_index], doubleVariables[i][MetPhi_index], doubleVariables[i][PFMuonPt_index], doubleVariables[i][MetPt_index]) ,doubleVariables[i][k][doubleArrayIndex[k]]);
                    //std::cout << doubleVariables[i][k][doubleArrayIndex[k]] << std::endl;
                    //

                }   
//                    std::cout << MetPhi_index <<std::endl;
//                   std::cout << PFAK4JetPhi_index <<std::endl;
//                    std::cout << doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] << std::endl;
//                    std::cout << doubleVariables[i][PFAK4JetPhi_index][doubleArrayIndex[PFAK4JetPhi_index]] << std::endl;
//                    std::cout << doubleVariables[i][9][doubleArrayIndex[9]] << std::endl;
                            
    std::cout << "transverse mass: " <<transverseMass(doubleVariables[i][PFMuonPhi_index], doubleVariables[i][MetPhi_index], doubleVariables[i][PFMuonPt_index], doubleVariables[i][MetPt_index]) << std::endl;
                            TH1F_transverseMass[i]->Fill(transverseMass(doubleVariables[i][PFMuonPhi_index], doubleVariables[i][MetPhi_index], doubleVariables[i][PFMuonPt_index], doubleVariables[i][MetPt_index]));

                            //MetPhiMinusPFAK4JetPhi[i]->Fill(abs(doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - (doubleVariables[i][PFAK4JetPhi_index][doubleArrayIndex[PFAK4JetPhi_index]]+doubleVariables[i][PFMuonPhi_index][doubleArrayIndex[PFMuonPhi_index]])/2)-3.14159265359);
                              MetPhiMinusPFAK4JetPhi[i]->Fill((doubleVariables[i][PFMuonPhi_index][doubleArrayIndex[PFMuonPhi_index]] - doubleVariables[i][PFAK4JetPhi_index][doubleArrayIndex[PFAK4JetPhi_index]]));
                            if(1==1){// doubleVariables[i][MetPt_index][doubleArrayIndex[MetPt_index]] < 277){
                            x_variable = doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - doubleVariables[i][PFMuonPhi_index][doubleArrayIndex[PFMuonPhi_index]];
                            y_variable = doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - doubleVariables[i][PFAK4JetPhi_index][doubleArrayIndex[PFAK4JetPhi_index]];

                              DeltaMuonPhiVDeltaJetPhi[i]->Fill(x_variable*cos(3.1415926*.25) - y_variable*sin(3.1415926*.25), x_variable*sin(3.1415926*.25) + y_variable*cos(3.1415926*.25));
                              //DeltaMuonPhiVDeltaJetPhi[i]->Fill(x_variable, y_variable);
                            }
                              //MetPhiMinusPFAK4JetPhi[i]->Fill((abs(doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - doubleVariables[i][PFMuonPhi_index][doubleArrayIndex[PFMuonPhi_index]])-3.14159265359));
                   //           MetPhiMinusPFAK4JetPhi[i]->Fill((abs(doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - doubleVariables[i][PFMuonPhi_index][doubleArrayIndex[PFMuonPhi_index]])-3.14159265359)+(abs(doubleVariables[i][MetPhi_index][doubleArrayIndex[MetPhi_index]] - doubleVariables[i][PFAK4JetPhi_index][doubleArrayIndex[PFAK4JetPhi_index]])-3.14159265359));
            }
        }
    }
    TCanvas *cst = new TCanvas("cst", "TH1F histograms", 10, 10, 700, 700);
    for(i = 0; i < nSamples; i++){
        name = "MetPhiMinusPFAK4JetPhi_" + std::to_string(i) + ".png";
        MetPhiMinusPFAK4JetPhi[i]->Draw();
        cst->SaveAs(name.c_str());
        name = "TH1F_transverseMass_" + std::to_string(i) + ".png";
        TH1F_transverseMass[i]->Draw();
        cst->SaveAs(name.c_str());
         name = "DeltaMuonPhiVDeltaJetPhi_" + std::to_string(i) + ".png";
         DeltaMuonPhiVDeltaJetPhi[i]->Draw("colz");
         cst->SaveAs(name.c_str());
        for(int k = 0; k < nDoubleVariables; k++){
             name = "TH2F_transverseMass_v_" +doubleVariableNames[k] + "_"+ std::to_string(i)+".png";
             TH2F_transverseMass_v[i][k]->Draw("colz");
            cst->SaveAs(name.c_str());
        }
    }
    integrateCenterTH1F(nSamples, MetPhiMinusPFAK4JetPhi, cst);





}




void plotTH2FHistograms(TTree *sampleTrees[], const int nSamples, double ***doubleVariables, const int nDoubleVariables,std::vector < std::vector < std::vector< TH2F* > > > TH2F_Vector, int  sampleLength[], const int maxEvents, const int doubleArrayIndex[], int ***intVariables, const int nIntVariables, const int intArrayIndex[], const double minDoubleCutValue[], const double maxDoubleCutValue[], const int minIntCutValue[], const int maxIntCutValue[], std::vector< std::vector< TH1F* > > TH1F_cutRatios, const std::string doubleVariableNames[], const std::string intVariableNames[], std::vector < std::vector< cuts* > > cutFlow){

    int numberOfIterations = 0;
    for(int i = 0; i < nSamples; i++){
        numberOfIterations = max(numberOfIterations,sampleLength[i]);
    }
    numberOfIterations = min(numberOfIterations, maxEvents);

    for(int i = 0; i < nSamples; i++){
        for(int j = 0; j < min(numberOfIterations, sampleLength[i]); j++){
        sampleTrees[i]->GetEntry(j);
            if(precut(i, nDoubleVariables, doubleVariables, doubleArrayIndex, intVariables, nIntVariables, intArrayIndex, minDoubleCutValue,  maxDoubleCutValue, minIntCutValue,  maxIntCutValue, TH1F_cutRatios, doubleVariableNames, intVariableNames, cutFlow)){
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

void drawTH1FHistograms(std::vector < std::vector< TH1F* > > TH1F_Vector, const std::string sampleNames[], const int nSamples, const std::string doubleVariableNames[], const int nDoubleVariables,const std::string intVariableNames[], const int nIntVariables, std::vector < std::vector< TH1F* > > TH1F_cutRatios){
            std::string fileName;
            TCanvas *cst = new TCanvas("cst", "TH1F histograms", 10, 10, 700, 700);
            cst->SetLogy();
            TH1F *SandB;
            TH1F *Signal;
            Double_t Norm[] = {.000007207, 61.5267, .83176};

            for(int i = 0; i < nSamples; i++){
                TH1F_cutRatios[i][0]->GetXaxis()->SetBinLabel(1,"All events" ); 
                for(int j = 0; j < nIntVariables; j++){
                    TH1F_cutRatios[i][0]->GetXaxis()->SetBinLabel(j + 2, intVariableNames[j].c_str()); 
                 //   std::cout << intVariableNames[j] << ": " << TH1F_cutRatios[i][0]->GetBinContent(j + 2) << std::endl; 
                }
                for(int j = nIntVariables; j < nIntVariables + nDoubleVariables + 1; j++){
                    TH1F_cutRatios[i][0]->GetXaxis()->SetBinLabel(j + 2, doubleVariableNames[j - nIntVariables].c_str()); 
                 //   std::cout << doubleVariableNames[j - nIntVariables] <<  ": " << TH1F_cutRatios[i][0]->GetBinContent(j + 2) << std::endl; 
                }
                TH1F_cutRatios[i][0]->GetXaxis()->SetBinLabel(nIntVariables + nDoubleVariables + 2, "No other isolated Muons"); 
                //std::cout << "No other Isolate Muons" << ":" << TH1F_cutRatios[i][0]->GetBinContent(nIntVariables + nDoubleVariables + 2) << std::endl; 
                fileName = "TH1F_" + sampleNames[i] +"_cutValues.png";
                TH1F_cutRatios[i][0]->Draw();
                cst->SaveAs(fileName.c_str());
                if(i > 0){
                    fileName = "TH1F_" + sampleNames[i] +"_cutValues_SoverB.png";
                    Signal = (TH1F*)TH1F_cutRatios[0][0]->Clone();       
                    Signal->SetTitle(fileName.c_str());
                    SandB = (TH1F*)TH1F_cutRatios[i][0]->Clone();       
                    SandB->SetTitle(fileName.c_str());
                    //SandB->Add(TH1F_cutRatios[0][0]);
////                    integrateTH1F(Signal);
////                    integrateTH1F(SandB);
//                    sqrtTH1F(SandB);
                    Signal->Divide(SandB);
                    Signal->Draw();
                    cst->SaveAs(fileName.c_str());
                }
            }


            for(int j = 0; j < nDoubleVariables; j++){
                for(int i = 0; i < nSamples; i++){
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

            for(int j = 0; j < nDoubleVariables; j++){
                for(int k = 0; k< nDoubleVariables; k++){
                    for(int i = 0; i < nSamples; i++){
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

    const std::string intVariableNames[] =      {"NPFMuon", "NPFAK4Jets"};
    const int minIntValue[] =                   {1,         1}; 
    const int maxIntValue[] =                   {13,        1};
    const int minIntCutValue[] =                {1,         1}; 
    const int maxIntCutValue[] =                {10,        100};
    const int intArraySize[] =                  {1,         1};
    const int intArrayIndex[] =                 {0,         0};

    const std::string doubleVariableNames[] =   {"MetPt",   "PFMuonPt", "MetPhi",   "PFMuonPhi",    "MetSumEt",     "PFMuonEta",    "PFMuonChargedHadronIso",   "PFAK4JetPt",   "PFAK4JetEta",  "PFAK4JetPhi",  "PFAK4JetBTagCSVv2", "PFMuonPhotonIso", "PFMuonNeutralHadronIso" };
    const double minDoubleValue[] =                {10,        10,         -4,         -4,             10,             -3,             -20,                          10,             -3,             -4,           0,                      -20,              -20};
    const double maxDoubleValue[] =                {700,       700,        4,          4,              4000,           3,              200,                         700,            3,              4,            20,                     200,             200};
    const double minDoubleCutValue[] =             {40,        33,         -4,         -4,             10,             -2.4,             -20,                          30,             -2.4,             -4,           -20000,                      -20,              -20};
    const double maxDoubleCutValue[] =             {10000,       10000,        4,          4,              4000,           2.4,              5000,                         700,            2.4,              4,            20,                     5000,             5000};
    const int doubleArraySize[] =               {40,        13,         40,         13,             40,             13,             13,                         16,             16,             16,             16,                     13,             13};
    const int doubleArrayIndex[] =              {10,        0,          10,         0,              10,             0,              0,                          0,              0,              0,              0,                      0,              0};

    const int nBins = 100;
    const int nDoubleVariables = sizeof(doubleVariableNames)/sizeof(doubleVariableNames[0]);
    const int nIntVariables = sizeof(intVariableNames)/sizeof(intVariableNames[0]);
    const int maxEvents = 10000; 

    ///variables
    int sampleLength[nSamples];
    double ***doubleVariables = new double**[nSamples];
    int ***intVariables = new int**[nSamples];
    std::string histName;
    std::string hist2DName;

    //root objects
    TTree *sampleTrees[nSamples];

    //TH1F **TH1F_histograms = new TH1F*[nSamples];
    TH1F *TH1F_histograms[nSamples][nDoubleVariables];
    std::vector < std::vector< TH1F* > > TH1F_Vector;
    TH1F_Vector.resize(nSamples);

    std::vector < std::vector < TH1F*> > TH1F_cutRatios;
    TH1F_cutRatios.resize(nSamples);

    std::vector < std::vector < std::vector < TH2F* > > > TH2F_Vector;
    TH2F_Vector.resize(nSamples);


    std::vector < std::vector < cuts* > > cutFlow;
    
    cutFlow.resize(nSamples);
 //   TH1F **ptrTo_TH1F_histograms = TH1F_histograms;

    //TH1F_histograms[0][0] = new TH1F("test","test",100,10, 100);
    

    loadRootFiles(sampleNames, nSamples, sampleTrees, sampleLength);

    for(int i = 0; i < nSamples; i++){
        
        TH2F_Vector[i].resize(nDoubleVariables);
        histName = "TH1F_" + sampleNames[i] + "_numbersCut";
        //cutFlow[i].resize(nDoubleVariables + nIntVariables + 2);
        TH1F_cutRatios[i].push_back(new TH1F(histName.c_str(),histName.c_str(), nDoubleVariables + nIntVariables + 5, -1, nDoubleVariables + nIntVariables+1));
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
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j <  nDoubleVariables + nIntVariables + 7; j++) {
            cuts *cutFlow_ptr = new cuts;
            cutFlow[i].push_back(cutFlow_ptr);//
          //  std::cout << "what" << std::endl;
            cutFlow[i][j]->cutVariable = "test";
            cutFlow[i][j]->count = 0;
            std::cout << i << " " << cutFlow[i][j]->cutVariable << cutFlow[i][j]->count << std::endl;
        }
    }

    //makeTH1FHistograms(sampleNames, nSamples, doubleVariableNames, nDoubleVariables, nBins, minDoubleValue, maxDoubleValue,  TH1F_Vector);

    getBranches(doubleVariableNames, nDoubleVariables, doubleArraySize, doubleVariables, sampleTrees, nSamples);    

    getIntBranches(intVariableNames, nIntVariables, intArraySize, intVariables, sampleTrees, nSamples);    

    plotTH1FHistograms(sampleTrees, nSamples, doubleVariables,  nDoubleVariables, TH1F_Vector, sampleLength, maxEvents, doubleArrayIndex, intVariables, nIntVariables, intArrayIndex, minDoubleCutValue, maxDoubleCutValue, minIntCutValue, maxIntCutValue, TH1F_cutRatios, doubleVariableNames, intVariableNames, cutFlow);

    plotTH2FHistograms(sampleTrees, nSamples, doubleVariables,  nDoubleVariables, TH2F_Vector, sampleLength, maxEvents, doubleArrayIndex, intVariables, nIntVariables, intArrayIndex, minDoubleCutValue, maxDoubleCutValue, minIntCutValue, maxIntCutValue, TH1F_cutRatios, doubleVariableNames, intVariableNames, cutFlow);




    drawTH1FHistograms(TH1F_Vector, sampleNames, nSamples, doubleVariableNames, nDoubleVariables, intVariableNames, nIntVariables, TH1F_cutRatios);
//
    drawTH2FHistograms(TH2F_Vector, sampleNames, nSamples, doubleVariableNames, nDoubleVariables);

            Double_t Norm[] = {.000007207, 61.5267, .83176};
    for(std::vector<std::vector < cuts* > >::size_type i = 0; i != cutFlow.size(); i++) {
        for(std::vector < cuts* >::size_type j = 0; j != cutFlow[i].size(); j++) {
            std::cout << i << '\t' << cutFlow[i][j]->cutVariable<< '\t' << (cutFlow[i][j]->count) << std::endl;
        }
    }

//    for(k = 0; k < nDoubleVariables + nIntVariables + 2; k++){
//        std::cout << i << " " << k << std::endl;
//        cutFlow[i][k].count = 0;
//    }

}
