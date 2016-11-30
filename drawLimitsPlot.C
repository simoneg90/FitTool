

#define nMasses 20
void drawLimitsPlot(std::string tau21="020"){

  std::cout<<"Entering in drawLimitsPlot program"<<std::endl;

  std::ifstream inputList;
  inputList.open("/cmshome/gellisim/CMSSW_7_4_14/src/FitTool/limitFileList.txt");
  if(inputList==NULL){
        std::cout<<"ERROR! File: limitFileList.txt doesn't exist!"<<std::endl;
        exit(-1);
  }

  TTree *tree[nMasses];
  TFile *inputFiles[nMasses];
  int counter=0;
  int nentries;
  double limit, mh;
  double expected, observed, oneSigmaLow, oneSigmaHigh, twoSigmaLow, twoSigmaHigh;
  int mass;
  TGraph *limitMassObs = new TGraph(0);
  TGraph *limitMassExp = new TGraph(0);
  TGraphAsymmErrors *limitMass1Sigma = new TGraphAsymmErrors(0);
  TGraphAsymmErrors *limitMass2Sigma = new TGraphAsymmErrors(0);

  std::string fileName="";

  std::ofstream outputList;
  outputList.open("/cmshome/gellisim/CMSSW_7_4_14/src/FitTool/limitOutputList.txt", std::fstream::app);
  while(true){

    if(inputList.eof()) break;
    
    if(counter>=nMasses){
      std::cout<<"ATTENTION, too many files added in the list"<<std::endl;
      exit(-1);
    }
    
    inputList>>fileName;
    std::cout<<"Trying to open: "<<fileName.c_str()<<std::endl;
    inputFiles[counter] = TFile::Open(fileName.c_str());
    if(inputFiles[counter]==NULL){
      std::cout<<"ERROR! File: "<<fileName.c_str()<<std::endl;
      exit(-1);
    }
    tree[counter] = (TTree *)inputFiles[counter]->Get("limit");
    nentries = tree[counter]->GetEntries();
    std::cout<<"Entries: "<<nentries<<std::endl;
    tree[counter]->SetBranchAddress("limit",&limit);
    tree[counter]->SetBranchAddress("mh",&mh);
    for(int i=0; i<nentries; ++i){
      tree[counter]->GetEntry(i);
      //std::cout<<"limit: "<<limit<<std::endl;
      //std::cout<<"mass: "<<mh<<std::endl;
      if(i==0) oneSigmaLow=limit;
      if(i==1) twoSigmaLow=limit;
      if(i==2) expected=limit;
      if(i==3) oneSigmaHigh=limit;
      if(i==4) twoSigmaHigh=limit;
      if(i==5) observed=limit;
      mass=mh;

    }
    std::cout<<"Mass: "<<mh<<std::endl;
    std::cout<<observed<<" "<<twoSigmaLow<<" "<<oneSigmaLow<<" "<<expected<<" "<<oneSigmaHigh<<" "<<twoSigmaHigh<<std::endl;
    outputList<<mh<<" "<<tau21.c_str()<<" "<<observed<<" "<<twoSigmaLow<<" "<<oneSigmaLow<<" "<<expected<<" "<<oneSigmaHigh<<" "<<twoSigmaHigh<<std::endl;
    
    limitMassExp->SetPoint(counter, mass, expected);
    limitMassObs->SetPoint(counter, mass, observed);
    limitMass1Sigma->SetPoint(counter, mass, expected);
    limitMass1Sigma->SetPointError(counter,0,0, oneSigmaLow, oneSigmaHigh);
    limitMass2Sigma->SetPoint(counter, mass, expected);
    limitMass2Sigma->SetPointError(counter,0,0, twoSigmaLow, twoSigmaHigh);

    ++counter;
  }


  TMultiGraph *limitGraph = new TMultiGraph();
  limitGraph->Add(limitMass2Sigma);
  limitGraph->Add(limitMass1Sigma);

  limitMass2Sigma->SetFillColor(5);
  limitMass1Sigma->SetFillColor(kGreen-7);
  limitMassExp->SetLineStyle(2);
  limitMassObs->SetLineWidth(2);
  std::cout<<"Creating canvas"<<std::endl;
  TCanvas *c = new TCanvas("c","Limits Mass Spectrum", 1);
  limitMassExp->GetYaxis()->SetRangeUser(.1,10000);
  c->SetLogy();

  //limitGraph->GetYaxis()->SetRangeUser(.1,10000);
  limitMassExp->GetXaxis()->SetTitle("mass [GeV]");
  limitMassExp->GetYaxis()->SetTitle("limit");
  //limitMass2Sigma->Draw("ZAP");
  //limitMassObs->Draw("PSAME");
  //limitMass2Sigma->Draw("a3SAME");
  //limitMass1Sigma->Draw("a3SAME");
  limitMassExp->Draw("ZAP");
  limitGraph->Draw("a3");
  limitMassExp->Draw("LSAME");
  limitMassObs->Draw("CSAME");
  std::cout<<"Saving output"<<std::endl;
  replace( tau21.begin(), tau21.end(), '.', 'p' );
  c->SaveAs(Form("limits_%s.pdf",tau21.c_str()));



  return;

}
