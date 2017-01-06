

#define nMasses 20
#define lumi 20000 #luminosity pb-1
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
  TGraph *nEvtsObs = new TGraph(0);
  TGraph *nEvtsExp = new TGraph(0);
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
    if (nentries<6){
      continue;
    }
    tree[counter]->SetBranchAddress("limit",&limit);
    tree[counter]->SetBranchAddress("mh",&mh);
    for(int i=0; i<nentries; ++i){
      tree[counter]->GetEntry(i);
      //std::cout<<"limit: "<<limit<<std::endl;
      //std::cout<<"mass: "<<mh<<std::endl;
      if(i==0) twoSigmaLow=limit;
      if(i==1) oneSigmaLow=limit;
      if(i==2) expected=limit;
      if(i==3) oneSigmaHigh=limit;
      if(i==4) twoSigmaHigh=limit;
      if(i==5) observed=limit;
      mass=mh;

    }
    std::cout<<"Mass: "<<mh<<std::endl;
    std::cout<<"obs: "<<observed<<" 2sig: "<<twoSigmaLow<<" 1sig: "<<oneSigmaLow<<" exp: "<<expected<<" 1sig: "<<oneSigmaHigh<<" 2sig: "<<twoSigmaHigh<<std::endl;
    outputList<<mh<<" "<<tau21.c_str()<<" "<<observed<<" "<<twoSigmaLow<<" "<<oneSigmaLow<<" "<<expected<<" "<<oneSigmaHigh<<" "<<twoSigmaHigh<<std::endl;
    
    limitMassExp->SetPoint(counter, mass, expected);
    limitMassObs->SetPoint(counter, mass, observed);
    limitMass1Sigma->SetPoint(counter, mass, expected);
    limitMass1Sigma->SetPointError(counter,0,0, expected-oneSigmaLow, oneSigmaHigh-expected);
    limitMass2Sigma->SetPoint(counter, mass, expected);
    limitMass2Sigma->SetPointError(counter,0,0, expected-twoSigmaLow, twoSigmaHigh-expected);

    nEvtsExp->SetPoint(counter, mass, expected*lumi);
    nEvtsObs->SetPoint(counter, mass, observed*lumi);

    ++counter;
  }


  TMultiGraph *limitGraph = new TMultiGraph();
  limitGraph->Add(limitMass2Sigma);
  limitGraph->Add(limitMass1Sigma);

  limitMass2Sigma->SetFillColor(5);
  limitMass1Sigma->SetFillColor(kGreen-7);
  limitMassExp->SetLineStyle(2);
  limitMass1Sigma->SetLineStyle(2);
  limitMass2Sigma->SetLineStyle(2);
  limitMass1Sigma->SetLineWidth(2);
  limitMass2Sigma->SetLineWidth(2);
  limitMassObs->SetLineWidth(2);
  TCanvas *c = new TCanvas("c","Limits Mass Spectrum", 1);
  limitMassExp->GetXaxis()->SetTitle("resonance mass [GeV]");
  limitMass1Sigma->GetXaxis()->SetTitle("resonance mass [GeV]");
  limitMass2Sigma->GetXaxis()->SetTitle("resonance mass [GeV]");
  std::cout<<"Creating canvas"<<std::endl;
  
  limitMassExp->GetYaxis()->SetRangeUser(.1,10000);
  c->SetLogy();
  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry((TObject*)0, "95\% C.L.", "");
  leg->AddEntry(limitMassExp, "Exp.", "l");
  leg->AddEntry(limitMass1Sigma, "1#sigma", "fl");
  leg->AddEntry(limitMass2Sigma, "2#sigma", "fl");

  //limitMassExp->GetXaxis()->SetTitle("mass [GeV]");
  limitMassExp->GetYaxis()->SetTitle("limit");
  //limitMass2Sigma->Draw("ZAP");
  //limitMassObs->Draw("PSAME");
  //limitMass2Sigma->Draw("a3SAME");
  //limitMass1Sigma->Draw("a3SAME");
//  limitMassExp->GetYaxis()->SetRangeUser(-1.,.01);
  limitMassExp->Draw("AP");
  limitGraph->Draw("a3");
  limitMassExp->Draw("LSAME");
  limitMassObs->Draw("CSAME");
  leg->Draw("SAME");
  limitGraph->GetXaxis()->SetTitle("resonance mass [GeV]");
  limitGraph->GetYaxis()->SetTitle("#sigma[pb]");

  TLatex latex;
  TString lumiText = "36 fb^{-1} (13 TeV)";//"#it{L}=12.9 fb^{-1} (2016) (13 TeV)";//"#bf{CMS Preliminary} #it{L}=12.9 fb^{-1}";
  TString cmsLogo = "#bf{CMS}  #it{Preliminary}";
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    
  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  float lumiTextOffset = 0.2;
  float lumiTextSize = 0.035;
  latex.SetTextSize(lumiTextSize);    
  latex.DrawLatex(.84,.96,lumiText);//1-r,1-t+lumiTextOffset*t,lumiText);
  latex.DrawLatex(.44/*.25*/,.96,cmsLogo);
  latex.Draw();

  TGaxis::SetMaxDigits(2);

  gPad->Modified();
  gPad->Update();
  c->Modified();
  c->Update();
  std::cout<<"Saving output"<<std::endl;
  replace( tau21.begin(), tau21.end(), '.', 'p' );
  c->SaveAs(Form("limits_%s.pdf",tau21.c_str()));



  TCanvas *c1 = new TCanvas("c1","# events Spectrum", 1);
  nEvtsExp->GetXaxis()->SetTitle("resonance mass [GeV]");
  std::cout<<"Creating canvas for nEvents"<<std::endl;
  
  nEvtsExp->GetYaxis()->SetRangeUser(.1,10000);
  c1->SetLogy();
  TLegend *leg1 = new TLegend(0.7,0.7,0.9,0.9);
  leg1->AddEntry(nEvtsExp, "Exp.", "l");
  leg1->AddEntry(nEvtsObs, "Obs.", "l");

  nEvtsExp->GetYaxis()->SetTitle("#Evts");
  nEvtsExp->Draw("AP");
  nEvtsObs->Draw("CSAME");
  leg1->Draw("SAME");
  latex.Draw();
  TGaxis::SetMaxDigits(2);
  gPad->Modified();
  gPad->Update();
  c1->Modified();
  c1->Update();
  std::cout<<"Saving output"<<std::endl;
  replace( tau21.begin(), tau21.end(), '.', 'p' );
  c1->SaveAs(Form("nEvents_%s.pdf",tau21.c_str()));

  return;

}
