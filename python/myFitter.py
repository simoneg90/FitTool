from optparse import OptionParser
import time
from ROOT import *
import rootTools
from framework import Config
from array import *
from itertools import *
from operator import *
import os
import random
import sys
import math

#gSystem.Load("/cmshome/gellisim/CMSSW_7_4_14/src/FitTool/src/RooDCBShape_cxx.so")
def drawBiasPlot(biasRootFile, outDir, funcName1, funcName2, suff=""):
  fileBias = TFile.Open(biasRootFile)
  histoBias = TH1F("histoBias", "histoBias", 20,-4,4)
  fileBias.tree_fit_sb.Project(histoBias.GetName(),"(mu-1)/muErr")
  gausFunction = TF1("gausFunction", "gaus(0)", -4.,4.)
  print "Setting gaus parameters: mean ", histoBias.GetMean(), " sigma ", histoBias.GetRMS(), "norm ", histoBias.Integral() 
  gausFunction.SetParameters(histoBias.Integral(), histoBias.GetMean(), histoBias.GetRMS())
  histoBias.Fit(gausFunction, "R+")
  counter = 0
  thisNorm = -999.
  thisMean = -999.
  thisSigma = -999.
  maxIterations = 3

  while ( (abs((gausFunction.GetParameter(1)- thisMean)/(gausFunction.GetParameter(1)+1.))>0.001) and (abs((gausFunction.GetParameter(2) - thisSigma ) / (gausFunction.GetParameter(2)+1) ) > 0.001) and (counter < maxIterations) ): #+1 at denominator prevent crashes when runnin on the asimov with no signal
    thisNorm = gausFunction.GetParameter(0)
    thisMean = gausFunction.GetParameter(1)
    thisSigma = gausFunction.GetParameter(2)
    gausFunction.SetParameters(thisNorm, thisMean, thisSigma)
    gausFunction.SetRange(thisMean - 2.*thisSigma, thisMean + 2.*thisSigma)
    histoBias.Fit(gausFunction, "R+")
    counter += 1

  print "Gauss fit exited after %i iterations"%(counter)
  c_bias = TCanvas("c_bias", "c_bias", 1)
  histoBias.Draw()
  c_bias.SaveAs("%s/bias_histogram_%s_%s%s.pdf"%(outDir, funcName1, funcName2, suff))
  histoBias.Delete()

  return

def prepareBiasDatacard (workspace, fit_result_fileName, outDir, mass, dataset, realVar, rangeMin, rangeMax):
  
  
  nFunctions = 4  
  dijet_p1 = RooRealVar("dijet_p1", "dijet_p1", -13., -15., -9.)
  dijet_p2 = RooRealVar("dijet_p2", "dijet_p2", -1.4, -2., -1.)
  dijet_p3 = RooRealVar("dijet_p3", "dijet_p3", .2,0.01,10.)
  dijet_function = RooGenericPdf("dijet_function","pow((1-(@0/13000)),@3)*pow(@0/13000,(@1+@2*log(@0/13000)))", RooArgList(realVar,dijet_p1,dijet_p2, dijet_p3))

  expLaw1_p1 = RooRealVar("expLaw1_p1", "expLaw1_p1", -13., -15., -9.)
  expLaw1_p2 = RooRealVar("expLaw1_p2", "expLaw1_p2", -1.4, -2., -1.)
  expLaw1_function = RooGenericPdf("expLaw1_function","pow((@0/13000),@2)*exp(@1*(@0/13000))", RooArgList(realVar,expLaw1_p1,expLaw1_p2))

  expLaw2_p1 = RooRealVar("expLaw2_p1", "expLaw2_p1", -13., -15., -9.)
  expLaw2_p2 = RooRealVar("expLaw2_p2", "expLaw2_p2", -1.4, -2., -1.)
  expLaw2_p3 = RooRealVar("expLaw2_p3", "expLaw2_p3", .2,0.01,10.)
  expLaw2_p4 = RooRealVar("expLaw2_p4", "expLaw2_p4", .2,0.01,10.)
  expLaw2_function = RooGenericPdf("expLaw2_function","pow((@0/13000),-(@2*@2)+@4)*exp(@1*(@0/13000)+@3*(@0/13000)**2)", RooArgList(realVar,expLaw2_p1,expLaw2_p2, expLaw2_p3,expLaw2_p4))
 
  #params.insert(pair<string,RooRealVar*>(logc, new RooRealVar(logc.c_str(),logc.c_str(),60.,-200.0,500.)));
  #}else if(i==2 && order==3){
  #params.insert(pair<string,RooRealVar*>(logc, new RooRealVar(logc.c_str(),logc.c_str(),30.,-200.0,500.)));
  #}else if(i==3 && order==3){
  #params.insert(pair<string,RooRealVar*>(logc, new RooRealVar(logc.c_str(),logc.c_str(),3.,-200.0,500.)));
  #}
  atlas1_p1 = RooRealVar("atlas1_p1", "atlas1_p1", 10., -200., 500.)
  atlas1_p2 = RooRealVar("atlas1_p2", "atlas1_p2", -10., -200., 500.)
  atlas1_p3 = RooRealVar("atlas1_p3", "atlas1_p3", -2., -200., 500.)
  atlas1_function = RooGenericPdf("atlas1_function","@1*(pow((1-pow((@0/13000),1/3)),@2)*pow(@0/13000,@3))", RooArgList(realVar,atlas1_p1,atlas1_p2, atlas1_p3))

  atlas2_p1 = RooRealVar("atlas2_p1", "atlas2_p1", 10., -200., 500.) 
  atlas2_p2 = RooRealVar("atlas2_p2", "atlas2_p2", -10., -200., 500.)
  atlas2_p3 = RooRealVar("atlas2_p3", "atlas2_p3", -2., -200., 500.)
  atlas2_p4 = RooRealVar("atlas2_p4", "atlas2_p3", 60., -200., 500.)
  atlas2_function = RooGenericPdf("atlas2_function","@1*(pow((1-pow((@0/13000),1/3)),@2)*pow(@0/13000,@3+(@4*(log(@0/13000))**2)))", RooArgList(realVar,atlas2_p1,atlas2_p2, atlas2_p3, atlas2_p4))

  
  dijet_function.fitTo(dataset, RooFit.Range(rangeMin, rangeMax))
  expLaw1_function.fitTo(dataset, RooFit.Range(rangeMin, rangeMax))
  expLaw2_function.fitTo(dataset, RooFit.Range(rangeMin, rangeMax))
  atlas1_function.fitTo(dataset, RooFit.Range(rangeMin, rangeMax))
  #atlas2_function.fitTo(dataset, RooFit.Range(rangeMin, rangeMax))

  frame = realVar.frame()
  dataset.plotOn(frame)
  dijet_function.plotOn(frame, RooFit.LineColor(kRed), RooFit.Name("dijet"))
  expLaw1_function.plotOn(frame, RooFit.LineColor(kGreen), RooFit.Name("exp1"))
  expLaw2_function.plotOn(frame, RooFit.LineColor(kBlue), RooFit.Name("exp2"))#, RooFit.LineStyle(2))
  atlas1_function.plotOn(frame, RooFit.LineColor(kViolet), RooFit.Name("atlas"), RooFit.LineStyle(2), RooFit.LineWidth(4))
  #atlas2_function.plotOn(frame, RooFit.LineColor(kBlack))#, RooFit.LineStyle(2), RooFit.LineWidth(3))

  leg = TLegend(0.7,0.7,0.9,0.9)
  leg.AddEntry(frame.findObject("dijet"), "dijet", "l")
  leg.AddEntry(frame.findObject("exp1"), "exp1", "l")
  leg.AddEntry(frame.findObject("exp2"), "exp2", "l")
  leg.AddEntry(frame.findObject("atlas"), "atlas", "l")
  c_bias = TCanvas("c_bias", "c_bias", 1)
  c_bias.SetLogy()
  frame.Draw()
  leg.Draw("SAME")
  c_bias.SaveAs("%s/bias_fitPlot_m%s.pdf"%(outDir, str(mass)))
  
  cat = RooCategory ("pdf_index","Index of Pdf which is active")
  multiPdf = RooMultiPdf ("multiPdf","All Pdfs for bias study",cat,RooArgList(dijet_function,expLaw1_function,expLaw2_function,atlas1_function))
  norm = RooRealVar ("multiPdf_norm","Number of background events",0,20000)

  dijet_p1.setConstant(True) 
  dijet_p2.setConstant(True) 
  dijet_p3.setConstant(True) 
  expLaw1_p1.setConstant(True)
  expLaw1_p2.setConstant(True)
  expLaw2_p1.setConstant(True)
  expLaw2_p2.setConstant(True)
  expLaw2_p3.setConstant(True)
  expLaw2_p4.setConstant(True)
  atlas1_p1.setConstant(True)
  atlas1_p2.setConstant(True)
  atlas1_p3.setConstant(True)
  norm.setConstant(True)


  WS_bias = RooWorkspace("WS_bias")#, kTRUE)
  rootTools.Utils.importToWS(WS_bias, cat)
  rootTools.Utils.importToWS(WS_bias, norm)
  rootTools.Utils.importToWS(WS_bias, multiPdf)

  WS_bias.Print("v")
  root_biasName = "%s/WS_bias_%s.root"%(outDirName,str(options.res_mass))
#  bias_fileName = "%s/datacard_bias_%s.root"%(outDirName,str(options.res_mass))
  WS_bias.SaveAs("%s/WS_bias_%s.root"%(outDirName,str(options.res_mass)))

  bias_dataCardName = writeCard(workspace, mass, fit_result_fileName, outDirName, False, True, WS_bias.GetName(), root_biasName, suffix="_bias")


  return bias_dataCardName, nFunctions, root_biasName

def doBias(workspace, fit_result_fileName, outDir, datacardFileName, mass, rangeVarLow, rangeVarHigh, mVg, dataSet, nToys=100, expectSignal=1.):

  #note: to do a 0 nExpected bias study, use nExpected=0.000001
  bias_datacard, nFs, bias_rootFile = prepareBiasDatacard(workspace, fit_result_fileName, outDir, mass, dataSet, mVg, rangeVarLow, rangeVarHigh)
  myPrint("Doing Bias test", "*")
  print "Datacard bias name: ", bias_datacard
  function_list = ['dijet', 'expLaw1', 'expLaw2', 'atlas']
  
  for i in range(0, nFs):

    os.system("combine %s -M GenerateOnly --setPhysicsModelParameters pdf_index=%s --toysFrequentist -t %s --expectSignal %s --saveToys -m %s -n %s --freezeNuisances pdf_index"%(bias_datacard, str(i), str(nToys), str(expectSignal), str(mass), "bias_"+function_list[i]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys)))
    os.system("mv higgsCombine%s.GenerateOnly.mH%s.123456.root %s/higgsCombine%s.GenerateOnly.mH%s.123456.root"%("bias_"+function_list[i]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys), str(mass), outDir, "bias_"+function_list[i]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys), str(mass)))
    for j in range(0, nFs):
      os.system("combine %s -M MaxLikelihoodFit  --setPhysicsModelParameters pdf_index=%s --toysFile %s/higgsCombine%s.GenerateOnly.mH%s.123456.root  -t %s --rMin -10 --rMax 10 -n %s -m %s --expectSignal %s --freezeNuisances pdf_index --plots"%(bias_datacard, str(j), outDir, "bias_"+function_list[i]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys), str(mass), str(nToys), "toyF_"+function_list[i]+"_fitF_"+function_list[j]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys), str(mass), str(expectSignal)))
      os.system("mv higgsCombine%s.MaxLikelihoodFit.mH%s.123456.root %s/higgsCombine%s.GenerateOnly.mH%s.123456.root"%("toyF_"+function_list[i]+"_fitF_"+function_list[j]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys), str(mass), outDir, "toyF_"+function_list[i]+"_fitF_"+function_list[j]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys), str(mass)))
      os.system("mv mlfit%s.root %s/mlfit%s.root"%("toyF_"+function_list[i]+"_fitF_"+function_list[j]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys), outDir, "toyF_"+function_list[i]+"_fitF_"+function_list[j]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys)))
      drawBiasPlot("%s/mlfit%s.root"%(outDir, "toyF_"+function_list[i]+"_fitF_"+function_list[j]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys)), outDir, function_list[i], function_list[j], "_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys))

      os.system("mv ./WS_Vg_prefit.png %s/WS_Vg_prefit_%s.png"%(outDir, "toyF_"+function_list[i]+"_fitF_"+function_list[j]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys)))
      os.system("mv ./covariance_fit_b.png %s/covariance_fit_b_%s.png"%(outDir, "toyF_"+function_list[i]+"_fitF_"+function_list[j]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys)))
      os.system("mv ./WS_Vg_fit_b.png %s/WS_Vg_fit_b_%s.png"%(outDir, "toyF_"+function_list[i]+"_fitF_"+function_list[j]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys)))
      os.system("mv ./WS_Vg_fit_s.png %s/WS_Vg_fit_s_%s.png"%(outDir, "toyF_"+function_list[i]+"_fitF_"+function_list[j]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys)))
      os.system("mv ./covariance_fit_s.png %s/covariance_fit_s_%s.png"%(outDir, "toyF_"+function_list[i]+"_fitF_"+function_list[j]+"_expected"+str(expectSignal).replace('.','p')+"_nToys"+str(nToys)))
#    os.system("combine -M GenerateOnly %s -n %s --toysFrequentist --saveToys --expectSignal %s -t %s -m %s"%(datacardFileName, "bias", str(expectSignal), str(nToys), str(mass)))
#    os.system("combine -M MaxLikelihoodFit %s -n %s --toysFile higgsCombine%s.GenerateOnly.mH120.123456.root -t %s -m %s --minimizerTolerance 0.01 --minimizerStrategy 2 --saveWorkspace"%(datacardFileName, "biasMaxLikelihood", "bias", str(nToys), str(mass)))


  return


def doManualLimitScan(nSteps, rMin, rMax, mass, outDir, datacardName):

  step = (rMax - rMin)/nSteps
  #cont = 0.1
  i=0
  #for i in range(0, nSteps+nSteps/10):
  myPrint("Doing limit scan")
  scanPoint = 0
  median = rMin + (rMax-rMin)/4
  sigma = (rMax-rMin)/3
  print "Range %f - %f, median %f, sigma %f" %(rMin, rMax, median, sigma)
  #while (scanPoint<(rMax+rMax/10)): #scan above rMax of 10%


  while (i<nSteps):

    scanPoint = abs(random.lognormvariate(median, sigma)-1)
    print "Scan step: ", i, " r: ", scanPoint
    combineFileName = "higgsCombine%s.Asymptotic.mH%s.root"%("_scan_"+str(scanPoint).replace('.','p'), str(mass))
    os.system("combine -M Asymptotic -m %s --singlePoint %s -n %s %s"%(str(mass), str(scanPoint), "_scan_"+str(scanPoint).replace('.','p'),datacardName))
    print "Moving %s in the right directory"%(combineFileName)
    os.system("mv %s %s/%s"%(combineFileName, outDir, combineFileName))
    i += 1
                
  os.system("hadd limitScan.root %s/higgsCombine_scan_*.Asymptotic.mH%s.root"%(outDir, str(mass)))

  nScansFileSuff = ""
  if nScansFile>0:
    #nScansFile+=1
    nScansFileSuff = "_"+str(nScansFile)
  print "Suffix: ", nScansFileSuff
  nScansFile = int(os.popen('ls %s/limitScan* | wc -l'%(outDir)).read())
  print "Renaming the old file: %s/limitScan_%s.root"%(outDir, str(nScansFile))
  os.system("mv limitScan.root %s/limitScan%s.root"%(outDir, nScansFileSuff))
  os.system("rm %s/higgsCombine_scan_*.Asymptotic.mH%s.root"%(outDir, str(mass)))

  return



def getGoodRange(histo, realVar=None, isRooDataHist=False):

  histogram = TH1F()
  histogram.SetName("histogram")
  print isRooDataHist
  if isRooDataHist:
    print "Extracting histogram from rooDataHist"
    histogram = histo.createHistogram("histogram", realVar)
  else:
    histogram = histo
  binStep = (histogram.GetXaxis().GetXmax()-histogram.GetXaxis().GetXmin())/histogram.GetNbinsX()
  counter_zeros=0
  SR_lo = -1.
  SR_hi = -1.
  for b in range(0,histogram.GetNbinsX()):
    if(histogram.GetBinContent(b)>0 and SR_lo<0):
      low_tmp =((histogram.GetXaxis().GetXmax()-histogram.GetXaxis().GetXmin())/histogram.GetNbinsX()*b) + histogram.GetXaxis().GetXmin()
      if(low_tmp>600):
        SR_lo=low_tmp
        counter_zeros=0
    if(histogram.GetBinContent(b)>0):
      val_tmp=b
      counter_zeros=0
    if(histogram.GetBinContent(b)==0):
      counter_zeros+=1
    if(counter_zeros==4 and SR_lo>0):
      break

  SR_hi = ((histogram.GetXaxis().GetXmax()-histogram.GetXaxis().GetXmin())/histogram.GetNbinsX()*val_tmp) + histogram.GetXaxis().GetXmin()
  print "Fit range: %f - %f"%(SR_lo, SR_hi)
  histogram.Delete()
  return SR_lo, SR_hi

def writeRootDatacard(fileName, mass):
  myPrint("Launching text2workspace command", "*")
  os.system("text2workspace.py %s -o %s -m %f" % (fileName, fileName.replace('.txt','.root'), mass))
  return

def searchCombineRLimits(combineFileName, iteration):
  myPrint("Looking for more constrained range in r for combine","*")
  print "Combine file: ", combineFileName
  combineFile = TFile.Open(combineFileName, "r")
#  os.system("lsof %s"%(combineFileName) )
  if not combineFile:
    sys.exit("Warning, failed to open %s"%(combineFileName))
#  combineTree = combineFile.Get("limit")
#  combineTree.SetBranchAddress("limit", limit)
  #combineTree.GetEntry(0)
  for entry in combineFile.limit:#combineTree:
    limitObs = combineFile.limit.limit
    print combineFile.limit.limit
    #break #to read only the first entry
  print "Observed limit: ", limitObs #to take the obsLimit. Last entry
  N = combineFile.limit.GetEntries()
  rmin=0
  rmax=limitObs
  if iteration>=1:
    rmin = limitObs/10
  elif iteration>=2:
    rmax = limitObs*10
  elif iteration>=3:
    rmin = limitObs/1e2
  elif iteration>=4:
    rmax = limitObs*1e2
  elif iteration>=5:
    rmin = limitObs/1e3
  elif iteration>=6:
    rmax = limitObs*1e3
  elif iteration>=7:
    rmin = limitObs/1e4
  elif iteration>=8:
    rmax = limitObs*1e4
  elif iteration>=9:
    rmin = limitObs/1e5
  elif iteration>=10:
    rmax = limitObs*1e5

  combineFile.Close()
  print "Range for r: %f - %f"%(rmin, rmax)
  return N, rmin, rmax#*100.

def runCombine(mass, fileName, outDir, doAsimov=False, isBlinded=False):
  blind = ""
  nMax = 6
  if doAsimov:
    myPrint("Warning: asimov option activated! Combine running on one toy!", "-")
    asimov=-1
    combineFileName = "higgsCombineTest.Asymptotic.mH%s.123456.root"%(str(mass))
  else:
    asimov=0
    combineFileName = "higgsCombineTest.Asymptotic.mH%s.root"%(str(mass))
  if isBlinded:
    myPrint("Blind analysis: no run on observed","-")
    blind = "--run expected"
    nMax = 5
  myPrint("Running combine command", "*")
  print "file for combine: ", fileName
#  os.system("combine -M Asymptotic %s"%(fileName)) #--verbose 9"%(str(mass), fileName)) #-m mass not used as for high resonances does not seem to work
#  os.system("combine -M Asymptotic -m %s -t 0 --rMin=0 --rMax=0.004 %s"%(str(mass), fileName))
  myPrint("Throwing asimov with MaxLikelihood fit")
  os.system("combine -M MaxLikelihoodFit -m %s -t -1 --expectSignal=1 %s --plots"%(str(mass), fileName)) #throw and use the Asimov to see the signal in the bkg distribution and test the fit
  print "Moving asimov output in right directory"
  os.system("mv ./WS_Vg_prefit.png %s/WS_Vg_prefit.png"%(outDir))
  os.system("mv ./covariance_fit_b.png %s/covariance_fit_b.png"%(outDir))
  os.system("mv ./WS_Vg_fit_b.png %s/WS_Vg_fit_b.png"%(outDir))
  os.system("mv ./WS_Vg_fit_s.png %s/WS_Vg_fit_s.png"%(outDir))
  os.system("mv ./covariance_fit_s.png %s/covariance_fit_s.png"%(outDir))
  os.system("mv ./mlfit.root %s/mlfit.root"%(outDir))
  os.system("mv ./higgsCombineTest.MaxLikelihoodFit.mH%s.123456.root %s/higgsCombineTest.MaxLikelihoodFit.mH%s.123456.root"%(str(mass), outDir, str(mass)))
  myPrint("First combine attempt")
  os.system("combine -M Asymptotic -m %s -t %s %s %s"%(str(mass), str(asimov), fileName, blind))
#  os.system("combine -M Asymptotic -m %s -H ProfileLikelihood %s --run expected"%(str(mass), fileName))
#  os.system("mv higgsCombineTest.Asymptotic.mH120.root %s/higgsCombineTest.Asymptotic.mH%s.root" % (outDir, str(mass)))

 
  for i in range(0, 11):
    N, rMin, rMax = searchCombineRLimits(combineFileName, i)
    if N==nMax:
      print "Correctly extracted the limit"
      break
#    os.system("rm ./%s"%(combineFileName))
#    os.system("lsof %s"%(combineFileName) )
    os.system("combine -M Asymptotic -m %s -t %s --rMin=%f --rMax=%f %s %s"%(str(mass), str(asimov), rMin, rMax, fileName, blind))
    print rMax, " ", rMin

  print "Moving %s in the right directory"%(combineFileName)
  os.system("mv %s %s/%s"%(combineFileName, outDir, combineFileName))
  return

def writeCard(workspaceName, mass, fit_result_fileName, outDirName, setConst, isBiasCard=False, biasWS=None, bias_fileName=None, suffix=""):
#  obsRate = workspace.data("data_obs").sumEntries()
  #dataCardRootFileName = "datacard_m%d_VGamma.root" % (mass)
  dataCardFileName = "%s/datacard_m%d_VGamma%s.txt" % (outDirName,mass, suffix)
  

  ws_file = TFile.Open(fit_result_fileName)
  workspace = ws_file.Get(workspaceName)
  bkg_fileName = ""
  backgroundWS = RooWorkspace()
  bkg_functionName = ""
  pdf_index = ' '
  pdf_par = ' '
  last_bracket = [' ', ' ', ' ', ' ']
  last_bracket1 = [' ', ' ', ' ', ' ']
  last_bracket2 = [' ', ' ', ' ', ' ']
  if isBiasCard:
    bias_file = TFile.Open(bias_fileName)
    bkg_fileName = bias_fileName
    backgroundWS = bias_file.Get(biasWS)
    bkg_functionName = "multiPdf"
    pdf_index = 'pdf_index'
    pdf_par = 'discrete'
    last_bracket = [pdf_index, pdf_par, ' ', ' ']
  else:
    bkg_fileName = fit_result_fileName
    backgroundWS = workspace
    bkg_functionName = "bg_function"
    if not setConst:
      last_bracket  = [workspace.var("bg_p1").GetName(), 'flatParam', ' ', ' ']
      last_bracket1 = [workspace.var("bg_p2").GetName(), 'flatParam', ' ', ' ']
      last_bracket2 = [workspace.var("bg_p3").GetName(), 'flatParam', ' ', ' ']
                     
    else:
      last_bracket = [' ', ' ', ' ', ' ']

  
  signalRate = workspace.data("signal_exp_%d"%(mass)).sumEntries()
  divider = "-------------"
  dataCard = [['imax', '1', 'number of channels', ' ', ' '], 
            ['jmax', '*', 'number of backgrounds', ' ', ' '],                                                                                                                                                      #  dataCard = [['imax', '1', 'number of channels', ' ', ' '], 
            ['kmax', '*', 'number of systematic uncertainty sources', ' ', ' '],                                                                                                                                   #            ['jmax', '*', 'number of backgrounds', ' ', ' '], 
            [divider, ' ', ' ',' ', ' '],                                                                                                                                                                          #            ['kmax', '*', 'number of systematic uncertainty sources', ' ', ' '],
            ['shapes sig',workspace.GetName(),fit_result_fileName,workspace.GetName()+':'+workspace.pdf("signal_f").GetName(), workspace.GetName()+':'+workspace.pdf("signal_f").GetName()+"_$SYSTEMATIC"],        #            [divider, ' ', ' ',' ', ' '],
            ['shapes bkg', workspace.GetName(),bkg_fileName, backgroundWS.GetName()+':'+backgroundWS.pdf(bkg_functionName).GetName(),' '],                                                                         #            #['shapes bkg', workspace.GetName(),fit_result_fileName, workspace.GetName()+':'+workspace.pdf("bg_function").GetName(),' '],
            ['shapes data_obs', workspace.GetName(), fit_result_fileName, workspace.GetName()+':'+workspace.data("data_obs").GetName(), ' '],                                                                      #            #['shapes background', workspace.GetName(),fit_result_fileName, workspace.GetName()+':'+workspace.pdf("bg_function_normalised").GetName()],
            [divider, ' ', ' ',' ',' '],                                                                                                                                                                           #            ['shapes data_obs', workspace.GetName(), fit_result_fileName, workspace.GetName()+':'+workspace.data("data_obs").GetName(), ' '],
            ['### Observation', ' ', ' ',' ',' '],                                                                                                                                                                 #            [divider, ' ', ' ',' ',' '],
            ['bin', workspace.GetName(), ' ', ' ',' '],                                                                                                                                                            #            ['### Observation', ' ', ' ',' ',' '],
    #        ['observation',str(obsRate), ' ', ' '],                                                                                                                                                               #            ['bin', workspace.GetName(), ' ', ' ',' '],
            ['observation','-1', ' ', ' ',' '],                                                                                                                                                                    #    #        ['observation',str(obsRate), ' ', ' '],
            [divider, ' ', ' ',' ',' '],                                                                                                                                                                           #            ['observation','-1', ' ', ' ',' '],
            ['bin', ' ', workspace.GetName(), workspace.GetName()],                                                                                                                            #            [divider, ' ', ' ',' ',' '],
            ['process', ' ','sig', 'bkg'],                                                                                                                                                                 #            ['bin', ' ', workspace.GetName(), workspace.GetName(),workspace.GetName()],
            ['process', ' ','0', '1'],                                                                                                                                                                         #            ['process', ' ','sig', 'allgJ','bkg'],
            ['rate', ' ', '1','1'],                                                                                                                                                                            #            ['process', ' ','0', '1','2'],
            ['cms_lumi_13TeV', 'lnN', '1.027', '-'],                                                                                                                                                       #            ['rate', ' ', '1','1','1'],
            last_bracket,
            last_bracket1,
            last_bracket2,
            #[pdf_index, pdf_par, ' ', ' '],                                                                                                                                                                   #            ['cms_lumi_13TeV', 'lnN', '1.027', '1.027','-'],
    ]                                                                                                                                                                                                              #            #[workspace.function("dcb_norm").GetName(),  'flatParam',' ', ' ']
#    ]
#  dataCard = [['imax', '1', 'number of channels', ' '],
#          ['jmax', '1', 'number of backgrounds', ' '], 
#          ['kmax', '*', 'number of systematic uncertainty sources', ' '],
#          [divider, ' ', ' ',' '],
#          ['shapes signal',workspace.GetName(),fit_result_fileName,workspace.GetName()+':'+workspace.pdf("dcb").GetName()],
##          ['shapes mc',workspace.GetName(),fit_result_fileName,workspace.GetName()+':'+workspace.pdf("mc_function").GetName()],
#          #['shapes background', workspace.GetName(),fit_result_fileName, workspace.GetName()+':'+workspace.pdf("bg_function").GetName()],
#          #['shapes signal',workspace.GetName(),fit_result_fileName,workspace.GetName()+':'+workspace.pdf("dcb_normalised").GetName()],
#          ['shapes background', workspace.GetName(),fit_result_fileName, workspace.GetName()+':'+workspace.pdf("bg_function").GetName()],
#          ['shapes data_obs', workspace.GetName(), fit_result_fileName, workspace.GetName()+':'+workspace.data("data_obs").GetName()],#("data_obs").GetName()],
#          [divider, ' ', ' ',' '],
#          ['### Observation', ' ', ' ',' '],
#          ['bin', workspace.GetName(), ' ', ' '],
#  #        ['observation',str(obsRate), ' ', ' '],
#          ['observation','-1', ' ', ' '],
#          [divider, ' ', ' ',' '],
#          ['bin', ' ', workspace.GetName(), workspace.GetName()],
#          ['process', ' ','signal', 'background'],
#          ['process', ' ','0', '1'],
#          ['rate', ' ',str(signalRate) ,'1'],
#          ['cms_lumi_13TeV', 'lnN', '1.027', '-'],
#          #[workspace.var("bg_function_norm").GetName(),  'flatParam',' ', ' '],
#          #[workspace.function("dcb_norm").GetName(),  'flatParam',' ', ' ']
#  ]

  dataCardFileTxt = open(dataCardFileName,'w+')
  col_width = max(len(word) for row in dataCard for word in row) + 2  # padding
  for row in dataCard:
    print >> dataCardFileTxt,"".join(word.ljust(col_width) for word in row)

  dataCardFileTxt.close()
  return dataCardFileName

def getHistoInfo(histogram):

  xLow = histogram.GetBinCenter(1) - histogram.GetBinWidth(1)/2
  xHigh = histogram.GetBinCenter(histogram.GetNbinsX()) + histogram.GetBinWidth(histogram.GetNbinsX())/2

  return histogram.GetNbinsX(), xLow, xHigh


def convertToTH1F(rooH, nBins, xLow, xMax):
  x = 0
  y = 0
  histo = TH1F("histo", "histo", nBins, xLow, xMax)
  for i in range(0, nBins):
    #rooH.GetPoint(i,x,y)
    y = rooH.GetY()[i]#rooH.Eval(((xMax-xLow)*i/nBins)+ xLow)
    rightBin = (rooH.GetX()[i]-xLow)/((xMax-xLow)/nBins)+.5
    if rightBin not in range(0, nBins):#righBin>nBins || righBin<0 || rightBin==nan:
      continue
    histo.SetBinContent((int)(rightBin),y)#i+1,y)
    #print "bin: ", rightBin, " y ", y
  return histo

def myPrint(string, symbol="+"):
  frame=symbol
  for j in range(0, len(string)+7):
    frame=frame+symbol
  print "\n",frame
  print symbol," ",string," ",symbol
  print frame,"\n"


def importRooDataset(inFile, realVar, outName, treeName):
  print "Entering in importRooDataset function"
  ttree= inFile.Get(treeName)
  rooDS = RooDataSet(outName, outName, RooArgSet(realVar), RooFit.Import(ttree))

  return rooDS

def importRooDatahist(inFile, realVar, histoName, outName, rebin):
    print "Entering in importRooDatahist function"
    histogram = inFile.Get(histoName)
    histogram.Rebin(rebin)
    print "Getting: ", histoName, "histogram"
    dataHist = RooDataHist(outName, outName, RooArgList(realVar), RooFit.Import(histogram))
    nBinsNew, xLowNew, xHighNew = getHistoInfo(histogram)

    return dataHist, nBinsNew, xLowNew, xHighNew

def doUnbinnedFit(pdf, rooDS):
  print "Entering in doUnbinnedFit function"

def fillWorkspace(w, realVar):
  bg_p0 = RooRealVar ("bg_p0", "bg_p0", 0., -1000., 200.)
#  bg_p1 = RooRealVar ("bg_p1", "bg_p1", -13, -15, -11.)
#  bg_p2 = RooRealVar ("bg_p2", "bg_p2", -1.4, -2, -1.)
#  print "Dijet function chosen"
#  pdf = RooGenericPdf("bg_function","(pow(@0/13000,@1+@2*log(@0/13000)))",RooArgList(realVar,bg_p1,bg_p2))
#
#  #w.RooWorkspace::import("bg_function")#pdf)
#  rootTools.Utils.importToWS(w,pdf)




def doBinnedFit(pdf, rooDH, minFit, maxFit):
  print "Entering in doBinnedFit function"
  fr = pdf.fitTo(rooDH)
  ####nll = pdf.createNLL(rooDH,RooFit.Range(minFit, maxFit),RooFit.Extended(True),RooFit.Offset(True))
  ####m2 = RooMinimizer(nll)
  ####m2.setStrategy(2)
  ####m2.setMaxFunctionCalls(100000)
  ####m2.setMaxIterations(100000)
  ####migrad_status = m2.minimize('Minuit2','migrad')
  ####improve_status = m2.minimize('Minuit2','improve')
  ####hesse_status = m2.minimize('Minuit2','hesse')
  ####minos_status = m2.minos()
  ####if hesse_status != 3 :
  ####    hesse_status = m2.minimize('Minuit2','hesse')
  ####fr = m2.save()

  return fr

def createFitFuntion(value, realVar):
  #Dictionary
  # 1: dijet function
  # 2: diphoton function
  
  print "Entering in creatingFitFuntion function"
  if value==1:
    bg_p0 = RooRealVar ("bg_p0", "bg_p0", 0., -1000., 200.)
    bg_p1 = RooRealVar ("bg_p1", "bg_p1", -13, -15, -11.)
    bg_p2 = RooRealVar ("bg_p2", "bg_p2", -1.4, -2, -1.)
    print "Dijet function chosen"
    bg_function = RooGenericPdf("bg_function","(pow(@0/13000,@1+@2*log(@0/13000)))",RooArgList(realVar,bg_p1,bg_p2))

  return bg_function
  


if __name__ == '__main__':
  
  startProgram = time.time()
  parser = OptionParser()
#  parser.add_option('-u','--unbinned',dest="unbinned",type="string",default="input.root",
#      help="input file name of the unbinned dataset")
  parser.add_option('-b','--binned',dest="binned",default=False, action='store_true',
      help="if binned or unbinned fit")
  parser.add_option('-H', '--histoName', dest="histoName", type="string", default="dataHisto",
      help="binned histogram name")
  parser.add_option('-i', '--inputFileName', dest="inputFileName", type="string",default="input.root",
      help="input dataset file name")
  parser.add_option('--mh', dest="res_mass", type=int, default= 120,
      help="resonance mass")
  parser.add_option('--is', dest="signalFileName", type="string",default="input.root",
      help="signal input file name")
  parser.add_option('-o', '--outDir', dest='outDir', type="string", default="outDir",
      help="output directory")
  parser.add_option('--runCombine', dest="runCombine", default=False, action='store_true',
      help="second step of the program. After the fit, relaunch the script with this option to run combine")
  parser.add_option('--asimov', dest="doAsimovOnly", default=False, action='store_true',
      help="throw only the asimov to extract the limit")
  parser.add_option('--blind', dest="isBlinded", default=False, action='store_true',
      help="run blinded or unblinded")
  parser.add_option('--doLimitScan', dest="doLimitScan", default=False, action='store_true',
      help="run a manual limit scan")
  parser.add_option('--lowL', dest="lowL", type="float", default=0,
      help="low limit for manual scan")
  parser.add_option('--highL', dest="highL", type="float", default=10,
      help="high limit for manual scan")
  parser.add_option('--nSteps', dest="nSteps", type="int", default=10,
      help="steps for manual limit scan")
  parser.add_option('--nToys', dest="nToys", type="int", default=100,
      help="n toys for bias test")
  parser.add_option('--doBias', dest="doBias", default=False, action='store_true',
      help="do bias test")
  parser.add_option('--nExpected', dest="nExpt", type="float", default=1.,
      help="n expected for bias toys")
  parser.add_option('--treeName', dest="treeName", type="string", default="mio",
      help="roodataset tree name")
  parser.add_option('--setConstant', dest="setConst", default=False, action='store_true',
      help="set all the fit parameters constant. If false, the datacard will have the bkg parameters set as FLATPARAM")
  parser.add_option('--mcHisto', dest="mcHisto", type="string", default="allBkgHisto",
      help="histogram name of MC background")
  parser.add_option('--sigHisto', dest="sigHisto", type="string", default="signalHisto_",
      help="histogram name of signal shape")
  parser.add_option('--rVarName', dest="rVarName", type="string", default="ak08_photon_mass",
       help="tree leaf/branch name of the invariant mass")
  parser.add_option('--rangeVarHigh', dest="rangeVarHigh", type="float", default="4000.",
       help="realvar high limit range")
  parser.add_option('--rangeVarLow', dest="rangeVarLow", type="float", default="600.",
       help="reavar low limit range")
  parser.add_option('--nDivisions', dest="nDivisions", type="int", default="40",
       help="bins width for plots")

  (options,args) = parser.parse_args()

  
  sigHistoName = options.sigHisto+str(options.res_mass)
  outDirName=options.outDir+'_m'+str(options.res_mass)
  if options.runCombine:
    dataCardName = "%s/datacard_m%d_VGamma.txt" % (outDirName,options.res_mass)
    runCombine(options.res_mass, dataCardName.replace('.txt','.root'), outDirName, options.doAsimovOnly, options.isBlinded)
    endProgram = time.time()
    myPrint("Time elapsed: %f s"%(endProgram-startProgram))
    sys.exit("Completing combine command")

  if options.doLimitScan:
    dataCardName = "%s/datacard_m%d_VGamma.txt" % (outDirName,options.res_mass)
    doManualLimitScan(options.nSteps, options.lowL, options.highL, options.res_mass, outDirName, dataCardName)
    endProgram = time.time()
    myPrint("Time elapsed: %f s"%(endProgram-startProgram))
    sys.exit("Completing manual scan")

    
  fileInput = TFile.Open(options.inputFileName)

  os.system('mkdir -p %s'%(outDirName))
  bins = (options.rangeVarHigh-options.rangeVarLow)/options.nDivisions

  myPrint("Analysing background")
  if not fileInput:
    print "File ", options.inputFileName, " not found!"
    sys.exit("File not found")

  mVg = RooRealVar(options.rVarName, options.rVarName, options.rangeVarLow, options.rangeVarHigh)

  if options.binned:
    #print "Found binned file"
    dataHist = importRooDatahist(fileInput, mVg, options.histoName)
    data_val = dataHist
  else:
    dataSet = importRooDataset(fileInput, mVg, "data_obs", options.treeName)
    data_val = dataSet


  if options.doBias:
    dataCardName = "%s/datacard_m%d_VGamma.txt" % (outDirName,options.res_mass)
    doBias("WS_Vg", "%s/WS_Vg_%i.root"%(outDirName,options.res_mass), outDirName, dataCardName, options.res_mass, options.rangeVarLow, options.rangeVarHigh, mVg, data_val, options.nToys, options.nExpt)
    endProgram = time.time()
    myPrint("Time elapsed: %f s"%(endProgram-startProgram))
    sys.exit("Completing bias test")

  #WS_Vg.Print("v")
  ##fillWorkspace(WS_Vg, mVg)

#  pdf = createFitFuntion(1, mVg)

  bg_function_norm = RooRealVar ("bg_function_norm", "bg_function_norm", 10., 1., 2000000.)
  #bg_p1 = RooRealVar ("bg_p1", "bg_p1", -13, -15, -11.)
  #bg_p2 = RooRealVar ("bg_p2", "bg_p2", -1.4, -2, -1.)
  bg_p1 = RooRealVar("bg_p1", "bg_p1", -13., -15., -9.)
  bg_p2 = RooRealVar("bg_p2", "bg_p2", -1.4, -2., -1.)
  bg_p3 = RooRealVar("bg_p3", "bg_p3", .2,0.01,10.)
  print "Dijet function chosen"
#  pdf = RooGenericPdf("bg_function","@3*(pow(@0/13000,@1+@2*log(@0/13000)))",RooArgList(mVg,bg_p1,bg_p2, bg_function_norm))
  bg_function = RooGenericPdf("bg_function","(pow((1-(@0/13000)),@3)*pow(@0/13000,(@1+@2*log(@0/13000))))", RooArgList(mVg,bg_p1,bg_p2, bg_p3))
  #fr = binnedFit(extDijetPdf,dataHist,sideband,options.useWeight)       
  bg_function_normalised = RooExtendPdf("bg_function_normalised","bg_function_normalised", bg_function, bg_function_norm)

  dataBinned = RooDataHist()
  dataBinned.SetName("dataBinned")#,"dataBinned", mVg)
  if options.binned:
    dataBinned = data_val
  else:
    dataBinned = data_val.binnedClone("dataBinned", mVg.GetName())#"dataBinned")
  histoData = TH1F()
  histoData.SetName("histoData")#,"histoData")
  nBins, xLow, xHigh = getHistoInfo(dataBinned.createHistogram("histoData", mVg))

  rangeLow, rangeHigh = getGoodRange(dataBinned, mVg, True)

  fr = bg_function_normalised.fitTo(data_val, RooFit.Range(rangeLow, rangeHigh))#, RooFit.Strategy(2))
  bg_function.fitTo(data_val, RooFit.Range(rangeLow, rangeHigh))#, RooFit.Strategy(2))


#  pdf.Print("v")
  frameBkg = mVg.frame(RooFit.Bins(bins))
  frameBkgPull = mVg.frame(RooFit.Bins(bins))
  
  #rootTools.Utils.importToWS(WS_Vg,fr)
#  fr = doBinnedFit(pdf, dataSet, 600, 3000)
  data_val.plotOn(frameBkg)
  bg_function.plotOn(frameBkg)#, RooFit.VisualizeError(fr,2), RooFit.FillColor(kRed))
  nset =  RooArgSet(mVg) 
  #print "----------> ", pdf.getVal(nset)
  #print "----------> ", pdf.getVal()
  bkgPull = frameBkg.pullHist()
  pullBkgHisto = convertToTH1F(bkgPull, bins, xLow, xHigh)#nBins, xLow, xHigh)
  frameBkgPull.addPlotable(bkgPull,"H")#, SetFillColor(kOrange), SetLineWidth(2), SetLineColor(kBlack))

  c_bkg = TCanvas("c_bkg","c_bkg",1)
  xPad = 0.3
  p_1Bkg = TPad("p_1Bkg", "Bkg plot LOG", 0, xPad-.01, 1, 1)
  p_1Bkg.SetFillStyle(4000)
  p_1Bkg.SetFrameFillColor(0)
  p_1Bkg.SetBottomMargin(0.02)
  p_1Bkg.SetTopMargin(0.06)

  p_2Bkg = TPad("p_2Bkg", "Bkg pull",0,0,1,xPad)
  p_2Bkg.SetBottomMargin((1.-xPad)/xPad*0.13)
  p_2Bkg.SetTopMargin(0.03)
  p_2Bkg.SetFillColor(0)
  p_2Bkg.SetBorderMode(0)
  p_2Bkg.SetBorderSize(2)
  p_2Bkg.SetFrameBorderMode(0)
  p_2Bkg.SetFrameBorderMode(0)

  p_1Bkg.Draw()
  p_2Bkg.Draw()

  p_1Bkg.cd()
  gPad.SetLogy()
  frameBkg.Draw()
  p_2Bkg.cd()
  #frameBkgPull.Draw()
  gPad.SetGridy()
  pullBkgHisto.SetLineColor(kBlack)
  pullBkgHisto.SetLineWidth(1)
  pullBkgHisto.SetFillColor(kOrange+8)
  bkgPull.SetFillColor(kWhite)
  bkgPull.SetLineColor(kWhite)
  bkgPull.SetMarkerColor(kWhite)
  bkgPull.GetYaxis().SetLabelSize(.08)
  bkgPull.GetXaxis().SetLabelSize(.08)
  bkgPull.Draw("ABX")
  pullBkgHisto.Draw("SAME")
  c_bkg.SaveAs("%s/bkg_fit_Log.png"%(outDirName))

  c_bkg_1 = TCanvas("c_bkg_1","c_bkg_1",1)
  p_1Bkg_1 = TPad("p_1Bkg_1", "Bkg plot", 0, xPad-.01, 1, 1)
  p_1Bkg_1.SetFillStyle(4000)
  p_1Bkg_1.SetFrameFillColor(0)
  p_1Bkg_1.SetBottomMargin(0.02)
  p_1Bkg_1.SetTopMargin(0.06)
  
  p_2Bkg_1 = TPad("p_2Bkg_1", "Bkg pull",0,0,1,xPad)
  p_2Bkg_1.SetBottomMargin((1.-xPad)/xPad*0.13)
  p_2Bkg_1.SetTopMargin(0.03)
  p_2Bkg_1.SetFillColor(0)
  p_2Bkg_1.SetBorderMode(0)
  p_2Bkg_1.SetBorderSize(2)
  p_2Bkg_1.SetFrameBorderMode(0)
  p_2Bkg_1.SetFrameBorderMode(0)
  
  p_1Bkg_1.Draw()
  p_2Bkg_1.Draw()
  
  p_1Bkg_1.cd()
  frameBkg.Draw()
  p_2Bkg_1.cd()
  #frameBkgPull.Draw()
  gPad.SetGridy()
  pullBkgHisto.SetLineColor(kBlack)
  pullBkgHisto.SetLineWidth(1)
  pullBkgHisto.SetFillColor(kOrange+8)
  bkgPull.SetFillColor(kWhite)
  bkgPull.SetLineColor(kWhite)
  bkgPull.SetMarkerColor(kWhite)
  bkgPull.GetYaxis().SetLabelSize(.08)
  bkgPull.GetXaxis().SetLabelSize(.08)
  bkgPull.Draw("ABX")
  pullBkgHisto.Draw("SAME")
  c_bkg_1.SaveAs("%s/bkg_fit.png"%(outDirName))

  
  myPrint("Analysing signal")
  print "Opening ", options.signalFileName, " signal file"
  signalFileInput = TFile.Open(options.signalFileName)
  if not signalFileInput:
    print "Warning, file: ", options.signalFileName, " not opened!"
  signalHist, nBins, xLow, xHigh = importRooDatahist(signalFileInput, mVg, sigHistoName, "signal_exp_"+str(options.res_mass), 20)
  rangeLow=TMath.Max(600., options.res_mass-0.3*options.res_mass)
  rangeHigh=TMath.Min(3600., options.res_mass+0.3*options.res_mass)

#  sg_p0 = RooRealVar("sg_p0", "sg_p0", options.res_mass, options.res_mass-0.1*options.res_mass, options.res_mass+0.1*options.res_mass)
#  MH      = RooRealVar("MH", "Double CB Bias", options.res_mass, options.res_mass-0.1*options.res_mass, options.res_mass+0.1*options.res_mass)
#  sg_p1 = RooRealVar("sg_p1", "sg_p1", 0.03*options.res_mass, 5., 400.)
#  sg_p2 = RooRealVar("sg_p2", "sg_p2", 1.3, 0., 200.)
#  sg_p3 = RooRealVar("sg_p3", "sg_p3", 5, 0., 300.)
#  sg_p4 = RooRealVar("sg_p4", "sg_p4", options.res_mass, 500., 800.)
#  sg_p5 = RooRealVar("sg_p5", "sg_p5", 150, 0., 3000.)
#  sg_p6 = RooRealVar("sg_p6", "sg_p6", 0.99, 0.,1.)
#    
#  signalCore = RooCBShape("signalCore", "signalCore", mVg, MH, sg_p1,sg_p2, sg_p3)
#  signalComb = RooGaussian("signalComb", "Combinatoric", mVg, MH, sg_p5)
#  norm1     = RooRealVar("norm1", "norm parameters", .1, .01, 10.)
#  norm2     = RooRealVar("norm2", "norm parameters", .1, .01, 10.)
#  norm3     = RooRealVar("norm3", "norm parameters", .1, .01, 10.)
#  signal_normalised_norm = RooFormulaVar("signal_normalised_norm","@0+@1*@2+@3*@3*@2", RooArgList(norm1, norm2, MH, norm3))
#  signal = RooAddPdf("signal", "signal", RooArgList(signalCore, signalComb), RooArgList(sg_p6))
#  signal_normalised = RooExtendPdf("signal_normalised", "signal_normalised", signal, signal_normalised_norm)

  mean      = RooRealVar("mean", "Double CB Bias", options.res_mass, options.res_mass-0.1*options.res_mass, options.res_mass+0.1*options.res_mass)
#  MH      = RooRealVar("MH", "Double CB Bias", options.res_mass, options.res_mass-0.1*options.res_mass, options.res_mass+0.1*options.res_mass)
  sigma     = RooRealVar("sigma", "Double CB Width", 0.03*options.res_mass, 5., 400.)
  dCBCutL   = RooRealVar("dCBCutL", "Double CB Cut left", 1., 0.1, 50.)
  dCBCutR   = RooRealVar("dCBCutR", "Double CB Cut right",1., 0.1, 50.)
  dCBPowerL = RooRealVar("dCBPowerL", "Double CB Power left", 2., 0.2, 50.)
  dCBPowerR = RooRealVar("dCBPowerR", "Double CB Power right",2., 0.2, 50.)

  dcb = RooDCBShape("dcb", "double crystal ball", mVg, mean, sigma, dCBCutL, dCBCutR, dCBPowerL, dCBPowerR)
  #dcb_normalised_norm = RooRealVar("dcb_normalised_norm", "Double Crystal Ball normalization", 10, 1, 20000)
  norm1     = RooRealVar("norm1", "norm parameters", .1, .01, 10.)
  norm2     = RooRealVar("norm2", "norm parameters", .1, .01, 10.)
  norm3     = RooRealVar("norm3", "norm parameters", .1, .01, 10.)
  dcb_norm = RooFormulaVar("dcb_norm","@0+@1*@2+@3*@3*@2", RooArgList(norm1, norm2, mean, norm3))
#  dcb_norm = RooRealVar("dcb_norm", "dcb_norm",10, 1, 20000000)

  dcb_normalised = RooExtendPdf("dcb_normalised", "dcb_normalised", dcb, dcb_norm)

  #dcb_normalised = RooAddPdf("dcb_normalised",dcb, dcb_norm)

  fr1 = dcb_normalised.fitTo(signalHist, RooFit.Range(rangeLow, rangeHigh))#, RooFit.Strategy(2))
  dcb.fitTo(signalHist, RooFit.Range(rangeLow, rangeHigh))#, RooFit.Strategy(2))

  frameSgn = mVg.frame(RooFit.Bins(bins))
#  frameSgn.SetMaximum(3000)
#  frameSgn.SetMinimum(0.01)
  frameSgnPull = mVg.frame(RooFit.Bins(bins))
  
#  print dcb.getVal()
#  dcb_norm = dcb.createIntegral(mVg)
  MH      = RooRealVar("MH", "Double CB Bias", mean.getVal())
#  dcb.fitTo(signalHist, RooFit.Range(rangeLow, rangeHigh))
#  fr1 = dcb_normalised.fitTo(signalHist, RooFit.Range(rangeLow, rangeHigh))

  sigma_signal     = RooRealVar("sigma_signal", "Double CB Width", sigma.getVal())
  dCBCutL_signal   = RooRealVar("dCBCutL_signal", "Double CB Cut left", dCBCutL.getVal())
  dCBCutR_signal   = RooRealVar("dCBCutR_signal", "Double CB Cut right",dCBCutR.getVal())
  dCBPowerL_signal = RooRealVar("dCBPowerL_signal", "Double CB Power left", dCBPowerL.getVal())
  dCBPowerR_signal = RooRealVar("dCBPowerR_signal", "Double CB Power right",dCBPowerR.getVal())
  norm1_signal     = RooRealVar("norm1_signal", "norm parameters", norm1.getVal())
  norm2_signal     = RooRealVar("norm2_signal", "norm parameters", norm2.getVal())
  norm3_signal     = RooRealVar("norm3_signal", "norm parameters", norm3.getVal())
  signal_f_norm = RooFormulaVar("signal_f_norm","@0+@1*@2+@3*@3*@2", RooArgList(norm1_signal, norm2_signal, MH, norm3_signal))

  signal_f = RooDCBShape("signal_f", "signal double crystal ball", mVg, MH, sigma_signal, dCBCutL_signal, dCBCutR_signal, dCBPowerL_signal, dCBPowerR_signal)

  signalHist.plotOn(frameSgn)
  #dcb.plotOn(frameSgn)#, RooFit.VisualizeError(fr1))
  signal_f.plotOn(frameSgn, RooFit.LineColor(2))
  sgnPull = frameSgn.pullHist()
  frameSgnPull.addPlotable(sgnPull,"H")#, SetFillColor(kOrange), SetLineWidth(2), SetLineColor(kBlack))

  c_sgn = TCanvas("c_sgn","c_sgn",1)
  p_1Sgn = TPad("p_1Sgn", "Sgn plot Log", 0, xPad-.01, 1, 1)
  p_1Sgn.SetFillStyle(4000)
  p_1Sgn.SetFrameFillColor(0)
  p_1Sgn.SetBottomMargin(0.02)
  p_1Sgn.SetTopMargin(0.06)

  p_2Sgn = TPad("p_2Sgn", "Sgn pull",0,0,1,xPad)
  p_2Sgn.SetBottomMargin((1.-xPad)/xPad*0.13)
  p_2Sgn.SetTopMargin(0.03)
  p_2Sgn.SetFillColor(0)
  p_2Sgn.SetBorderMode(0)
  p_2Sgn.SetBorderSize(2)
  p_2Sgn.SetFrameBorderMode(0)
  p_2Sgn.SetFrameBorderMode(0)
                                                                                                                          
  p_1Sgn.Draw()
  p_2Sgn.Draw()

  p_1Sgn.cd()
  gPad.SetLogy() 
  fakeH = TH1D("fakeH","",bins, xLow, xHigh) #when having negative bins in the signal histogram (it's reweighted)
  fakeH.GetYaxis().SetRangeUser(0.1,3000)
  fakeH.Draw()
  frameSgn.Draw("SAME")

  p_2Sgn.cd()
  #frameSgnPull.Draw()
  gPad.SetGridy()
  myHisto = convertToTH1F(sgnPull, bins, xLow, xHigh)
  myHisto.SetLineColor(kBlack)
  myHisto.SetLineWidth(1)
  myHisto.SetFillColor(kOrange+8)
  sgnPull.Draw("ABX")
#  sgnPull.GetXaxis().SetTitle("inv. mass [GeV]")
#  sgnPull.GetYaxis().SetTitle("pull")
#  sgnPull.SetFillColor(kWhite)
#  sgnPull.SetLineColor(kWhite)
#  sgnPull.SetMarkerColor(kWhite)
#  sgnPull.SetMarkerStyle(21)
#  sgnPull.Draw("ABX")
#  myHisto.Draw("histosame")
#  
  c_sgn.SaveAs("%s/sgn_fit_Log.png"%(outDirName))

  c_sgn_1 = TCanvas("c_sgn_1","c_sgn",1)
  p_1Sgn_1 = TPad("p_1Sgn_1", "Sgn plot", 0, xPad-.01, 1, 1)
  p_1Sgn_1.SetFillStyle(4000)
  p_1Sgn_1.SetFrameFillColor(0)
  p_1Sgn_1.SetBottomMargin(0.02)
  p_1Sgn_1.SetTopMargin(0.06)
  
  p_2Sgn_1 = TPad("p_2Sgn_1", "Sgn pull",0,0,1,xPad)
  p_2Sgn_1.SetBottomMargin((1.-xPad)/xPad*0.13)
  p_2Sgn_1.SetTopMargin(0.03)
  p_2Sgn_1.SetFillColor(0)
  p_2Sgn_1.SetBorderMode(0)
  p_2Sgn_1.SetBorderSize(2)
  p_2Sgn_1.SetFrameBorderMode(0)
  p_2Sgn_1.SetFrameBorderMode(0)
  
  p_1Sgn_1.Draw()
  p_2Sgn_1.Draw()
  
  p_1Sgn_1.cd()
  frameSgn.Draw()
  p_2Sgn_1.cd()
  #frameSgnPull.Draw()
  gPad.SetGridy()
  sgnPull.Draw("ABX")
  myHisto.Draw("histosame")
  
  c_sgn_1.SaveAs("%s/sgn_fit.png"%(outDirName))


  mc_function_norm = RooRealVar ("mc_function_norm", "mc_function_norm", 10., 1., 2000000.)
  #mc_p1 = RooRealVar ("mc_p1", "mc_p1", -13, -15, -11.)
  #mc_p2 = RooRealVar ("mc_p2", "mc_p2", -1.4, -2, -1.)
  mc_p1 = RooRealVar("mc_p1", "mc_p1", -13., -15., -9.)
  mc_p2 = RooRealVar("mc_p2", "mc_p2", -1.4, -2., -1.)
  mc_p3 = RooRealVar("mc_p3", "mc_p3", .2,0.01,10.)
  print "Dijet function chosen"
#  mc = RooGenericPdf("mc_function","@3*(pow(@0/13000,@1+@2*log(@0/13000)))",RooArgList(mVg,mc_p1,mc_p2, mc_function_norm))
  mc = RooGenericPdf("mc_function","(pow((1-(@0/13000)),@3)*pow(@0/13000,(@1+@2*log(@0/13000))))", RooArgList(mVg,mc_p1,mc_p2, mc_p3))
  #fr = binnedFit(extDijetPdf,dataHist,sideband,options.useWeight)       
  mc_normalised = RooExtendPdf("mc_function_normalised","mc_function_normalised", mc, mc_function_norm)#alised_norm)

  mcHist, nBins, xLow, xHigh = importRooDatahist(signalFileInput, mVg, options.mcHisto, "mc_exp", 20)
  rangeLow, rangeHigh = getGoodRange(mcHist, mVg, True)

  fr2 = mc_normalised.fitTo(mcHist, RooFit.Range(rangeLow, rangeHigh))#, RooFit.Strategy(2))
  if options.setConst:
    mc_p1.setConstant(True)
    mc_p2.setConstant(True)
    mc_p3.setConstant(True)
  mc.fitTo(mcHist, RooFit.Range(rangeLow, rangeHigh))#, RooFit.Strategy(2))

  frameMc = mVg.frame(RooFit.Bins(bins))
  frameMcPull = mVg.frame(RooFit.Bins(bins))
   
  mcHist.plotOn(frameMc)
  mc.plotOn(frameMc)#, RooFit.VisualizeError(fr2))
  mcPull = frameMc.pullHist()
  frameMcPull.addPlotable(mcPull,"H")#, SetFillColor(kOrange), SetLineWidth(2), SetLineColor(kBlack))


  c_mc = TCanvas("c_mc","c_mc",1)
  p_1Mc = TPad("p_1Mc", "Mc plot Log", 0, xPad-.01, 1, 1)
  p_1Mc.SetFillStyle(4000)
  p_1Mc.SetFrameFillColor(0)
  p_1Mc.SetBottomMargin(0.02)
  p_1Mc.SetTopMargin(0.06)

  p_2Mc = TPad("p_2Mc", "Mc pull",0,0,1,xPad)
  p_2Mc.SetBottomMargin((1.-xPad)/xPad*0.13)
  p_2Mc.SetTopMargin(0.03)
  p_2Mc.SetFillColor(0)
  p_2Mc.SetBorderMode(0)
  p_2Mc.SetBorderSize(2)
  p_2Mc.SetFrameBorderMode(0)
  p_2Mc.SetFrameBorderMode(0)
                                                                                                                          
  p_1Mc.Draw()
  p_2Mc.Draw()

  p_1Mc.cd()
#  gPad.SetLogy()
  frameMc.Draw()
  p_2Mc.cd()
  #frameMcPull.Draw()
  gPad.SetGridy()
  mcPullHisto = convertToTH1F(mcPull, nBins, xLow, xHigh)
  mcPullHisto.SetLineColor(kBlack)
  mcPullHisto.SetLineWidth(1)
  mcPullHisto.SetFillColor(kOrange+8)
  mcPull.GetXaxis().SetTitle("inv. mass [GeV]")
  mcPull.GetYaxis().SetTitle("pull")
  mcPull.SetFillColor(kWhite)
  mcPull.SetLineColor(kWhite)
  mcPull.SetMarkerColor(kWhite)
  mcPull.SetMarkerStyle(21)
  mcPull.Draw("ABX")
  mcPullHisto.Draw("histosame")
  
  c_mc.SaveAs("%s/mc_fit_Log.png"%(outDirName))

  c_mc_1 = TCanvas("c_mc_1","c_mc",1)
  p_1Mc_1 = TPad("p_1Mc_1", "Mc plot", 0, xPad-.01, 1, 1)
  p_1Mc_1.SetFillStyle(4000)
  p_1Mc_1.SetFrameFillColor(0)
  p_1Mc_1.SetBottomMargin(0.02)
  p_1Mc_1.SetTopMargin(0.06)
  
  p_2Mc_1 = TPad("p_2Mc_1", "Mc pull",0,0,1,xPad)
  p_2Mc_1.SetBottomMargin((1.-xPad)/xPad*0.13)
  p_2Mc_1.SetTopMargin(0.03)
  p_2Mc_1.SetFillColor(0)
  p_2Mc_1.SetBorderMode(0)
  p_2Mc_1.SetBorderSize(2)
  p_2Mc_1.SetFrameBorderMode(0)
  p_2Mc_1.SetFrameBorderMode(0)
  
  p_1Mc_1.Draw()
  p_2Mc_1.Draw()
  
  p_1Mc_1.cd()
  frameMc.Draw()
  p_2Mc_1.cd()
  #frameMcPull.Draw()
  gPad.SetGridy()
  mcPull.Draw("ABX")
  mcPullHisto.Draw("histosame")
  
  c_mc.SaveAs("%s/mc_fit.png"%(outDirName))


#  print "----------- ", dcb_norm.getVal()
#  dcb_norm.setVal(10.)
  bg_function_norm.setVal(dataBinned.sumEntries())
  mc_function_norm.setVal(mcHist.sumEntries())

  if options.setConst:
  #  MH.setConstant(True)
  #  sigma.setConstant(True)
  #  dCBCutL.setConstant(True)
  #  dCBCutR.setConstant(True)
  #  dCBPowerL.setConstant(True)
  #  dCBPowerR.setConstant(True)
  #  norm1.setConstant(True)
  #  norm2.setConstant(True)
  #  norm3.setConstant(True)
  ##  bg_function_norm.setConstant(True)
  ##  mc_function_norm.setConstant(True)
    bg_p1.setConstant(True)
    bg_p2.setConstant(True)
    bg_p3.setConstant(True)
    mc_p1.setConstant(True)
    mc_p2.setConstant(True)
    mc_p3.setConstant(True)




  WS_Vg = RooWorkspace("WS_Vg")#, kTRUE)
  rootTools.Utils.importToWS(WS_Vg, mVg)
  rootTools.Utils.importToWS(WS_Vg, data_val)
  rootTools.Utils.importToWS(WS_Vg, bg_function)#_normalised)
  rootTools.Utils.importToWS(WS_Vg, signalHist)
  rootTools.Utils.importToWS(WS_Vg, dcb)#signal_normalised)
  rootTools.Utils.importToWS(WS_Vg, dcb_norm)
  rootTools.Utils.importToWS(WS_Vg, bg_function_norm)
  rootTools.Utils.importToWS(WS_Vg, signal_f)
  rootTools.Utils.importToWS(WS_Vg, signal_f_norm)
#  rootTools.Utils.importToWS(WS_Vg, mcHist)
#  rootTools.Utils.importToWS(WS_Vg, mc)
#  rootTools.Utils.importToWS(WS_Vg, mc_function_norm)
  WS_Vg.Print("v")

#  WS_Vg_outFile = TFile.Open("WS_Vg_%i.root"%(options.res_mass), "recreate")
  WS_Vg.SaveAs("%s/WS_Vg_%i.root"%(outDirName,options.res_mass))#Write()


  myPrint("Closing input and workspace files")
  fileInput.Close()
  signalFileInput.Close()
#  WS_Vg_outFile.Close()

  dataCardName = writeCard(WS_Vg.GetName(), options.res_mass, "%s/WS_Vg_%i.root"%(outDirName,options.res_mass), outDirName, options.setConst) 
  writeRootDatacard(dataCardName, options.res_mass)
#  runCombine(options.res_mass, dataCardName)

  endProgram = time.time()

  myPrint("Time elapsed: %f s"%(endProgram-startProgram))
