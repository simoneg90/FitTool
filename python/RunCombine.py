from optparse import OptionParser
import os
import ROOT as rt
from array import *
from framework import Config
import sys
import glob
import rootTools
import time

NSIGMA = 10.0

def massIterable(massList):    
    if len(massList.split(','))==1:
        massIterableList = [massList]
    else:
        massIterableList = list(eval(massList))
    print "Masses list: ", massIterableList
    return massIterableList
        
def exec_me(command,dryRun=True):
    print command
    if not dryRun: os.system(command)
        
def writeBashScript(options,massPoint,iJob=0):
        
    lumiFloat = [float(lumiStr) for lumiStr in options.lumi.split('_')]    
    lumiTotal = sum(lumiFloat)
    
    submitDir = options.outDir
    massPoint = str(massPoint)

    
    signalSys = ''
    if options.noSignalSys:
        signalSys = '--no-signal-sys'
    elif options.noSys:
        signalSys = '--no-sys'
        
    penaltyString = ''
    if options.penalty:
        penaltyString = '--penalty'

    decoString = ''
    if options.deco:
        decoString  ='--deco'
        
    bayesString = ''
    if options.bayes:
        bayesString  ='--bayes'
        
    toyString = ''
    if options.toys>-1:
        toyString  ='--toys %i'%options.toys

    xsecString = '--xsec %f'%options.xsec

    signifString = ''
    if options.signif:
        signifString = '--signif'
        
    # prepare the script to run
    outputname = submitDir+"/submit_"+options.model+"_"+massPoint+"_lumi-%.3f_"%(lumiTotal)+options.box+"_%i"%(iJob)+".src"
        
    ffDir = submitDir+"/logs_"+options.model+"_"+massPoint+"_"+options.box+"_%i"%(iJob)
    user = os.environ['USER']
    pwd = os.environ['PWD']
        
    if options.noSys:
        combineDir = "/afs/cern.ch/work/%s/%s/DIJET/Limits/%s_nosys/"%(user[0],user,options.model) # directory where combine output files will be copied
    else:        
        combineDir = "/afs/cern.ch/work/%s/%s/DIJET/Limits/%s/"%(user[0],user,options.model) # directory where combine output files will be copied
    cmsswBase = "/afs/cern.ch/work/%s/%s/DIJET/CMSSW_7_4_14"%(user[0],user) # directory where 'cmsenv' will be run (needs to have combine setup)
    
    script =  '#!/usr/bin/env bash -x\n'
    script += 'mkdir -p %s\n'%combineDir        
    script += 'echo $SHELL\n'
    script += 'pwd\n'
    script += 'cd %s/src/CMSDIJET/DijetRootTreeAnalyzer \n'%(cmsswBase)
    script += 'pwd\n'
    script += "export SCRAM_ARCH=slc6_amd64_gcc491\n"
    script += "export CMSSW_BASE=%s\n"%(cmsswBase)
    script += 'eval `scramv1 runtime -sh`\n'
    script += 'cd - \n'
    script += "export TWD=${PWD}/%s_%s_lumi-%.3f_%s\n"%(options.model,massPoint,lumiTotal,options.box)
    script += "mkdir -p $TWD\n"
    script += "cd $TWD\n"
    script += 'pwd\n'
    script += 'git clone git@github.com:CMSDIJET/DijetRootTreeAnalyzer CMSDIJET/DijetRootTreeAnalyzer\n'
    script += 'cd CMSDIJET/DijetRootTreeAnalyzer\n'
    script += 'git checkout -b Limits %s\n'%(options.tag)
    script += 'mkdir -p %s\n'%submitDir    
    if 'CaloDijet2015' in options.box.split('_') or options.box=='CaloDijet20152016':
        script += 'wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_%s_13TeV_CaloScouting_Spring15.root -P inputs/\n'%(options.model)
        for sys in ['JERUP','JERDOWN','JESUP','JESDOWN']:
            script += 'wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_%s_13TeV_CaloScouting_Spring15_%s.root -P inputs/\n'%(options.model,sys)            
    if 'CaloDijet2016' in options.box.split('_'):        
        script += 'wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_%s_13TeV_CaloScouting_Spring16.root -P inputs/\n'%(options.model)
        for sys in ['JERUP','JERDOWN','JESUP','JESDOWN']:
            script += 'wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_%s_13TeV_CaloScouting_Spring16_%s.root -P inputs/\n'%(options.model,sys)            
    if 'PFDijet2016' in options.box.split('_'):
        script += 'wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_%s_13TeV_Spring16.root -P inputs/\n'%(options.model)
        for sys in ['JERUP','JESUP','JESDOWN']:
            script += 'wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/master/ResonanceShapes_%s_13TeV_Spring16_%s.root -P inputs/\n'%(options.model,sys)        
    script += 'python python/RunCombine.py -i %s -m %s --mass %s -c %s --lumi %s -d %s -b %s %s %s --min-tol %e --min-strat %i --rMax %f %s %s %s %s %s\n'%(options.inputFitFile,
                                                                                                                                                         options.model,
                                                                                                                                                         massPoint,
                                                                                                                                                         options.config,
                                                                                                                                                         options.lumi,
                                                                                                                                                         submitDir,
                                                                                                                                                         options.box,
                                                                                                                                                         penaltyString,
                                                                                                                                                         signalSys,
                                                                                                                                                         options.min_tol,
                                                                                                                                                         options.min_strat,
                                                                                                                                                         options.rMax,
                                                                                                                                                         decoString,
                                                                                                                                                         bayesString,
                                                                                                                                                         toyString,
                                                                                                                                                         xsecString,
                                                                                                                                                         signifString)
    script += 'cp %s/higgsCombine* %s/\n'%(submitDir,combineDir)
    script += 'cd ../..\n'
    script += 'rm -rf $TWD\n'
        
    outputfile = open(outputname,'w')
    outputfile.write(script)
    outputfile.close
    
    return outputname,ffDir

def submit_jobs(options,args):    
     
    for massPoint in massIterable(options.mass):

        for iJob in range(0,options.jobs):
            outputname,ffDir = writeBashScript(options,massPoint,iJob)

            pwd = os.environ['PWD']
            os.system("mkdir -p "+pwd+"/"+ffDir)
            os.system("echo bsub -q "+options.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)      
            if not options.dryRun:
                time.sleep(3)
                os.system("bsub -q "+options.queue+" -o "+pwd+"/"+ffDir+"/log.log source "+pwd+"/"+outputname)
    
def main(options,args):
    
    boxes = options.box.split('_')

    signif = options.signif
    
    model = options.model

    lumiFloat = [float(lumiStr) for lumiStr in options.lumi.split('_')]

    rRangeStringList = []
    sysStringList = []
    for box,lumi in zip(boxes,lumiFloat):

        paramDict = {}
        if options.inputFitFile is not None and options.bayes:
            inputRootFile = rt.TFile.Open(options.inputFitFile,"r")
            wIn = inputRootFile.Get("w"+box).Clone("wIn"+box)            
            if wIn.obj("fitresult_extDijetPdf_data_obs") != None:
                frIn = wIn.obj("fitresult_extDijetPdf_data_obs")
            elif wIn.obj("nll_extDijetPdf_data_obs") != None:
                frIn = wIn.obj("nll_extDijetPdf_data_obs")
            elif wIn.obj("fitresult_extDijetPdf_data_obs_with_constr") != None:
                fr = wIn.obj("fitresult_extDijetPdf_data_obs_with_constr")
            elif wIn.obj("nll_extDijetPdf_data_obs_with_constr") != None:
                frIn = wIn.obj("nll_extDijetPdf_data_obs_with_constr")
            elif wIn.obj("simNll") != None:
                frIn = wIn.obj("simNll")
            paramDict = {}
            for p in rootTools.RootIterator.RootIterator(frIn.floatParsFinal()):
                paramDict[p.GetName()] = [p.getVal(), p.getError()]
            print "grabbing parameter ranges +-%gsigma for bayesian"%NSIGMA
    

        signalSys = ''
        if options.noSignalSys or options.noSys:
            signalSys = '--no-signal-sys'
        else:
            if box=='CaloDijet2015' or box=='CaloDijet20152016':
                signalSys  =   '--jesUp inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring15_JESUP.root --jesDown inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring15_JESDOWN.root'%(model,model)
                signalSys += ' --jerUp inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring15_JERUP.root --jerDown inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring15_JERDOWN.root'%(model,model)
            elif box=='CaloDijet2016':
                signalSys  =   '--jesUp inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring16_JESUP.root --jesDown inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring16_JESDOWN.root'%(model,model)
                signalSys += ' --jerUp inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring16_JERUP.root --jerDown inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring16_JERDOWN.root'%(model,model)
            elif box=='PFDijet2016':
                signalSys  =   '--jesUp inputs/ResonanceShapes_%s_13TeV_Spring16_JESUP.root --jesDown inputs/ResonanceShapes_%s_13TeV_Spring16_JESDOWN.root'%(model,model)
                signalSys += ' --jerUp inputs/ResonanceShapes_%s_13TeV_Spring16_JERUP.root'%(model)
        
        penaltyString = ''
        if options.penalty:
            penaltyString = '--penalty'
        elif options.noSys:
            penaltyString = '--fixed'
    
        xsecString = '--xsec %f'%(options.xsec)    

        
        if box=='CaloDijet2015' or box=='CaloDijet20152016':
            signalDsName = 'inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring15.root'%model
        elif box=='CaloDijet2016':
            signalDsName = 'inputs/ResonanceShapes_%s_13TeV_CaloScouting_Spring16.root'%model
        elif box=='PFDijet2016':
            signalDsName = 'inputs/ResonanceShapes_%s_13TeV_Spring16.root'%model
        elif 'PFDijetbb2016' in box:
            signalDsName = 'inputs/ResonanceShapes_%s_13TeV_Spring16.root'%model
        else:
            print "Signal Files: ", options.resFile
            signalDsName = options.resFile
            
        backgroundDsName = {'CaloDijet2015':'inputs/data_CaloScoutingHT_Run2015D_BiasCorrected_CaloDijet2015.root',
                            'CaloDijet2016':'inputs/data_CaloScoutingHT_Run2016BCD_NewBiasCorrectedFlat_Golden12910pb_CaloDijet2016.root',
                            #'PFDijet2016':'inputs/data_PFRECOHT_Run2016BCD_Golden12910pb_PFDijet2016.root',
                            'CaloDijet20152016':'inputs/data_CaloScoutingHT_Run2015D2016B_CaloDijet20152016.root',
                            'PFDijet2016':'inputs/20fb_btag.root',
                            'PFDijetbb20160mt':'inputs/20fb_btag.root',
                            'PFDijetbb20161mt':'inputs/20fb_btag.root',
                            'PFDijetbb20162mt':'inputs/20fb_btag.root',
                            'newVGamma':'inputs/data_CaloScoutingHT_Run2015D_BiasCorrected_CaloDijet2015.root'
                            }

        blindString = ''
        if options.blind:
            blindString = '--noFitAsimov --run expected'

        sysString = ''
        if options.noSys and options.deco:
            sysString = '-S 0 --freezeNuisances=shapeBkg_%s_bkg_deco_%s__norm,deco_%s_eig1,deco_%s_eig2,deco_%s_eig3,jes,jer,lumi'%(box,box,box,box,box)
        elif options.noSys:
            sysString = '-S 0 --freezeNuisances=shapeBkg_%s_bkg_%s__norm,p1_%s,p2_%s,p3_%s,jes,jer,lumi'%(box,box,box,box,box)
        sysStringList.append(sysString)

        decoString = ''
        if options.deco:
            decoString  ='--deco'

        massPoints = options.mass.split("_")
        for massPoint in massPoints:#massIterable(options.mass):
            print "Mass points: ", massPoint, " in ", massPoints
            print model
            print massPoint
            print options.inputFitFile
            print 1000*lumi
            print options.config
            print box
            print options.outDir
            print signalDsName
            print backgroundDsName[box]
            print penaltyString
            print signalSys
            print xsecString
            print decoString
            print "Command: python python/WriteDataCard.py -m %s --mass %s -i %s -l %f -c %s -b %s -d %s %s %s %s %s %s %s"%(model, massPoint, options.inputFitFile,1000*lumi,options.config,box,options.outDir,signalDsName,backgroundDsName[box],penaltyString,signalSys,xsecString,decoString)
            luminosity = ("%.3f"%(lumi)).replace('.','p')
            exec_me('python python/WriteDataCard.py -m %s --mass %s -i %s -l %f -c %s -b %s -d %s %s %s %s %s %s %s'%(model, massPoint, options.inputFitFile,1000*lumi,options.config,box,options.outDir,signalDsName,backgroundDsName[box],penaltyString,signalSys,xsecString,decoString),options.dryRun)    
            if options.bayes:
                rRangeString =  '--setPhysicsModelParameterRanges '
                if options.deco:
                    rRangeString += 'shapeBkg_%s_bkg_deco_%s__norm=%f,%f'%(box,box,1-NSIGMA*paramDict['Ntot_bkg_%s'%box][1]/paramDict['Ntot_bkg_%s'%box][0],1+NSIGMA*paramDict['Ntot_bkg_%s'%box][1]/paramDict['Ntot_bkg_%s'%box][0])
                    rRangeString += ':deco_%s_eig1=%f,%f'%(box,-1.0*NSIGMA,NSIGMA)
                    rRangeString += ':deco_%s_eig2=%f,%f'%(box,-1.0*NSIGMA,NSIGMA)
                    rRangeString += ':deco_%s_eig3=%f,%f'%(box,-1.0*NSIGMA,NSIGMA)
                else:
                    rRangeString += 'shapeBkg_%s_bkg_%s__norm=%f,%f'%(box,box,1-NSIGMA*paramDict['Ntot_bkg_%s'%box][1]/paramDict['Ntot_bkg_%s'%box][0],1+NSIGMA*paramDict['Ntot_bkg_%s'%box][1]/paramDict['Ntot_bkg_%s'%box][0])
                    rRangeString += ':p1_%s=%f,%f'%(box,paramDict['p1_%s'%box][0]-NSIGMA*paramDict['p1_%s'%box][1],paramDict['p1_%s'%box][0]+NSIGMA*paramDict['p1_%s'%box][1])
                    rRangeString += ':p2_%s=%f,%f'%(box,paramDict['p2_%s'%box][0]-NSIGMA*paramDict['p2_%s'%box][1],paramDict['p2_%s'%box][0]+NSIGMA*paramDict['p2_%s'%box][1])
                    rRangeString += ':p3_%s=%f,%f'%(box,paramDict['p3_%s'%box][0]-NSIGMA*paramDict['p3_%s'%box][1],paramDict['p3_%s'%box][0]+NSIGMA*paramDict['p3_%s'%box][1])            
                if options.rMax>-1:
                    rRangeString += ':r=0,%f'%(options.rMax)
                rRangeStringList.append(rRangeString)
                toyString = ''
                if options.toys>-1:
                    toyString = '-t %i -s -1'%options.toys
                if len(boxes)==1:
                    print "Before combine line 279 ", luminosity
                    exec_me('combine -M MarkovChainMC -H Asymptotic %s/dijet_combine_%s_%s_lumi-%%s.txt -n %s_%s_lumi-%s_%s --tries 20 --proposal ortho --burnInSteps 200 --iteration 30000 --propHelperWidthRangeDivisor 10 %s %s %s %s'%(options.outDir,model,massPoint,luminosity,box,model,massPoint,luminosity,box,rRangeString,blindString,sysString,toyString),options.dryRun)
                    exec_me('mv higgsCombine%s_%s_lumi-%.3f_%s.MarkovChainMC.mH120*root %s/'%(model,massPoint,lumi,box,options.outDir),options.dryRun)  
            else:
                if signif:
                    rRangeString = ''               
                    if options.rMax>-1:
                        rRangeString = '--setPhysicsModelParameterRanges r=0,%f'%(options.rMax)
                        rRangeStringList.append(rRangeString)
                    if len(boxes)==1:
                        print "Before combine line 288"
                        exec_me('combine -M ProfileLikelihood -m %s --signif %s/dijet_combine_%s_%s_lumi-%s_%s.txt -n %s_%s_lumi-%s_%s %s %s'%(massPoint,options.outDir,model,massPoint,luminosity,box,model,massPoint,luminosity,box,rRangeString,sysString),options.dryRun)
                        exec_me('mv higgsCombine%s_%s_lumi-%s_%s.ProfileLikelihood.mH120.root %s/'%(model,massPoint,luminosity,box,options.outDir),options.dryRun)
                else:
                    rRangeString = ''
                    if options.rMax>-1:                
                        rRangeString =  '--setPhysicsModelParameterRanges r=0,%f'%(options.rMax)                        
                        rRangeStringList.append(rRangeString)
                    if len(boxes)==1:
                        print "Before combine line 296 ", ("%.3f"%(lumi)).replace('.','p'), " ", lumi
                        exec_me('combine -M Asymptotic -m %s -H ProfileLikelihood %s/dijet_combine_%s_%s_lumi-%s_%s.txt -n %s_%s_lumi-%s_%s --minimizerTolerance %f --minimizerStrategy %i %s --saveWorkspace %s %s'%(massPoint,options.outDir,model,massPoint,luminosity,box,model,massPoint,luminosity,box,options.min_tol,options.min_strat,rRangeString,blindString,sysString),options.dryRun)
                        exec_me('mv higgsCombine%s_%s_lumi-%s_%s.Asymptotic.mH%s.root %s/'%(model,massPoint,luminosity,box,massPoint,options.outDir),options.dryRun)
    if len(boxes)>1:
        lumiTotal = sum(lumiFloat)
        for box,lumi in zip(boxes,lumiFloat): 
          luminosity = ("%.3f"%(lumi/1000.)).replace('.','p')
          exec_me('cp %s/dijet_combine_%s_%s_lumi-%s_%s.txt .'%(options.outDir,model,massPoint,luminosity,box),options.dryRun)        
        cmds = ['%s=dijet_combine_%s_%s_lumi-%s_%s.txt'%(box,model,massPoint,("%.3f"%(lumi/1000.)).replace('.','p'),box) for box,lumi in zip(boxes,lumiFloat)] 
        exec_me('combineCards.py %s > %s/dijet_combine_%s_%s_lumi-%s_%s.txt'%(' '.join(cmds),options.outDir,model,massPoint,("%.3f"%(lumiTotal)).replace('.','p'),options.box),options.dryRun)
        exec_me('cat %s/dijet_combine_%s_%s_lumi-%s_%s.txt'%(options.outDir,model,massPoint,("%.3f"%(lumiTotal)).replace('.','p'),options.box),options.dryRun)
        if options.bayes:
            rRangeStringListMod = [rRangeString.replace('--setPhysicsModelParameterRanges ','') for rRangeString in rRangeStringList ]
            paramRangeList = []
            for listMod in rRangeStringListMod:
                paramRangeList.extend(listMod.split(':'))
            paramRangeList = list(set(paramRangeList))
            rRangeStringTotal = ''
            if options.deco or rMax>=-1:
                rRangeStringTotal = '--setPhysicsModelParameterRanges ' + ','.join(paramRangeList)
                
            sysStringListMod = [sysString.replace('-S 0 --freezeNuisances=','') for sysString in sysStringList ]
            paramFreezeList = []
            for listMod in sysStringListMod:
                paramFreezeList.extend(listMod.split(','))
            paramFreezeList = list(set(paramFreezeList))
            sysStringTotal = ''
            if options.noSys:
                sysStringTotal = '-S 0 --freezeNuisances=' + ','.join(paramFreezeList)
            exec_me('combine -M MarkovChainMC -H Asymptotic %s/dijet_combine_%s_%s_lumi-%s_%s.txt -n %s_%s_lumi-%s_%s --tries 30 --proposal ortho --burnInSteps 1000 --iteration 40000 --propHelperWidthRangeDivisor 10 %s %s %s %s'%(options.outDir,model,massPoint,("%.3f"%(lumiTotal)).replace('.','p'),options.box,model,massPoint,("%.3f"%(lumiTotal)).replace('.','p'),options.box,rRangeStringTotal,blindString,sysStringTotal,toyString),options.dryRun)
            exec_me('mv higgsCombine%s_%s_lumi-%s_%s.MarkovChainMC.mH120*root %s/'%(model,massPoint,("%.3f"%(lumiTotal)).replace('.','p'),options.box,options.outDir),options.dryRun)
        else:
            if signif:
                rRangeString = ''
                if options.rMax>-1:
                    rRangeString = '--setPhysicsModelParameterRanges r=0,%f'%(options.rMax)
                print "Before combine line 331 ", luminosity
                exec_me('combine -M ProfileLikelihood --signif %s/dijet_combine_%s_%s_lumi-%s_%s.txt -n %s_%s_lumi-%s_%s %s %s'%(options.outDir,model,massPoint,("%.3f"%(lumiTotal)).replace('.','p'),options.box,model,massPoint,("%.3f"%(lumiTotal)).replace('.','p'),options.box,rRangeString,sysString),options.dryRun)
                exec_me('mv higgsCombine%s_%s_lumi-%s_%s.ProfileLikelihood.mH120.root %s/'%(model,massPoint,("%.3f"%(lumiTotal)).replace('.','p'),options.box,options.outDir),options.dryRun)
            else:
                rRangeString = ''
                if options.rMax>-1:
                    rRangeString =  '--setPhysicsModelParameterRanges r=0,%f'%(options.rMax)
                print "Before combine line 337 ", luminosity
                exec_me('combine -M Asymptotic -H ProfileLikelihood %s/dijet_combine_%s_%s_lumi-%s_%s.txt -n %s_%s_lumi-%s_%s --minimizerTolerance %f --minimizerStrategy %i %s --saveWorkspace %s %s'%(options.outDir,model,massPoint,("%.3f"%(lumiTotal)).replace('.','p'),options.box,model,massPoint,("%.3f"%(lumiTotal)).replace('.','p'),options.box,options.min_tol,options.min_strat,rRangeString,blindString,sysString),options.dryRun)
                exec_me('mv higgsCombine%s_%s_lumi-%s_%s.Asymptotic.mH120.root %s/'%(model,massPoint,("%.3f"%(lumiTotal)).replace('.','p'),options.box,options.outDir),options.dryRun)
            for box,lumi in zip(boxes,lumiFloat): exec_me('rm dijet_combine_%s_%s_lumi-%s_%s.txt'%(model,massPoint,("%.3f"%(lumi/1000.)).replace('.','p'),box),options.dryRun)
    
if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-b','--box',dest="box", default="CaloDijet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="gg",type="string",
                  help="signal model name")
    parser.add_option('--mass',dest="mass", default='750',type="string",
                  help="mass of resonance")
    parser.add_option('-l','--lumi',dest="lumi", default="1.918",type="string",
                  help="lumi in fb^-1, possibly for different channels e.g.: 1.918_2.590")
    parser.add_option('--signif',dest="signif",default=False,action='store_true',
                  help="calculate significance instead of limit")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
    parser.add_option('--min-tol',dest="min_tol",default=0.001,type="float",
                  help="minimizer tolerance (default = 0.001)")
    parser.add_option('--min-strat',dest="min_strat",default=2,type="int",
                  help="minimizer strategy (default = 2)")
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                  help="Just print out commands to run")
    parser.add_option('--penalty',dest="penalty",default=False,action='store_true',
                  help="penalty terms on background parameters")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default='FitResults/BinnedFitResults.root',type="string",
                  help="input fit file")
    parser.add_option('--no-signal-sys',dest="noSignalSys",default=False,action='store_true',
                  help="do not create signal shape systematic histograms / uncertainties")
    parser.add_option('--no-sys',dest="noSys",default=False,action='store_true',
                  help="no systematic uncertainties when running combine")
    parser.add_option('--blind',dest="blind",default=False,action='store_true',
                  help="run only blinded expected limits")
    parser.add_option('--rMax',dest="rMax",default=-1,type="float",
                  help="maximum r value (for better precision)")
    parser.add_option('--xsec',dest="xsec",default=1,type="float",
                  help="xsec for signal in pb (r = 1)")
    parser.add_option('-j','--jobs',dest="jobs",default=0,type="int",
                  help="number of jobs to submit when running toys for each mass point (just set to 1 for observed limits)")
    parser.add_option('--bayes',dest="bayes",default=False,action='store_true',
                  help="bayesian limits")
    parser.add_option('--deco',dest="deco",default=False,action='store_true',
                  help="decorrelate shape parameters")
    parser.add_option('--tag',dest="tag", default='master',type="string",
                  help="tag for repository")
    parser.add_option('-q','--queue',dest="queue",default="1nh",type="string",
                  help="queue: 1nh, 8nh, 1nd, etc.")
    parser.add_option('-t','--toys',dest="toys",default=-1,type="int",
                  help="number of toys per job(for bayesian expected limits)")
    parser.add_option('--resFile', dest="resFile", default="file.root", type="string",
                  help="file containing the resonances distributions")


    (options,args) = parser.parse_args()


    if options.jobs:
        submit_jobs(options,args)
    else:            
        main(options,args)
