from MyAnalysis_1 import MyAnalysis
from ROOT import TTree, TFile, Double
from Plotter_1 import plotVar, plotVarNorm, plotShapes, getHisto, getBkgHisto, getSigHisto
import math as  m
import numpy as np
### Instantiation of an object of kind MyAnalysis for each single sample
__myanalysis = False

if __myanalysis:  
  TT = MyAnalysis("ttbar")
  TT.processEvents()

  DY = MyAnalysis("dy")
  DY.processEvents()

  QCD = MyAnalysis("qcd")
  QCD.processEvents()

  SingleTop = MyAnalysis("single_top")
  SingleTop.processEvents()

  WJets = MyAnalysis("wjets")
  WJets.processEvents()

  WW = MyAnalysis("ww")
  WW.processEvents()

  ZZ = MyAnalysis("zz")
  ZZ.processEvents()

  WZ = MyAnalysis("wz")
  WZ.processEvents()

  Data = MyAnalysis("data")
  Data.processEvents()




samples = ["qcd", "zz", "wz", "ww", "single_top", "dy","wjets", "ttbar"]

vars = ["NIsoMu", "Muon_Pt","Muon_Iso","Jet_Pt","NJet","NBtag","M_muon_blept"]



for v in vars:
  

  print ("Variable: ", v)
  if v in vars:
    plotVar(v, samples, True, False)  
    ### plotShapes (variable, samples,logScale )
    #plotShapes(v, samples,  True)
    ### plotVar(variable, samples,isData, logScale )
    
      
  else:
    plotVar(v, samples, True, True)

def getIntegral(histo):
  minBin = histo.GetXaxis().GetFirst()
  maxBin = histo.GetXaxis().GetLast()
  nDataErr = Double()
  nData = histo.IntegralAndError(minBin, maxBin, nDataErr)
  ### nDataErr will be overwritten with the correct value for the uncertainty
  #print ("Integral %f +/- %f"%(nData,nDataErr))
  return(nData,nDataErr)


def SB_Ratio(name):
  
  samples = ["qcd", "zz", "wz", "ww", "single_top", "dy","wjets"] ## ttbar in getSighisto,..
  ### simulated events 
  S = getIntegral(getSigHisto(name))[0] 
  D_S = getIntegral(getSigHisto(name))[1]
  
  B = getIntegral(getBkgHisto(name,samples))[0] 
  D_B = getIntegral(getBkgHisto(name,samples))[1]

  ### measured events 
  D = getIntegral(getHisto(name,"data"))[0]  
  D_D = getIntegral(getHisto(name,"data"))[1]

  return (S/m.sqrt(S+B), S, D_S, B, D_B, D, D_D)
    

def cross_section(name):  
  
  S = SB_Ratio(name)[1]
  D_S = SB_Ratio(name)[2]
  
  B = SB_Ratio(name)[3]
  D_B = SB_Ratio(name)[4]
  
  D = SB_Ratio(name)[5]
  D_D = SB_Ratio(name)[6]  
  
  #### integrated luminosity ####
  
  L = 50.
  D_L =L/10.
  
  #### acceptance ####
  A = S/7928.61
  D_A = m.sqrt((D_S/7928.61)**2 +  + (S / 7928.61**(3/2))**2 )
    
  #### number of ttbar events ####
  N_ev = (D-B) 
  D_N_ev = m.sqrt(D_B**2 + D_D**2)
  
  #### epsilon trigger ####
  trigger_eff = S/410.5962505340576
  D_S_No_trigg = 10.408828128763657
  D_trigger_eff = m.sqrt((D_S/410.5962505340576)**2 + (S*D_S_No_trigg/410.5962505340576**2)**2)
  
  
      
    
  #### cross section ####
  cross_section = N_ev/(A *trigger_eff*L)
  D_cross_section = m.sqrt((D_N_ev/(A*trigger_eff*L))**2 +(N_ev*D_A/(A**2*trigger_eff*L))**2 + (N_ev*D_trigger_eff/(A*trigger_eff**2*L))**2 + (N_ev*D_L/(A*trigger_eff*L**2))**2)

  
  return (cross_section, D_cross_section, trigger_eff, D_trigger_eff, N_ev, D_N_ev, A, D_A)

def chi_squared(M, D_M, T, D_T):
  chi2 = 1./2.*(M-T)**2/(D_M**2+ D_T**2)
  return chi2

order = ["Muon_Pt","NIsoMu","Jet_Pt","NJet","NBtag","M_muon_blept"]

SB_R = []
for o in order:
  SB_R.append(SB_Ratio(o)[0])

cross_sect_th = 173.60
D_cross_sect_th = 11.24

winner = order[SB_R.index(max(SB_R))]
print("########### %s ###########"%(winner))
print("Signal %s = %s +/- %s"%(winner,SB_Ratio(winner)[1], SB_Ratio(winner)[2]))
print("Data %s = %s +/- %s"%(winner,SB_Ratio(winner)[5], SB_Ratio(winner)[6]))
print("Background %s = %s +/- %s"%(winner,SB_Ratio(winner)[3], SB_Ratio(winner)[4]))
print("S/(sqrt(S+B)) = %s "%(SB_Ratio(winner)[0]))
print("Triggere efficency = %s +/- %s"%(cross_section(winner)[2], cross_section(winner)[3]))
print("Number measured events = %s +/- %s"%(cross_section(winner)[4], cross_section(winner)[5]))
print("Acceptance = %s +/- %s"%(cross_section(winner)[6], cross_section(winner)[7]))
print("Cross section = %s +/- %s"%(cross_section(winner)[0], cross_section(winner)[1]))
print("Chi square = %s"%chi_squared(cross_section(winner)[0], cross_section(winner)[1],cross_sect_th, D_cross_sect_th))
