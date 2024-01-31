from MyAnalysis_2 import MyAnalysis
from ROOT import TTree, TFile, Double
from Plotter_2 import plotVar, plotVarNorm, getSigHisto, getBkgHisto, getHisto
import math as m
import numpy as np


### Instantiation of an object of kind MyAnalysis for each single sample
__my_analysis = True
__print = True

if __my_analysis:
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
vars = ["W_p_MET", "W_m_MET"] 
#vars = ["M_W_p", "M_W_m"]

for v in vars:
   print ("Variable: ", v)
   plotVar(v, samples,  True, False)
    
def getIntegral(histo):
    minBin = histo.GetXaxis().GetFirst()
    maxBin = histo.GetXaxis().GetLast()
    nDataErr = Double()
    nData = histo.IntegralAndError(minBin, maxBin, nDataErr)
    ### nDataErr will be overwritten with the correct value for the uncertainty
    #print ("Integral %f +/- %f"%(nData,nDataErr))
    return(nData,nDataErr)



def S_B_D(name, make_print=False, make_return = False):
    
    samples = ["qcd", "zz", "wz", "ww", "single_top", "dy"]
    
    
    S = getIntegral(getHisto(name, "wjets"))[0] 
    D_S = getIntegral(getHisto(name, "wjets"))[1]
  
    B = getIntegral(getBkgHisto(name,samples))[0] 
    D_B = getIntegral(getBkgHisto(name,samples))[1]

    ### measured events 
    D = getIntegral(getHisto(name,"data"))[0]  
    D_D = getIntegral(getHisto(name,"data"))[1]

    if make_print:
        print("############################")
        print("#####%s################"% name)
        print("############################")
        print("Signal  = %s +/- %s"%(S, D_S))
        print("Data = %s +/- %s"%(D,D_D))
        print("Background  = %s +/- %s"%(B, D_B))
        print("D-B = %s +/- %s"%(D-B, m.sqrt(D_D**2 + D_B**2)))

    if make_return:
        
        return (D-B,m.sqrt(D_D**2 + D_B**2)) 
    


def asymmetry(name): # name is an array

    N_W_p_m = [] # [N_W^+, N_W^-]
    D_N_W_p_m = []
    
    for v in name:
        N_W_p_m.append(S_B_D(v, False, True)[0])
        D_N_W_p_m.append(S_B_D(v, False, True)[1])

    A = (N_W_p_m[0]-N_W_p_m[1])/(N_W_p_m[0]+N_W_p_m[1])
    D_A = m.sqrt( (2*N_W_p_m[1]*D_N_W_p_m[0]/(N_W_p_m[0]+N_W_p_m[1])**2)**2 + (2*N_W_p_m[0]*D_N_W_p_m[1]/(N_W_p_m[0]+N_W_p_m[1])**2)**2)
    
    return (A,D_A, N_W_p_m , D_N_W_p_m)


if __print: 
    print("#################################")
    print("########## Asymmetry ############")
    print("A = %s +/- %s"%(asymmetry(vars)[0], asymmetry(vars)[1]))
    print("#################################")
    print("#################################")
    print("####### Number of signal ########")
    print("N^+ = %s +/- %s"%(asymmetry(vars)[2][0], asymmetry(vars)[3][0]))
    print("N^- = %s +/- %s"%(asymmetry(vars)[2][1], asymmetry(vars)[3][1]))
    print("#################################")

    



        
    
    

    

