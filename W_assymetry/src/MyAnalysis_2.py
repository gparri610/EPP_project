from re import S
import ROOT 
import copy
from Samples_2 import samp
import math as m
import numpy as np


def D_R(p1, p2):
    D_eta = p1.PseudoRapidity() - p2.PseudoRapidity()
    D_phi = p1.Phi()-p2.Phi()
    return m.sqrt(D_eta**2 + D_phi**2)


def W_p_hist(self,title):

    name = ROOT.TH1F(title,"W^+ mass", 100, 0., 150.)
    name.SetXTitle("W^+ mass")
    self.histograms[title] = name
    
def W_m_hist(self,title):
    name = ROOT.TH1F(title,"W^- mass", 100, 0., 150.)
    name.SetXTitle("W^- mass")
    self.histograms[title] = name

def W_p_MET_hist(self,title):
    name = ROOT.TH1F(title,"MET W^+", 100, 0., 150.)
    self.histograms[title] = name

def W_m_MET_hist(self,title):
    name = ROOT.TH1F(title,"MET W^-", 100, 0., 150.)
    self.histograms[title] = name

class MyAnalysis(object):
   
    def __init__(self, sample):

        """ The Init() function is called when an object MyAnalysis is initialised
        The tree corresponding to the specific sample is picked up 
        and histograms are booked.
        """

        self._tree = ROOT.TTree()        
        if(sample not in samp.keys() and sample != "data"):
            print (RuntimeError("Sample %s not valid. please, choose among these: %s" % (sample, str(samp.keys())) ) )
            exit
        self.histograms = {}
        self.sample = sample
        self._file = ROOT.TFile("/home/uzh/gparri/EPP_ex/project/EXERCISE/files/"+sample+".root")
        self._file.cd()
        tree = self._file.Get("events")
        self._tree = tree
        self.nEvents = self._tree.GetEntries()
        print ("Number of entries for " + self.sample + ": " + str(self.nEvents))
        
        ### Book histograms
        self.bookHistos()

    def getTree(self):
        return self._tree

    def getHistos(self):
        return self.histograms

    def bookHistos(self):
        
        W_m_hist(self,"M_W_m")
        W_p_hist(self,"M_W_p")
        W_p_MET_hist(self,"W_p_MET")
        W_m_MET_hist(self,"W_m_MET")
        
    
    
    

    
    
    def saveHistos(self):
        outfilename = self.sample + "_histos.root"
        outfile = ROOT.TFile(outfilename, "RECREATE")
        outfile.cd()
        for h in self.histograms.values():
            h.Write()
        outfile.Close()

    

    ### processEvent function implements the actions to perform on each event
    ### This is the place where to implement the analysis strategy: study of most sensitive variables
    ### and signal-like event selection
    
    
    def processEvent(self, entry):
        tree = self.getTree()
        tree.GetEntry(entry)
        w = tree.EventWeight
        
        
        if tree.triggerIsoMu24:
        
            min_muon_Pt = 30.
            muon_I_rel = 0.12
            max_eta =0.4
            
            if tree.NMuon ==1:
                muon = ROOT.TLorentzVector(tree.Muon_Px[0],tree.Muon_Py[0],tree.Muon_Pz[0],tree.Muon_E[0])
                MET = m.sqrt(tree.MET_px**2 + tree.MET_py**2)
                ### W- histogram
                
                if tree.Muon_Charge[0] < 0:
                    if (muon.Pt() > min_muon_Pt and (tree.Muon_Iso[0]/muon.Pt()) < muon_I_rel):
                        muon_T = ROOT.TLorentzVector(tree.Muon_Px[0],tree.Muon_Py[0], 0, m.sqrt(-tree.Muon_Pz[0]**2+ tree.Muon_E[0]**2))
                        neutrino = ROOT.TLorentzVector(tree.MET_px,tree.MET_py, 0., MET) 
                        cosphi = (muon.Px()*tree.MET_px + muon.Py()*tree.MET_py)/(MET*m.sqrt(muon.Px()**2 + muon.Py()**2))             
                        M_T_m = m.sqrt(2. * muon.Pt() * MET * (1. - cosphi))
                        
                        if abs(muon.PseudoRapidity())< max_eta:
                            self.histograms["M_W_m"].Fill(M_T_m,w)
                            self.histograms["W_m_MET"].Fill(MET, w)
                
                ### W+ historgram
                
                if tree.Muon_Charge[0] > 0:
                    if (muon.Pt() > min_muon_Pt and (tree.Muon_Iso[0]/muon.Pt()) < muon_I_rel):
                        muon_T = ROOT.TLorentzVector(tree.Muon_Px[0],tree.Muon_Py[0], 0, m.sqrt(-tree.Muon_Pz[0]**2+ tree.Muon_E[0]**2))
                        neutrino = ROOT.TLorentzVector(tree.MET_px,tree.MET_py, 0., MET) 
                        cosphi = (muon.Px()*tree.MET_px + muon.Py()*tree.MET_py)/(MET*m.sqrt(muon.Px()**2 + muon.Py()**2))             
             
                        M_T_p = m.sqrt(2. * muon.Pt() * MET * (1. - cosphi ))
                        
                        
                        
                        if abs(muon.PseudoRapidity())< max_eta:
                            self.histograms["M_W_p"].Fill(M_T_p,w)
                            self.histograms["W_p_MET"].Fill(MET, w)
                            
                            
    ### processEvents run the function processEvent on each event stored in the tree    
    def processEvents(self):
        nevts = self.nEvents
        for i in range(nevts):
            self.processEvent(i)
        
        self.saveHistos()                     
                        