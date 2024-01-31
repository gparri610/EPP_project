import ROOT 
import copy
from Samples_1 import samp
import math as m


#no_trigger = True

def D_R(p1, p2):
    D_eta = p1.PseudoRapidity() - p2.PseudoRapidity()
    D_phi = p1.Phi()-p2.Phi()
    return m.sqrt(D_eta**2 + D_phi**2)

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
        h_nJet = ROOT.TH1F("NJet","#of jets", 6, 2.5, 8.5)
        h_nJet.SetXTitle("%# of jets")
        self.histograms["NJet"] = h_nJet 

        h_nJetFinal = ROOT.TH1F("NJetFinal","#of jets", 6, -0.5, 6.5)
        h_nJetFinal.SetXTitle("%# of jets")
        self.histograms["NJetFinal"] = h_nJetFinal 

        h_MuonIso = ROOT.TH1F("Muon_Iso","Muon Isolation", 25, 0., 3.)
        h_MuonIso.SetXTitle("Muon Isolation")
        self.histograms["Muon_Iso"] = h_MuonIso 

        h_NIsoMu = ROOT.TH1F("NIsoMu","Number of isolated muons", 5, 0.5, 5.5)
        h_NIsoMu.SetXTitle("Number of isolated muons")
        self.histograms["NIsoMu"] = h_NIsoMu 

        h_MuonPt = ROOT.TH1F("Muon_Pt","Muon P_T", 50, 0., 200.)
        h_MuonPt.SetXTitle("Muon P_T")
        self.histograms["Muon_Pt"] = h_MuonPt 

        h_METpt = ROOT.TH1F("MET_Pt","MET P_T", 25, 0., 300.)
        h_METpt.SetXTitle("MET P_T")
        self.histograms["MET_Pt"] = h_METpt 

        h_JetPt = ROOT.TH1F("Jet_Pt","Jet P_T", 50, 0., 200.)
        h_JetPt.SetXTitle("Jet P_T")
        self.histograms["Jet_Pt"] = h_JetPt 

        h_JetBtag = ROOT.TH1F("Jet_Btag","Jet B tag", 10, 1., 6.)
        h_JetBtag.SetXTitle("Jet B tag")
        self.histograms["Jet_btag"] = h_JetBtag 

        h_NBtag = ROOT.TH1F("NBtag","Jet B tag", 4, 0.5, 4.5)
        h_NBtag.SetXTitle("Number of B tagged jets")
        self.histograms["NBtag"] = h_NBtag 

        h_M_muon_blept = ROOT.TH1F("M_muon_blept","invariant mass of muon and leptonic b jet", 50, 0., 310.)
        h_M_muon_blept.SetXTitle("Invariant mass of muon + leptonic b-jet system")
        self.histograms["M_muon_blept"] = h_M_muon_blept

    
    
    

    
    
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

        ### Muon selection - Select events with at least 1 isolated muon 
        ### with pt>25 GeV to match trigger requirements
             
        
        
        if (tree.triggerIsoMu24): ### triggerIsoMu24 = True for exp data, for MC true/false (see histogram_MC)
        #if (no_trigger):
            ### initialize the bool values ###
            MET_ok = False
            muon_ok = False
            jet_ok = False

            ###  check if MET passes the requirements (MET > 30.)
            
            min_MET = 20.
            MET_value = (tree.MET_px**2 + tree.MET_py**2)**(0.5)
            if MET_value > min_MET:
                MET_ok = True
            else:
                MET_ok = False
            
            ###  check if muons pass the requirements
            
            if MET_ok:
                ### my strategy parameters #######
                muon_I_rel = 0.12
                max_eta = 2.1
        
                ##################################
                min_muon_Pt = 25.
                ##muonRelIsoCut = 0.05
                N_Iso_muons = 0
                
                if tree.NMuon >0:
                    for m in range(tree.NMuon):
                        muon = ROOT.TLorentzVector(tree.Muon_Px[m],tree.Muon_Py[m],tree.Muon_Pz[m],tree.Muon_E[m])
                        if(muon.Pt()>min_muon_Pt and (tree.Muon_Iso[m]/muon.Pt()) < muon_I_rel and abs(muon.PseudoRapidity()) < max_eta):
                            
                            muon_ok = True
                        else:
                            muon_ok = False
                else:
                    muon_ok = False
                    
                ###  check if jet observables pass the requirements
                
                if muon_ok:
                    ### my strategy parameters #######
                    min_pt_Jet = 40.  
                    max_eta_jet = 2.5
                    ##################################
                    
                    ### A jet is tagged as b-jet if b-tag discriminant is greater than 2.0
                    b_tag = 2.0
                    btagged = 0.0
                    
                    if (tree.Jet_ID):   ###Jet quality identifier - Used to reject jets from detector noise
                        if tree.NJet > 3: ### At least four jets are required with p_T>40 GeV  and |\eta|<2.5
                            for j in range(tree.NJet):
                                jet = ROOT.TLorentzVector(tree.Jet_Px[j],tree.Jet_Py[j],tree.Jet_Pz[j],tree.Jet_E[j])

                                if (jet.Pt() > min_pt_Jet and abs(jet.PseudoRapidity()) < max_eta_jet):
                                    jet_ok = True

                                    if (tree.Jet_btag[j] > b_tag):
                                        btagged += 1

                                else:
                                    jet_ok = False
                            
                            if (btagged == 0.0):
                                jet_ok = False

                        else:
                            jet_ok = False

                    else:
                        jet_ok = False
        
    ### MET, muons and jets pass all the requirements 
                
            N_Iso_muons = 0.
            Iso_muons = []
            N_b_b_tagged = 0.
            btagged_jets = []
            
            if jet_ok:                
                for m in range(tree.NMuon):
                    muon =ROOT.TLorentzVector(tree.Muon_Px[m],tree.Muon_Py[m],tree.Muon_Pz[m],tree.Muon_E[m])
                    self.histograms["Muon_Iso"].Fill(tree.Muon_Iso[m], w)                
                    self.histograms["Muon_Pt"].Fill(muon.Pt(), w)
                    N_Iso_muons += 1
                    Iso_muons.append(muon)
                    
                self.histograms["NIsoMu"].Fill(N_Iso_muons, w)
                    
                for n in range(tree.NJet):
                    jet = ROOT.TLorentzVector(tree.Jet_Px[n],tree.Jet_Py[n],tree.Jet_Pz[n],tree.Jet_E[n])
                    self.histograms["Jet_Pt"].Fill(jet.Pt(),w)
                    self.histograms["NJet"].Fill(tree.NJet,w)
                        
                    if (tree.Jet_btag[n] > b_tag):
                        N_b_b_tagged += 1
                        btagged_jets.append(jet)
                        
                        
                    
                self.histograms["NBtag"].Fill(N_b_b_tagged,w)
            
                    
                ### muon blepton  ###
                    
                if len(Iso_muons) == 1:
                        
                    p_muon = Iso_muons[0]
                    
                    if len(btagged_jets)==1:
                        p_lep_b = btagged_jets[0]
                        Mass_mu_b = (p_lep_b + p_muon).M()
                        self.histograms["M_muon_blept"].Fill(Mass_mu_b,w)                        
                    
                    elif len(btagged_jets) == 2:
                        
                        p_lep_b_1 = btagged_jets[0]
                        p_lep_b_2 = btagged_jets[1]

                        if D_R(p_lep_b_1,p_muon) < D_R(p_lep_b_2,p_muon):
                            p_lep_b = p_lep_b_1
                        else:
                            p_lep_b = p_lep_b_2

                        Mass_mu_b = (p_lep_b + p_muon).M()
                        self.histograms["M_muon_blept"].Fill(Mass_mu_b,w)         
                    
                    
                    elif len(btagged_jets) == 3:
                        
                        p_lep_b_1 = btagged_jets[0]
                        p_lep_b_2 = btagged_jets[1]
                        p_lep_b_3 = btagged_jets[2]

                        if D_R(p_lep_b_1,p_muon) < D_R(p_lep_b_2,p_muon) and D_R(p_lep_b_1,p_muon) < D_R(p_lep_b_3,p_muon) :
                            p_lep_b = p_lep_b_1
                        elif D_R(p_lep_b_2,p_muon) < D_R(p_lep_b_1,p_muon) and D_R(p_lep_b_2,p_muon) < D_R(p_lep_b_3,p_muon) :
                            p_lep_b = p_lep_b_2
                        else:
                            p_lep_b = p_lep_b_3

                        Mass_mu_b = (p_lep_b + p_muon).M()
                        self.histograms["M_muon_blept"].Fill(Mass_mu_b,w)          

    ### processEvents run the function processEvent on each event stored in the tree    
    def processEvents(self):
        nevts = self.nEvents
        for i in range(nevts):
            self.processEvent(i)
        
        self.saveHistos()
