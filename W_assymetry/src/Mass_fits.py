from sys import dont_write_bytecode
import ROOT 
from ROOT import TTree, TFile, Double, TF1, TMath, TStyle
import math as m
import numpy as np

__print = True
__avarage = True

vars = ["M_W_p", "M_W_m"]

def plot_fit(var):
	
	c = ROOT.TCanvas("c")
	data = ROOT.TFile("data_histos.root")
	sig = ROOT.TFile("wjets_histos.root")
	data_hist_pm = data.Get(var)
	sig_hist_pm = sig.Get(var)
	bkg = ROOT.TFile("dy_histos.root")
	bkg_hist_pm = bkg.Get(var)
	if var == "M_W_p":
		data_hist_pm.GetXaxis().SetTitle("M_{T,W^{+}} (GeV)")
		data_hist_pm.SetTitle("Transverse mass distribution of the #it{W}^{+}-boson")
		bkg1 = ROOT.TFile("qcd_histos.root")
		

  
	elif var == "M_W_m":
		data_hist_pm.GetXaxis().SetTitle("M_{T,W^{-}} (GeV)")
		data_hist_pm.SetTitle("Transverse mass distribution of the #it{W}^{-}-boson")
	fit_sig = ROOT.TF1("fit_sig","gaus(0)+gaus(3)",0,150) 
	fit_sig.SetParameters(500,100,10,300,50,40)
	fit_sig.SetLineColor(215)
	fit_sig.SetLineStyle(2)
	fit_bkg = ROOT.TF1("fit_bkg1","[0]*(TMath::Erf((x-[1])/[2])+1.)",0,150)
	fit_bkg.SetParameters(10,60,-10)
	fit_bkg.SetLineColor(210)
	fit_bkg.SetLineStyle(2)
	fit_bkg_sig = ROOT.TF1("fit_bkg_sig","[0]*(TMath::Erf((x-[1])/[2])+1)+gaus(3)+gaus(6)",0,150)
	fit_bkg_sig.SetLineWidth(2)
	sig_hist_pm.Fit("fit_sig")
	bkg_hist_pm.Fit("fit_bkg")
	for i in range(9):
		if i==0:
			fit_bkg_sig.SetParameter(i,fit_bkg.GetParameter(i))
		elif i<3: ## we fix background parameters, but the normalization
			fit_bkg_sig.FixParameter(i,fit_bkg.GetParameter(i))
		else:
			fit_bkg_sig.SetParameter(i,fit_sig.GetParameter(i-3))
   
	data_hist_pm.Fit(fit_bkg_sig)
	data_hist_pm.GetYaxis().SetTitle("Number of Events")
	data_hist_pm.Draw()
	fit_bkg_sig.Draw("Same")
	fit_sig.Draw("Same")
	fit_bkg.Draw("Same")


	leg = ROOT.TLegend(0.1, 0.6, 0.4, 0.89)
	leg.SetFillColor(0)
	leg.AddEntry(fit_sig, "fitted signal ")
	leg.AddEntry(fit_bkg, "fitted background ")
	leg.AddEntry(fit_bkg_sig, "fitting function #it{f}_{fit} ")
	leg.Draw("Same")

	c.SaveAs(var + "_fitted.pdf")

	I = fit_bkg_sig.Integral(0., 150.)
	D_I = fit_bkg_sig.IntegralError(0.,150)
	return fit_bkg_sig, I, D_I

"""
def make_a(fit, para = True):
	values = np.zeros(9) 
	if para:
		for i in range(len(values)):
			values[i] = fit.GetParameter(i) 
	else:
		for i in range(len(values)):
			values[i] = fit.GetParError(i)

	return values


values = np.array([np.zeros(9), np.zeros(9)])
values[0] = make_a(plot_fit(vars[0])[0], True)
values[1] = make_a(plot_fit(vars[1])[0], True)
"""
N_W_p = plot_fit(vars[0])[1]
N_W_m = plot_fit(vars[1])[1]

D_N_W_p = plot_fit(vars[0])[2]
D_N_W_m = plot_fit(vars[1])[2]


def asymmetry(N_W_p, N_W_m,D_N_W_p, D_N_W_m):
	A = (N_W_p - N_W_m)/(N_W_p + N_W_m)
	D_A = m.sqrt((2*N_W_m*D_N_W_p/(N_W_p+N_W_m)**2)**2 + (2*N_W_p*D_N_W_m/(N_W_p+N_W_m)**2)**2)

	return A, D_A

if __print: 
    print("#################################")
    print("########## Asymmetry ############")
    print("A = %s +/- %s"%(asymmetry(N_W_p, N_W_m,D_N_W_p, D_N_W_m )[0], asymmetry(N_W_p, N_W_m, D_N_W_p, D_N_W_m)[1]))
    print("#################################")
    print("#################################")
    print("####### Number of signal ########")
    print("N^+ = %s +/- %s"%(N_W_p, D_N_W_p))
    print("N^- = %s +/- %s"%(N_W_m,D_N_W_m))
    print("#################################")


if __avarage:
    def avar_A(A1, A2, D_A1, D_A2):
        n = A1/D_A1**2 + A2/D_A2**2
        d = 1./D_A1**2 + 1/D_A2**2
        err = m.sqrt(1./d)
        return n/d, err
    
    
    print("#################################")
    print("########## Asymmetry ############")
    print("A = %s +/- %s"%(avar_A(0.11504138219327226, 0.1301225666219243, 0.0057735063533351175, 0.013297464734210798 )))