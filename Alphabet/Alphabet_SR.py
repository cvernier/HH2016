import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import scipy

# Our functions:
import Alphabet_Header
from Alphabet_Header import *
import Plotting_Header
from Plotting_Header import *
import Converters
from Converters import *
import Distribution_Header
from Distribution_Header import *
import Alphabet
from Alphabet import *
import CMS_lumi, tdrstyle


CMS_lumi.lumi_13TeV = "2.7 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12

iPeriod =4


presel = "(triggerpassbb>0.&dijetmass_corr>800&jet1tau21<0.6&jet2tau21<0.6&jet2pmass<135&jet2pmass>95&jet2bbtag>0.6 && jet1ID>0.&&jet2ID>0)"
tag = "(triggerpassbb>0.&dijetmass_corr>800&jet1tau21<0.6&jet2tau21<0.6&jet2pmass<135&jet2pmass>105&&(jet1pmass<135&jet1pmass>105)&jet2bbtag>0.6 &jet1bbtag>0.6 && jet1ID>0.&&jet2ID>0)" 
antitag = "(triggerpassbb>0.&dijetmass_corr>800&jet1tau21<0.6&jet2tau21<0.6&jet2pmass<135&jet2pmass>105&&(jet1pmass<135&jet1pmass>105)&jet2bbtag>0.6 &jet1bbtag<0.6&& jet1ID>0.&&jet2ID>0.)"

#jet1ID>0.&&jet2ID>0



DATA = DIST("DATA", "/uscms_data/d3/asady1/HHcode/CMSSW_7_4_2/src/Analysis/HbbAnalysis/Jet_HT_76X.root","myTree","1.")
DistsWeWantToEstiamte = [DATA]
Hbb = Alphabetizer("Analysis_Unblinding", DistsWeWantToEstiamte, [])
var_array = ["jet1pmass", "jet1bbtag", 30,50,200, 200, -1.1, 1.1]
Hbb.SetRegions(var_array, presel) # make the 2D plot
Hbb.TwoDPlot.SetStats(0)
C1 = TCanvas("C1", "", 800, 600)
C1.cd()
Hbb.TwoDPlot.Draw("COLZ")
Hbb.TwoDPlot.GetXaxis().SetTitle("jet mass (GeV)")
Hbb.TwoDPlot.GetYaxis().SetTitle("bb-tag")
C1.SaveAs("CutPlane.root")

cut = [0.6, ">"]
bins = [[50,75],[75,105],[135,160],[160,200]]
truthbins = [[105,135]]
center = 120.
F = QuadraticFit([0.1,0.1,0.1], -75, 75, "quadfit", "EMRFNEX0")
Hbb.GetRates(cut, bins, truthbins, center, F)
leg = TLegend(0.4,0.65,0.6,0.89)
leg.SetLineColor(0)
leg.SetFillColor(4001)
leg.SetTextSize(0.03)
leg.AddEntry(Hbb.G, "events used in fit", "PLE")
leg.AddEntry(Hbb.Fit.fit, "fit", "L")
leg.AddEntry(Hbb.Fit.ErrUp, "fit errors", "L")
C2 = TCanvas("C2", "", 800, 800)
C2.cd()
Hbb.G.SetTitle("")
Hbb.G.Draw("AP")
Hbb.G.GetXaxis().SetTitle("m_{J} - m_{H} (GeV)")
Hbb.G.GetYaxis().SetTitle("R_{p/f}")
Hbb.G.GetYaxis().SetTitleOffset(1.3)
#Hbb.truthG.Draw("P same")
Hbb.Fit.fit.Draw("same")
Hbb.Fit.ErrUp.SetLineStyle(2)
Hbb.Fit.ErrUp.Draw("same")
Hbb.Fit.ErrDn.SetLineStyle(2)
Hbb.Fit.ErrDn.Draw("same")
leg.Draw()
CMS_lumi.CMS_lumi(C2, iPeriod, iPos)
C2.SaveAs("Fit.pdf")

variable = "dijetmass_corr"

binBoundaries = [800, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687,
	    1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019 ]


Hbb.MakeEstVariable(variable, binBoundaries, antitag, tag)
FILE = TFile("Hbb_Unblind.root", "RECREATE")
FILE.cd()
V = TH1F("data_obs", "", len(binBoundaries)-1, array('d',binBoundaries))
for i in Hbb.hists_MSR:
	V.Add(i,1.)
# the estimate is the sum of the histograms in self.hists_EST and self.hist_MSR_SUB
N = TH1F("EST", "", len(binBoundaries)-1, array('d',binBoundaries))
for i in Hbb.hists_EST:
	N.Add(i,1.)
# We can do the same thing for the Up and Down shapes:
NU = TH1F("EST_CMS_scale_13TeVUp", "", len(binBoundaries)-1, array('d',binBoundaries)) 
for i in Hbb.hists_EST_UP:
	NU.Add(i,1.)
ND = TH1F("EST_CMS_scale_13TeVDown", "", len(binBoundaries)-1, array('d',binBoundaries)) 
for i in Hbb.hists_EST_DN:
	ND.Add(i,1.)
A =  TH1F("EST_Antitag", "", len(binBoundaries)-1, array('d',binBoundaries)) 
for i in Hbb.hists_ATAG:
	    A.Add(i,1.)

for bin in range(0,len(binBoundaries)-1):
	if not A.GetBinContent(bin+1) > 0.:
		print A.GetBinError(bin+1)
		A.SetBinError(bin+1, 2.0)
		A.SetBinContent(bin+1, 0.001)
		N.SetBinContent(bin+1,0.0001)
		ND.SetBinContent(bin+1,0.00001)
		NU.SetBinContent(bin+1,0.001)

FILE.Write()
FILE.Save()

vartitle = "M_{jj}^{red} (GeV)"

NU.SetLineColor(kBlack)
ND.SetLineColor(kBlack)
NU.SetLineStyle(2)
ND.SetLineStyle(2)
N.SetLineColor(kPink+3)
N.SetLineWidth(3)



#V.SetStats(0)
V.Sumw2()
V.SetLineColor(1)
V.SetFillColor(0)
V.SetMarkerColor(1)
V.SetMarkerStyle(20)
N.GetYaxis().SetTitle("Events / Bin Width (GeV)")#/ "+str((hbins[2]-hbins[1])/hbins[0])+" GeV")
N.GetXaxis().SetTitle(vartitle)
N.GetXaxis().SetTitleSize(0.04)
N.GetXaxis().SetTitleOffset(1.1)

FindAndSetMax([V,N, NU, ND])
Pull = V.Clone("Pull")
Pull.Add(N, -1.)

Boxes = []
sBoxes = []
pBoxes = []
maxy = 0.
for i in range(1, N.GetNbinsX()+1):
	P = Pull.GetBinContent(i)
	Ve = V.GetBinError(i)
	if Ve > 1.:
		Pull.SetBinContent(i, P/Ve)
	Pull.SetBinError(i, 1.)
	a = A.GetBinError(i)*N.GetBinContent(i)/A.GetBinContent(i)
	u = NU.GetBinContent(i) - N.GetBinContent(i)
	d = N.GetBinContent(i) - ND.GetBinContent(i)
	x1 = Pull.GetBinCenter(i) - (0.5*Pull.GetBinWidth(i))
	y1 = N.GetBinContent(i) - math.sqrt((d*d) + (a*a))
	s1 = N.GetBinContent(i) - a
	if y1 < 0.:
		y1 = 0
	if s1 < 0:
		s1 = 0
	x2 = Pull.GetBinCenter(i) + (0.5*Pull.GetBinWidth(i))
	y2 = N.GetBinContent(i) + math.sqrt((u*u) + (a*a))
	s2 = N.GetBinContent(i) + a
	if maxy < y2:
		maxy = y2
	if Ve > 1.:
		yP1 = -math.sqrt((d*d) + (a*a))/Ve
		yP2 = math.sqrt((u*u) + (a*a))/Ve
	else:
		yP1 = -math.sqrt((d*d) + (a*a))
		yP2 = math.sqrt((u*u) + (a*a))
	tempbox = TBox(x1,y1,x2,y2)
	temppbox = TBox(x1,yP1,x2,yP2)
	tempsbox = TBox(x1,s1,x2,s2)
	Boxes.append(tempbox)
	sBoxes.append(tempsbox)
	pBoxes.append(temppbox)

Pull.GetXaxis().SetTitle("")
Pull.GetXaxis().SetLabelSize(0)
Pull.SetStats(0)
Pull.SetLineColor(1)
Pull.SetFillColor(0)
Pull.SetMarkerColor(1)
Pull.SetMarkerStyle(20)
#Pull.GetXaxis().SetNdivisions(0)
Pull.GetYaxis().SetNdivisions(4)
Pull.GetYaxis().SetTitle("#frac{Data - Bkg}{Stat. Unc.}")
Pull.GetYaxis().SetLabelSize(85/15*Pull.GetYaxis().GetLabelSize())
Pull.GetYaxis().SetTitleSize(5.2*Pull.GetYaxis().GetTitleSize())
Pull.GetYaxis().SetTitleOffset(0.205)
Pull.GetYaxis().SetRangeUser(-3.,3.)

for i in Boxes:
	i.SetFillColor(12)
	i.SetFillStyle(3244)
for i in pBoxes:
	i.SetFillColor(12)
	i.SetFillStyle(3144)
for i in sBoxes:
	i.SetFillColor(38)
	i.SetFillStyle(3002)

leg2 = TLegend(0.4,0.6,0.8,0.89)
leg2.SetLineColor(0)
leg2.SetFillColor(0)
leg2.SetTextSize(0.03)
leg2.AddEntry(V, "Data", "PLE")
leg2.AddEntry(N, "background prediction", "LF")
leg2.AddEntry(Boxes[0], "total uncertainty", "F")
leg2.AddEntry(sBoxes[0], "statistical uncertainty", "F")


T0 = TLine(800,0.,3000,0.)
T0.SetLineColor(kBlue)
T2 = TLine(800,2.,3000,2.)
T2.SetLineColor(kBlue)
T2.SetLineStyle(2)
Tm2 = TLine(800,-2.,3000,-2.)
Tm2.SetLineColor(kBlue)
Tm2.SetLineStyle(2)

T1 = TLine(800,1.,3000,1.)
T1.SetLineColor(kBlue)
T1.SetLineStyle(3)
Tm1 = TLine(800,-1.,3000,-1.)
Tm1.SetLineColor(kBlue)
Tm1.SetLineStyle(3)

C4 = TCanvas("C4", "", 800, 800)
#draw the lumi text on the canvas

plot = TPad("pad1", "The pad 80% of the height",0,0.15,1,1)
pull = TPad("pad2", "The pad 20% of the height",0,0,1.0,0.15)
plot.Draw()

gPad.SetTickx(1)
gPad.SetTicky(1) 
#pull.SetTickx(1)
pull.Draw()

plot.cd()
N.Draw("Hist")
for i in Boxes:
	i.Draw("same")
for i in sBoxes:
	i.Draw("same")
V.Draw("E0 same")
leg2.Draw()

l1 = TLatex()
l1.SetTextAlign(13)
l1.SetTextFont(42)
l1.SetNDC()
l1.SetTextSize(0.04)
l1.DrawLatex(0.54+0.03,0.25, "pp#rightarrow X #rightarrow HH #rightarrow b#bar{b}b#bar{b}")
l1.DrawLatex(0.54+0.03,0.3, "double-b tagging")

CMS_lumi.CMS_lumi(plot, iPeriod, iPos)
pull.cd()
gPad.SetTickx(1)
gPad.SetTicky(1)
Pull.Draw("")
for i in pBoxes:
	i.Draw("same")
T0.Draw("same")
T2.Draw("same")
Tm2.Draw("same")
T1.Draw("same")
Tm1.Draw("same")
Pull.Draw("same")
C4.Print("mx_unblind.pdf")

