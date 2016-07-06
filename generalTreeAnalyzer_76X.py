import os
import glob
import math
import ROOT
from ROOT import *

#ROOT.gROOT.Macro("rootlogon.C")

import FWCore.ParameterSet.Config as cms

import sys
from DataFormats.FWLite import Events, Handle

from array import *

ROOT.gSystem.Load('libCondFormatsBTagObjects') 

from optparse import OptionParser
parser = OptionParser()

parser.add_option("-f", "--pathIn", dest="inputFile",
                  help="inputFile path")

parser.add_option("-o", "--outName", dest="outName",
                  help="output file name")

parser.add_option("-i", "--min", dest="min", 
		  help="input index low end")

parser.add_option("-j", "--max", dest="max", 
		  help="input index high end")

parser.add_option("-l", "--file", dest="txt", 
		  help="input txt file")

parser.add_option("-t", "--trigger", dest="trigger", 
		  help="bool for trigger cut")

parser.add_option("-k", "--jets", dest="jets", 
		  help="bool for jet cuts")

parser.add_option("-d", "--deta", dest="deta", 
		  help="bool for delta eta cut")

parser.add_option("-m", "--isMC", dest="isMC", 
		  help="bool for is MC")

parser.add_option("-x", "--xsec", dest="xsec", 
		  help="cross section")

parser.add_option("-S", "--syst", dest="syst",
                  help="Systematic")


(options, args) = parser.parse_args()

inputfile = options.txt 

ff_n = 1000

num1 = int(options.min)
num2 = int(options.max)

d1 = options.outName 
d2 = '_'
print(options.outName)
outputfilename = d1 + d2 + options.syst + '.root'

print outputfilename

import copy
File_tr=ROOT.TFile.Open("trigger_objects.root", "R")
histo_efficiency=copy.copy(File_tr.Get("histo_efficiency"))
histo_efficiency_up=copy.copy(File_tr.Get("histo_efficiency_upper"))
histo_efficiency_down=copy.copy(File_tr.Get("histo_efficiency_lower"))
histo_efficiency_2up=copy.copy(File_tr.Get("histo_efficiency_upper_2sigma"))
histo_efficiency_2down=copy.copy(File_tr.Get("histo_efficiency_lower_2sigma"))
File_tr.Close()


def div_except(a, b):
    if b>0:
        return float(a)/b
    else:
        return 1


def btagging_efficiency_medium(pt):
   result = 0.898 + 0.0002254*pt -1.74e-6*pt*pt +2.71e-9*pt*pt*pt -1.39e-12*pt*pt*pt*pt
   return result	

def trigger_function(histo_efficiency,htJet30=700):
    result = histo_efficiency.GetBinContent(htJet30)
    return result


#defining functions
def ClosestJet(jets, fourvec): #returns the index of the jet (from a collection "jets") closest to the given four-vector
	DR = 9999.
	index = -1
	for j in range(0,len(jets)):
	    if jets[j].Pt() > 0 :
		dR = fourvec.DeltaR(jets[j])
		if dR < DR :
			DR = fourvec.DeltaR(jets[j])
			index = j
	return index

def MatchCollection(Col, jet): #matches a jet to a jet in a different collection
	j = -1
        dr = 0.4
	for i in range(len(Col)):
		C = ROOT.TLorentzVector()
                C.SetPtEtaPhiM( Col[i].Pt(), Col[i].Eta(), Col[i].Phi(), Col[i].M() )
		dr = abs(jet.DeltaR(C))
		if dr < 0.4 :
			#print "WOOHOO MATCH with index " + str(j) + " with dr " + str(dr)
			j = i
                        break
        if dr > 0.4:
	#	print "No Match :( for dr: " + str(dr)
#		print "index " + str(j)
		return -1
	return j

def MatchCollection2(Col, jet, index): #matches a jet to a jet in a different collection
	j = -1
        dr = 0.4
	for i in range(len(Col)):
            if i != index:
                C = ROOT.TLorentzVector()
                C.SetPtEtaPhiM( Col[i].Pt(), Col[i].Eta(), Col[i].Phi(), Col[i].M() )
		dr = abs(jet.DeltaR(C))
		if dr < 0.4 :
			#print "WOOHOO MATCH with index " + str(j) + " with dr " + str(dr)
			j = i
                        break
        if dr > 0.4:
	#	print "No Match :( for dr: " + str(dr)
#		print "index " + str(j)
		return -1
	return j

def MatchCollection3(Col, jet, index1, index2): #matches a jet to a jet in a different collection
	j = -1
        dr = 0.4
	for i in range(len(Col)):
            if i != index1 and i != index2:
                C = ROOT.TLorentzVector()
                C.SetPtEtaPhiM( Col[i].Pt(), Col[i].Eta(), Col[i].Phi(), Col[i].M() )
		dr = abs(jet.DeltaR(C))
		if dr < 0.4 :
			#print "WOOHOO MATCH with index " + str(j) + " with dr " + str(dr)
			j = i
                        break
        if dr > 0.4:
	#	print "No Match :( for dr: " + str(dr)
#		print "index " + str(j)
		return -1
	return j

def MatchCollection4(Col, jet, index1, index2, index3): #matches a jet to a jet in a different collection
	j = -1
        dr = 0.4
	for i in range(len(Col)):
            if i != index1 and i != index2 and i != index3:
                C = ROOT.TLorentzVector()
                C.SetPtEtaPhiM( Col[i].Pt(), Col[i].Eta(), Col[i].Phi(), Col[i].M() )
		dr = abs(jet.DeltaR(C))
		if dr < 0.4 :
			#print "WOOHOO MATCH with index " + str(j) + " with dr " + str(dr)
			j = i
                        break
        if dr > 0.4:
	#	print "No Match :( for dr: " + str(dr)
#		print "index " + str(j)
		return -1
	return j

def open_files(file_name) : #opens files to run on

    g = open(file_name)
    list_file = []
    final_list = []
    for i in range(ff_n):  # this is the length of the file
        list_file.append(g.readline().split())
    s = options.inputFile

    for i in range(len(list_file)):
        for j in range(len(list_file[i])) :
            final_list.append(s + list_file[i][j])
  #  print final_list
    return final_list


def deltaR( particle, jet ) : #gives deltaR between two particles
    DeltaPhiHere = math.fabs( particle.phi() - jet.phi() )
    if DeltaPhiHere > math.pi :
        DeltaPhiHere = 2*math.pi - DeltaPhiHere

    deltaRHere = math.sqrt( (particle.eta() - jet.eta() )**2 + ( DeltaPhiHere )**2  )
    return deltaRHere

print sys.argv[1]

f =  ROOT.TFile(outputfilename, 'recreate')

f.cd()


myTree =  ROOT.TTree('myTree', 'myTree')

#first tree (boosted case)
#creating the tree objects we need
jet1pt = array('f', [-100.0])
jet2pt = array('f', [-100.0])
jet1ID = array('f', [-100.0])
jet2ID = array('f', [-100.0])
jet1eta = array('f', [-100.0])
jet2eta = array('f', [-100.0])
etadiff = array('f', [-100.0])
dijetmass = array('f', [-100.0])
dijetmass_corr = array('f', [-100.0])
dijetmass_corr_punc = array('f', [-100.0])
jet1tau21 = array('f', [-100.0])
jet2tau21 = array('f', [-100.0])
jet1pmass = array('f', [-100.0])
jet2pmass = array('f', [-100.0])
jet1pmassunc = array('f', [-100.0])
jet2pmassunc = array('f', [-100.0])
#pjet1pmass = array('f', [-100.0])
#pjet2pmass = array('f', [-100.0])
jet1bbtag = array('f', [-100.0])
jet2bbtag = array('f', [-100.0])
jet1s1csv = array('f', [-100.0])
jet2s1csv = array('f', [-100.0])
jet1s2csv = array('f', [-100.0])
jet2s2csv = array('f', [-100.0])
jetSJfla = array('f', [-100.0]*4)
jetSJpt =  array('f', [-100.0]*4)
jetSJcsv = array('f', [-100.0]*4)
jetSJeta = array('f', [-100.0]*4)
triggerpassbb = array('f', [-100.0])
triggerpasssj = array('f', [-100.0])
nHiggsTags = array('f', [-100.0])
nTrueInt = array('f', [-100])
#PUWeight  = array('f', [-100.0])
vtype = array('f', [-100.0])
isData = array('f', [100.0])
jet1nbHadron = array('f', [-100.0])
jet2nbHadron = array('f', [-100.0])
jet1flavor = array('f', [-100.0])
jet2flavor = array('f', [-100.0])
jet1ncHadron = array('f', [-100.0])
jet2ncHadron = array('f', [-100.0])
gen1Pt = array('f', [-100.0])
gen1Phi = array('f', [-100.0])
gen1Eta = array('f', [-100.0])
gen1Mass = array('f', [-100.0])
gen1ID = array('f', [-100.0])
gen2Pt = array('f', [-100.0])
gen2Phi = array('f', [-100.0])
gen2Eta = array('f', [-100.0])
gen2Mass = array('f', [-100.0])
gen2ID = array('f', [-100.0])
jet1l1l2l3 = array('f', [-100.0])
jet1l2l3 = array('f', [-100.0])
jet2l1l2l3 = array('f', [-100.0])
jet2l2l3 = array('f', [-100.0])
jet1JER = array('f', [-100.0])
jet2JER = array('f', [-100.0])
puWeights = array('f', [-100.0])
puWeightsUp = array('f', [-100.0])
puWeightsDown = array('f', [-100.0])
json = array('f', [-100.0])
SF = array('f', [-100.0])
SFup = array('f', [-100.0])
SFdown = array('f', [-100.0])
SF4sj = array('f', [-100.0])
SF4sjUp = array('f', [-100.0])
SF4sjDown = array('f', [-100.0])
SF3sj = array('f', [-100.0])
SF3sjUp = array('f', [-100.0])
SF3sjDown = array('f', [-100.0])
trigWeight = array('f', [-100.0])
trigWeightUp = array('f', [-100.0])
trigWeightDown = array('f', [-100.0])
trigWeight2Up = array('f', [-100.0])
trigWeight2Down = array('f', [-100.0])
norm = array('f', [-100.0])
evt = array('f', [-100.0])
ht = array('f', [-100.0])
xsec = array('f', [-100.0])
sjSF = array('f', [-100.0])
sjSFup = array('f', [-100.0])
sjSFdown = array('f', [-100.0])

#creating the tree branches we need
myTree.Branch('jet1pt', jet1pt, 'jet1pt/F')
myTree.Branch('jet2pt', jet2pt, 'jet2pt/F')
myTree.Branch('jet1eta', jet1eta, 'jet1eta/F')
myTree.Branch('jet2eta', jet2eta, 'jet2eta/F')
myTree.Branch('etadiff', etadiff, 'etadiff/F')
myTree.Branch('dijetmass', dijetmass, 'dijetmass/F')
myTree.Branch('dijetmass_corr', dijetmass_corr, 'dijetmass_corr/F')
myTree.Branch('dijetmass_corr_punc', dijetmass_corr_punc, 'dijetmass_corr_punc/F')
myTree.Branch('jet1tau21', jet1tau21, 'jet1tau21/F')
myTree.Branch('jet2tau21', jet2tau21, 'jet2tau21/F')
myTree.Branch('jet1pmass', jet1pmass, 'jet1pmass/F')
myTree.Branch('jet2pmass', jet2pmass, 'jet2pmass/F')
myTree.Branch('jet1pmassunc', jet1pmassunc, 'jet1pmassunc/F')
myTree.Branch('jet2pmassunc', jet2pmassunc, 'jet2pmassunc/F')
myTree.Branch('jet1bbtag', jet1bbtag, 'jet1bbtag/F')
myTree.Branch('jet2bbtag', jet2bbtag, 'jet2bbtag/F')
myTree.Branch('jet1s1csv', jet1s1csv, 'jet1s1csv/F')
myTree.Branch('jet2s1csv', jet2s1csv, 'jet2s1csv/F')
myTree.Branch('jet1s2csv', jet1s2csv, 'jet1s2csv/F')
myTree.Branch('jet2s2csv', jet2s2csv, 'jet2s2csv/F')
myTree.Branch('nHiggsTags', nHiggsTags, 'nHiggsTags/F')
myTree.Branch('triggerpassbb', triggerpassbb, 'triggerpassbb/F')
myTree.Branch('triggerpasssj', triggerpasssj, 'triggerpasssj/F')
myTree.Branch('nTrueInt',nTrueInt,'nTrueInt/F')
myTree.Branch('puWeights',puWeights,'puWeights/F')
myTree.Branch('puWeightsUp',puWeightsUp,'puWeightsUp/F')
myTree.Branch('puWeightsDown',puWeightsDown,'puWeightsDown/F')
myTree.Branch('jet1ID', jet1ID, 'jet1ID/F')
myTree.Branch('jet2ID', jet2ID, 'jet2ID/F')
myTree.Branch('vtype', vtype, 'vtype/F') 
myTree.Branch('isData', isData, 'isData/F') 
myTree.Branch('jet1nbHadron', jet1nbHadron, 'jet1nbHadron/F')
myTree.Branch('jet2nbHadron', jet2nbHadron, 'jet2nbHadron/F')
myTree.Branch('jet1flavor', jet1flavor, 'jet1flavor/F') 
myTree.Branch('jet2flavor', jet2flavor, 'jet2flavor/F') 
myTree.Branch('jet1ncHadron', jet1ncHadron, 'jet1ncHadron/F')
myTree.Branch('jet2ncHadron', jet2ncHadron, 'jet2ncHadron/F')
myTree.Branch('gen1Pt', gen1Pt, 'gen1Pt/F')
myTree.Branch('gen1Phi', gen1Phi, 'gen1Phi/F')
myTree.Branch('gen1Eta', gen1Eta, 'gen1Eta/F')
myTree.Branch('gen1Mass', gen1Mass, 'gen1Mass/F')
myTree.Branch('gen1ID', gen1ID, 'gen1ID/F')
myTree.Branch('gen2Pt', gen2Pt, 'gen2Pt/F')
myTree.Branch('gen2Phi', gen2Phi, 'gen2Phi/F')
myTree.Branch('gen2Eta', gen2Eta, 'gen2Eta/F')
myTree.Branch('gen2Mass', gen2Mass, 'gen2Mass/F')
myTree.Branch('gen2ID', gen2ID, 'gen2ID/F')
myTree.Branch('jet1l1l2l3', jet1l1l2l3, 'jet1l1l2l3/F') 
myTree.Branch('jet1l2l3', jet1l2l3, 'jet1l2l3/F')
myTree.Branch('jet2l1l2l3', jet2l1l2l3, 'jet2l1l2l3/F') 
myTree.Branch('jet2l2l3', jet2l2l3, 'jet2l2l3/F')
myTree.Branch('jet1JER', jet1JER, 'jet1JER/F') 
myTree.Branch('jet2JER', jet2JER, 'jet2JER/F') 
myTree.Branch('json', json, 'json/F')
myTree.Branch('SF', SF, 'SF/F')
myTree.Branch('SFup', SFup, 'SFup/F')
myTree.Branch('SFdown', SFdown, 'SFdown/F')
myTree.Branch('SF4sj', SF4sj, 'SF4sj/F')
myTree.Branch('SF4sjUp', SF4sjUp, 'SF4sjUp/F')
myTree.Branch('SF4sjDown', SF4sjDown, 'SF4sjDown/F')
myTree.Branch('SF3sj', SF3sj, 'SF3sj/F')
myTree.Branch('SF3sjUp', SF3sjUp, 'SF3sjUp/F')
myTree.Branch('SF3sjDown', SF3sjDown, 'SF3sjDown/F')
myTree.Branch('trigWeight', trigWeight, 'trigWeight/F')
myTree.Branch('trigWeightUp', trigWeightUp, 'trigWeightUp/F')
myTree.Branch('trigWeightDown', trigWeightDown, 'trigWeightDown/F')
myTree.Branch('trigWeight2Up', trigWeight2Up, 'trigWeight2Up/F')
myTree.Branch('trigWeight2Down', trigWeight2Down, 'trigWeight2Down/F')
myTree.Branch('norm',norm,'norm/F')
myTree.Branch('evt',evt,'evt/F')
myTree.Branch('ht', ht, 'ht/F')
myTree.Branch('xsec', xsec, 'xsec/F')
myTree.Branch('sjSF', sjSF, 'sjSF/F')
myTree.Branch('sjSFup', sjSFup, 'sjSFup/F')
myTree.Branch('sjSFdown', sjSFdown, 'sjSFdown/F')
myTree.Branch('jetSJfla',jetSJfla,'jetSJfla[4]/F') 
myTree.Branch('jetSJpt', jetSJpt,'jetSJpt[4]/F')
myTree.Branch('jetSJcsv',jetSJcsv,'jetSJcsv[4]/F')
myTree.Branch('jetSJeta',jetSJeta,'jetSJeta[4]/F')

files_list	= open_files( inputfile )
#nevent = treeMine.GetEntries();

#list of histograms that may be useful
bbj = ROOT.TH1F("bbj", "Before any cuts", 3, -0.5, 1.5)
bb0 = ROOT.TH1F("bb0", "After Json", 3, -0.5, 1.5)
bb1 = ROOT.TH1F("bb1", "After Trigger", 3, -0.5, 1.5)
bb2 = ROOT.TH1F("bb2", "After jet cuts", 3, -0.5, 1.5)
bb3 = ROOT.TH1F("bb3", "After Delta Eta cuts", 3, -0.5, 1.5)

calib = ROOT.BTagCalibration("csvv2","/uscms_data/d3/cvernier/DiH_13TeV/optimization/Alphabet-76x/CSVv2_subjet.csv")
readerHF = ROOT.BTagCalibrationReader(calib,0, "lt","central")  # 0 is for loose op
readerHFup = ROOT.BTagCalibrationReader(calib, 0,"lt", "up")  # 0 is for loose op
readerHFdown = ROOT.BTagCalibrationReader(calib, 0,"lt", "down")  # 0 is for loose op
readerLF = ROOT.BTagCalibrationReader(calib,0, "incl","central")  # 0 is for loose op
readerLFup = ROOT.BTagCalibrationReader(calib, 0,"incl", "up")  # 0 is for loose op
readerLFdown = ROOT.BTagCalibrationReader(calib, 0,"incl", "down")  # 0 is for loose op




gSystem.Load("DrawFunctions_h.so")

count = 0
#loop over files
for i in range(num1, num2):
    files = files_list[i]
    print files
    f1 = ROOT.TFile(files, "READ")
    treeMine  = f1.Get('tree')
    nevent = treeMine.GetEntries();

    #getting the norm
    if options.isMC:
        histo_weight=f1.Get("CountWeighted")
        norm[0]=histo_weight.GetBinContent(1)
    else:
        norm[0] = 1
    #loop over events in file
    print "Start looping"
    for j in range(0,nevent):
        treeMine.GetEntry(j)
	count = count + 1
 	if count % 1000 == 0 :
	    print "processing events", count
	
	#variables we need from the heppy ntuple
	fJetPt  = treeMine.Jet_pt
        fJetEta  = treeMine.Jet_eta
        fNJets = treeMine.nJet
        genPt = treeMine.GenJet_pt
        genEta = treeMine.GenJet_eta
        genPhi = treeMine.GenJet_phi
        genMass = treeMine.GenJet_mass
	genBH = treeMine.GenJet_numBHadrons	
	genCH = treeMine.GenJet_numCHadrons
	fjUngroomedN = treeMine.nFatjetAK08ungroomed
        fjUngroomedPt = treeMine.FatjetAK08ungroomed_pt
	fjUngroomedEta = treeMine.FatjetAK08ungroomed_eta
	fjUngroomedPhi = treeMine.FatjetAK08ungroomed_phi
	fjUngroomedMass = treeMine.FatjetAK08ungroomed_mass
	fjUngroomedSDMass = treeMine.FatjetAK08ungroomed_msoftdrop
	fjUngroomedTau1 = treeMine.FatjetAK08ungroomed_tau1
	fjUngroomedTau2 = treeMine.FatjetAK08ungroomed_tau2
	fjUngroomedBbTag = treeMine.FatjetAK08ungroomed_bbtag
	fjUngroomedJetID = treeMine.FatjetAK08ungroomed_id_Tight
	fjUngroomedPrunedMass = treeMine.FatjetAK08ungroomed_mprunedcorr
        fjUngroomedPrunedMass_Unc = treeMine.FatjetAK08ungroomed_mpruned
        if options.isMC:
            fjUngroomedFlavour = treeMine.FatjetAK08ungroomed_Flavour
            fjUngroomedBHadron = treeMine.FatjetAK08ungroomed_BhadronFlavour
            fjUngroomedCHadron = treeMine.FatjetAK08ungroomed_ChadronFlavour
            fjUngroomedJER = treeMine.FatjetAK08ungroomed_GenPt
        fjL2L3 = treeMine.FatjetAK08ungroomed_JEC_L2L3
        fjL1L2L3 = treeMine.FatjetAK08ungroomed_JEC_L1L2L3
        if options.isMC:
            puweight = treeMine.puWeight 
            puweightUp = treeMine.puWeightUp
            puweightDown = treeMine.puWeightDown
	sjPrunedPt = treeMine.SubjetAK08softdrop_pt
	sjPrunedEta = treeMine.SubjetAK08softdrop_eta
	sjPrunedPhi = treeMine.SubjetAK08softdrop_phi
	sjPrunedMass = treeMine.SubjetAK08softdrop_mass
	sjPrunedBtag = treeMine.SubjetAK08softdrop_btag
        if options.isMC:
            hPt = treeMine.GenHiggsBoson_pt
            hEta = treeMine.GenHiggsBoson_eta
            hPhi = treeMine.GenHiggsBoson_phi
            hMass = treeMine.GenHiggsBoson_mass
	hltHT800 = treeMine.HLT_BIT_HLT_PFHT800_v
        hltHT650_MJJ900 = treeMine.HLT_BIT_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v
        hltHT650_MJJ950 = treeMine.HLT_BIT_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v
        hltAK8PFJet360 = treeMine.HLT_BIT_HLT_AK8PFJet360_TrimMass30_v
        hltAK8PFHT700 = treeMine.HLT_BIT_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v
        Data = treeMine.isData
        vType = treeMine.Vtype
        EVT = treeMine.evt 
        if options.isMC:
            nTInt = treeMine.nTrueInt
        genTopPts = treeMine.GenTop_pt
        JSON = treeMine.json
        
      
	#saving whether an event passes desired trigger (bb = HT800 pass, sj = pass any of the five saved triggers
        matched = 0    
	if hltHT800 > 0:
            matched += 1
        triggerpassbb[0] = matched
        
        matchedsj = 0
        if hltHT800 > 0:
            matchedsj += 1
        if hltHT650_MJJ900 > 0:
            matchedsj += 1
        if hltHT650_MJJ950 > 0:
            matchedsj += 1
        if hltAK8PFJet360 > 0:
            matchedsj += 1
        if hltAK8PFHT700 > 0:
            matchedsj += 1
        triggerpasssj[0] = matchedsj

        #trigger weights
	hT =0
        for i in range(0,fNJets):
                if abs(fJetEta[i])<3 and fJetPt[i] >40 :
                        hT=hT+fJetPt[i]

        ht[0] = hT 
        trigWeight[0] = trigger_function(histo_efficiency, int(round(hT)))
        trigWeightUp[0] = trigger_function(histo_efficiency_up, int(round(hT)))    
        trigWeightDown[0] = trigger_function(histo_efficiency_down, int(round(hT)))    
        trigWeight2Up[0] = trigger_function(histo_efficiency_2up, int(round(hT)))   
        trigWeight2Down[0] = trigger_function(histo_efficiency_2down, int(round(hT)))                     
        #json for data
        bbj.Fill(triggerpassbb[0])
        if Data and treeMine.json_silver < 1:
            continue		
        bb0.Fill(triggerpassbb[0])

        #requiring event pass trigger
#        if options.trigger and triggerpass[0] < 1:
#            continue

        bb1.Fill(triggerpassbb[0])

	#filling an array with jet 4-vectors for jets pt > 30 and |eta| < 2.5, an array of tau21s, and an array of bbtag values, pmass, id, nbhadrons, nchadrons, flavor, l1l2l3 corr, l2l3 corr, JER
        jets = []
	jet_tau = []
	jet_bbtag = []
        jet_pmass = []
        jet_pmassunc = []
        jet_id = []
        jet_nb = []
        jet_nc = []
        jet_flav = []
        jet_123 = []
        jet_23 = []
        jet_JER = []
        for j in range(len(fjUngroomedPt)):
            jettemp = ROOT.TLorentzVector()
            jettemp.SetPtEtaPhiM(fjUngroomedPt[j], fjUngroomedEta[j], fjUngroomedPhi[j], fjUngroomedMass[j])
	    if (options.syst=="FJEC_Up"):
                            correction_factor=1+(tree.FatjetAK08ungroomed_JEC_UP[j]-tree.FatjetAK08ungroomed_JEC_L1L2L3[j])
                            jettemp*=correction_factor
	    if (options.syst=="FJEC_Down"):
                            correction_factor=1-(tree.FatjetAK08ungroomed_JEC_UP[j]-tree.FatjetAK08ungroomed_JEC_L1L2L3[j])
                            jettemp*=correction_factor
	    if (options.syst=="FJER_Up"):
                            correction_factor=div_except(tree.FatjetAK08ungroomed_JER_UP_PT[j],tree.FatjetAK08ungroomed_pt[j])
                            jettemp*=correction_factor
	    if (options.syst=="FJER_Down"):
                            pJERDown=2*tree.FatjetAK08ungroomed_pt[j]-tree.FatjetAK08ungroomed_JER_UP_PT[j]
			    correction_factor=div_except((pJERDown),tree.FatjetAK08ungroomed_pt[j])
			    jettemp*=correction_factor
	
	
	    if jettemp.Pt() > 300. and abs(jettemp.Eta()) < 2.4: 	
                    jets.append(jettemp)
		    if fjUngroomedTau1[j] > 0:
			    jet_tau.append(fjUngroomedTau2[j]/fjUngroomedTau1[j])
		    else:
			    jet_tau.append(100)
		    mpruned_syst=fjUngroomedPrunedMass[j]
		    if (options.syst=="MJEC_Down"):
                            sigma=tree.FatjetAK08ungroomed_JEC_L2L3_UP[j]-tree.FatjetAK08ungroomed_JEC_L2L3[j]
                            mpruned_syst=tree.FatjetAK08ungroomed_mpruned[j]*(tree.FatjetAK08ungroomed_JEC_L2L3[j]-sigma)
                    if (options.syst=="MJEC_Up"): 
			mpruned_syst=tree.FatjetAK08ungroomed_mpruned[j]*tree.FatjetAK08ungroomed_JEC_L2L3_UP[j]

		    jet_bbtag.append(fjUngroomedBbTag[j])	
                    jet_pmass.append(mpruned_syst)
                    jet_pmassunc.append(fjUngroomedPrunedMass_Unc[j])
                    jet_id.append(fjUngroomedJetID[j])
                    if options.isMC:
                        jet_nb.append(fjUngroomedBHadron[j])
                        jet_nc.append(fjUngroomedCHadron[j])
                        jet_flav.append(fjUngroomedFlavour[j])
                        jet_JER.append(fjUngroomedJER[j])
                    jet_123.append(fjL1L2L3[j])
                    jet_23.append(fjL2L3[j])
                        
                        
                            
                        

	if options.jets and len(jets) < 2: # two jets with pt > 30 and |eta| < 2.5
		continue

        bb2.Fill(triggerpassbb[0])

	#dEta selection : selecting the two jets which minimizes the dEta requirement. (to find a better one?)
	idxH1 = -1
	idxH2=-1
	if abs(jets[0].Eta() - jets[1].Eta()) < 1.3:
		minDEta = abs(jets[0].Eta() - jets[1].Eta())
		idxH1 = 0
		idxH2 = 1
			
	if options.deta and (idxH1 < 0 or idxH2 <0) : continue
  
        bb3.Fill(triggerpassbb[0])

	#higgs tagging - matching higgs gen jet to the 1 and 2 pt jet
	if options.isMC:
            hjets = []
            for j in range(len(hPt)):
		jettemp = ROOT.TLorentzVector()
		jettemp.SetPtEtaPhiM(hPt[j], hEta[j], hPhi[j], hMass[j])
		hjets.append(jettemp)

            h1 = MatchCollection(hjets, jets[idxH1])
            h2 = MatchCollection2(hjets, jets[idxH2],h1)
            nHiggsTags[0] = 0
            if h1 > -1:
		nHiggsTags[0] += 1
            if h2 > -1:
		nHiggsTags[0] += 1

        #filling jet variables
        jet1pmass[0] = jet_pmass[idxH1]
        jet2pmass[0] = jet_pmass[idxH2]
        jet1pmassunc[0] = jet_pmassunc[idxH1]
        jet2pmassunc[0] = jet_pmassunc[idxH2]
        jet1ID[0] = jet_id[idxH1]
	jet2ID[0] = jet_id[idxH2]
	jet1tau21[0] = jet_tau[idxH1]
        jet2tau21[0] = jet_tau[idxH2]
        if options.isMC:
            jet1nbHadron[0] = jet_nb[idxH1]
            jet2nbHadron[0] = jet_nb[idxH2]
            jet1ncHadron[0] = jet_nc[idxH1]
            jet2ncHadron[0] = jet_nc[idxH2]
            jet1flavor[0] = jet_flav[idxH1]
            jet2flavor[0] = jet_flav[idxH2]
            jet1JER[0] = jet_JER[idxH1]
            jet2JER[0] = jet_JER[idxH2]
        jet1l1l2l3[0] = jet_123[idxH1]
        jet2l1l2l3[0] = jet_123[idxH2]
        jet1l2l3[0] = jet_23[idxH1]
        jet2l2l3[0] = jet_23[idxH2]

        #finding gen jets to match higgs jets
        if options.isMC:
            ujets = []
	    ujetsCH = []
	    ujetsBH = []		
            for j in range(len(genPt)):
                jettemp = ROOT.TLorentzVector()
                jettemp.SetPtEtaPhiM(genPt[j], genEta[j], genPhi[j], genMass[j])
                ujets.append(jettemp)
		ujetsCH.append(genCH[j])
		ujetsBH.append(genBH[j])

            j1 = MatchCollection(ujets, jets[idxH1])
            j2 = MatchCollection2(ujets, jets[idxH2],j1)
           
            #filling gen jet info
            gen1Pt[0] = ujets[j1].Pt()
            gen1Phi[0] = ujets[j1].Phi()
            gen1Eta[0] = ujets[j1].Eta()
            gen1Mass[0] = ujets[j1].M()
            gen1ID[0] = j1
            gen2Pt[0] = ujets[j2].Pt()
            gen2Phi[0] = ujets[j2].Phi()
            gen2Eta[0] = ujets[j2].Eta()
            gen2Mass[0] = ujets[j2].M()
            gen2ID[0] = j2

	#filling min subjet csv
	subjets = []
	jet1sj = []
	jet1sjcsv = []
	jet2sj = []
	jet2sjcsv = []
	samesj = 0
        for j in range(len(sjPrunedPt)):
            jettemp = ROOT.TLorentzVector()
            jettemp.SetPtEtaPhiM(sjPrunedPt[j], sjPrunedEta[j], sjPrunedPhi[j], sjPrunedMass[j])
            subjets.append(jettemp)

	for j in range(len(subjets)):
            dR1 = subjets[j].DeltaR(jets[idxH1])
	    dR2 = subjets[j].DeltaR(jets[idxH2])
	    if dR1 < 0.4 and dR2 < 0.4:
		    samesj += 1
	    elif dR1 < 0.4:
		    jet1sj.append(subjets[j])
		    jet1sjcsv.append(sjPrunedBtag[j])
	    elif dR2 < 0.4:
		    jet2sj.append(subjets[j])
		    jet2sjcsv.append(sjPrunedBtag[j])
	n1sj = len(jet1sj)
	n2sj = len(jet2sj)

	#Finding the subjet csvs
	jet1s1csv[0] = -1.
	jet2s1csv[0] = -1.
        jet1s2csv[0] = -1.
	jet2s2csv[0] = -1.
	
	for i in range(0,4):
	  jetSJfla[i] =-1
	  jetSJpt[i]  =-1
	  jetSJcsv[i] =-1
	  jetSJeta[i] =-1

        
        if len(jet1sjcsv) > 1:
            jet1s1csv[0] = jet1sjcsv[0]
            jet1s2csv[0] = jet1sjcsv[1]
        elif len(jet1sjcsv) == 1:
            jet1s1csv[0] = jet1sjcsv[0]

        if len(jet2sjcsv) > 1:
            jet2s1csv[0] = jet2sjcsv[0]
            jet2s2csv[0] = jet2sjcsv[1]
        elif len(jet2sjcsv) == 1:
            jet2s1csv[0] = jet2sjcsv[0]
        sfsj3 =-1
        sfsj4 =-1
	sfsj1 =-1
        sfsj2 =-1	
	sfsj3up =-1
        sfsj4up =-1
        sfsj1up =-1
        sfsj2up =-1
	sfsj3down =-1
        sfsj4down =-1
        sfsj1down =-1
        sfsj2down =-1
        #finding gen jets for subjets
        if options.isMC:
            if len(jet1sjcsv) > 1:
                sj1gen = MatchCollection(ujets, jet1sj[0])
                sj2gen = MatchCollection2(ujets, jet1sj[1],sj1gen)
		isL = False
		isL2 = False
                if sj1gen>0 and ujetsBH[sj1gen]>0 :
                    sj1flav = BTagEntry.FLAV_B
		    jetSJfla[0] = 1
		    jetSJpt[0] = jet1sj[0].Pt()
        	    jetSJcsv[0] = jet1sjcsv[0]
                    jetSJeta[0] = jet1sj[0].Eta()
                elif ujetsCH[sj1gen]>0 and sj1gen>0:
                    sj1flav = BTagEntry.FLAV_C
                else:
		    isL= True 
		    sj1flav = BTagEntry.FLAV_UDSG	
                if  ujetsBH[sj2gen]>0 and sj2gen>0:
                    sj2flav = BTagEntry.FLAV_B
		    jetSJpt[1] = jet1sj[1].Pt()
           	    jetSJeta[1] = jet1sj[1].Eta()
                    jetSJcsv[1] = jet1sjcsv[1]
		    jetSJfla[1] = 1
                elif ujetsCH[sj2gen]>0 and sj2gen>0:
                    sj2flav = BTagEntry.FLAV_C
                else:
		    sj2flav = BTagEntry.FLAV_UDSG
		    isL2 = True
	   #compute SF
	        if not isL:	
	         if(jet1sj[0].Pt()<670. and abs(jet1sj[0].Eta())<2.4 ) :

                  sfsj1 = readerHF.eval(sj1flav, jet1sj[0].Eta(), jet1sj[0].Pt())  # jet flavor, eta, pt
		  sfsj1up = readerHFup.eval(sj1flav, jet1sj[0].Eta(), jet1sj[0].Pt()) 
		  sfsj1down = readerHFdown.eval(sj1flav, jet1sj[0].Eta(), jet1sj[0].Pt())
		 elif abs(jet1sj[0].Eta())>2.4:
		  sfsj1 = readerHF.eval(sj1flav, 2.399, jet1sj[0].Pt()) 
		  sfsj1up = readerHFup.eval(sj1flav, 2.399, jet1sj[0].Pt())	
		  sfsj1down = readerHFdown.eval(sj1flav, 2.399, jet1sj[0].Pt())
		 elif jet1sj[0].Pt()>670.:
		  sfsj1 = readerHF.eval(sj1flav, jet1sj[0].Eta(), 669.)   
		  sfsj1up = readerHFup.eval(sj1flav, jet1sj[0].Eta(), 669.)
		  sfsj1down = readerHFdown.eval(sj1flav, jet1sj[0].Eta(), 669.)
		else:
		 if(jet1sj[0].Pt()<670. and abs(jet1sj[0].Eta())<2.4 ) :
                  sfsj1 = readerLF.eval(sj1flav, jet1sj[0].Eta(), jet1sj[0].Pt())  # jet flavor, eta, pt
                  sfsj1up = readerLFup.eval(sj1flav, jet1sj[0].Eta(), jet1sj[0].Pt())
                  sfsj1down = readerLFdown.eval(sj1flav, jet1sj[0].Eta(), jet1sj[0].Pt())
                 elif abs(jet1sj[0].Eta())>2.4:
                  sfsj1 = readerLF.eval(sj1flav, 2.399, jet1sj[0].Pt())
                  sfsj1up = readerLFup.eval(sj1flav, 2.399, jet1sj[0].Pt())
                  sfsj1down = readerLFdown.eval(sj1flav, 2.399, jet1sj[0].Pt())
                 elif jet1sj[0].Pt()>670.:
                  sfsj1 = readerLF.eval(sj1flav, jet1sj[0].Eta(), 669.)
                  sfsj1up = readerLFup.eval(sj1flav, jet1sj[0].Eta(), 669.)
                  sfsj1down = readerLFdown.eval(sj1flav, jet1sj[0].Eta(), 669.)
		if not isL2:
		 if(jet1sj[1].Pt()<670. and abs(jet1sj[1].Eta())<2.4 ) :
                  sfsj2 = readerHF.eval(sj2flav, jet1sj[1].Eta(), jet1sj[1].Pt())  # jet flavor, eta, pt
                  sfsj2up = readerHFup.eval(sj2flav, jet1sj[1].Eta(), jet1sj[1].Pt())
                  sfsj2down = readerHFdown.eval(sj2flav, jet1sj[1].Eta(), jet1sj[1].Pt())
                 elif abs(jet1sj[1].Eta())>2.4:
                  sfsj2 = readerHF.eval(sj2flav, 2.399, jet1sj[1].Pt())
                  sfsj2up = readerHFup.eval(sj2flav, 2.399, jet1sj[1].Pt())
                  sfsj2down = readerHFdown.eval(sj2flav, 2.399, jet1sj[1].Pt())
                 elif jet1sj[1].Pt()>670.:
                  sfsj2 = readerHF.eval(sj2flav, jet1sj[1].Eta(), 669.)
                  sfsj2up = readerHFup.eval(sj2flav, jet1sj[1].Eta(), 669.)
                  sfsj2down = readerHFdown.eval(sj2flav, jet1sj[1].Eta(), 669.)
	        else:	
		 if(jet1sj[1].Pt()<670. and abs(jet1sj[1].Eta())<2.4 ) :
                  sfsj2 = readerLF.eval(sj2flav, jet1sj[1].Eta(), jet1sj[1].Pt())  # jet flavor, eta, pt
                  sfsj2up = readerLFup.eval(sj2flav, jet1sj[1].Eta(), jet1sj[1].Pt())
                  sfsj2down = readerLFdown.eval(sj2flav, jet1sj[1].Eta(), jet1sj[1].Pt())
                 elif abs(jet1sj[1].Eta())>2.4:
                  sfsj2 = readerLF.eval(sj2flav, 2.399, jet1sj[1].Pt())
                  sfsj2up = readerLFup.eval(sj2flav, 2.399, jet1sj[1].Pt())
                  sfsj2down = readerLFdown.eval(sj2flav, 2.399, jet1sj[1].Pt())
                 elif jet1sj[1].Pt()>670.:
                  sfsj2 = readerLF.eval(sj2flav, jet1sj[1].Eta(), 669.)
                  sfsj2up = readerLFup.eval(sj2flav, jet1sj[1].Eta(), 669.)
                  sfsj2down = readerLFdown.eval(sj2flav, jet1sj[1].Eta(), 669.)
	    
	    if len(jet2sjcsv) > 1:
                sj3gen = MatchCollection(ujets, jet2sj[0])
                sj4gen = MatchCollection2(ujets, jet2sj[1],sj3gen)
		isL = False
                isL2 = False
                if  ujetsBH[sj3gen]>0 and sj3gen>0.:
                    sj3flav = BTagEntry.FLAV_B
		    jetSJpt[2] = jet2sj[0].Pt()
                    jetSJeta[2] = jet2sj[0].Eta()
                    jetSJcsv[2] = jet2sjcsv[0]
		    jetSJfla[2] = 1
                elif ujetsCH[sj3gen]>0 and sj3gen>0.:
                    sj3flav = BTagEntry.FLAV_C
                else:
                    sj3flav = BTagEntry.FLAV_UDSG
	 	    isL = True
                if  ujetsBH[sj4gen]>0 and sj4gen>0.:
                    sj4flav = BTagEntry.FLAV_B
		    jetSJpt[3] = jet2sj[1].Pt()
                    jetSJeta[3] = jet2sj[1].Eta()
                    jetSJcsv[3] = jet2sjcsv[1]
                    jetSJfla[3] = 1
                elif ujetsCH[sj4gen]>0 and sj4gen>0.:
                    sj4flav = BTagEntry.FLAV_C
                else:
                    sj4flav = BTagEntry.FLAV_UDSG
		    isL2 = True
           #compute SF
		if not isL:
                 if(jet2sj[0].Pt()<670. and abs(jet2sj[0].Eta())<2.4 ) :
                  sfsj3 = readerHF.eval(sj3flav, jet2sj[0].Eta(), jet2sj[0].Pt())  # jet flavor, eta, pt
                  sfsj3up = readerHFup.eval(sj3flav, jet2sj[0].Eta(), jet2sj[0].Pt())
                  sfsj3down = readerHFdown.eval(sj3flav, jet2sj[0].Eta(), jet2sj[0].Pt())
                 elif abs(jet2sj[0].Eta())>2.4:
                  sfsj3 = readerHF.eval(sj3flav, 2.399, jet2sj[0].Pt())
                  sfsj3up = readerHFup.eval(sj3flav, 2.399, jet2sj[0].Pt())
                  sfsj3down = readerHFdown.eval(sj3flav, 2.399, jet2sj[0].Pt())
                 elif jet2sj[0].Pt()>670.:
                  sfsj3 = readerHF.eval(sj3flav, jet2sj[0].Eta(), 669.)
                  sfsj3up = readerHFup.eval(sj3flav, jet2sj[0].Eta(), 669.)
                  sfsj3down = readerHFdown.eval(sj3flav, jet2sj[0].Eta(), 669.)
		else:
		 if(jet2sj[0].Pt()<670. and abs(jet2sj[0].Eta())<2.4 ) :
                  sfsj3 = readerLF.eval(sj3flav, jet2sj[0].Eta(), jet2sj[0].Pt())  # jet flavor, eta, pt
                  sfsj3up = readerLFup.eval(sj3flav, jet2sj[0].Eta(), jet2sj[0].Pt())
                  sfsj3down = readerLFdown.eval(sj3flav, jet2sj[0].Eta(), jet2sj[0].Pt())
                 elif abs(jet2sj[0].Eta())>2.4:
                  sfsj3 = readerLF.eval(sj3flav, 2.399, jet2sj[0].Pt())
                  sfsj3up = readerLFup.eval(sj3flav, 2.399, jet2sj[0].Pt())
                  sfsj3down = readerLFdown.eval(sj3flav, 2.399, jet2sj[0].Pt())
                 elif jet2sj[0].Pt()>670.:
                  sfsj3 = readerLF.eval(sj3flav, jet2sj[0].Eta(), 669.)
                  sfsj3up = readerLFup.eval(sj3flav, jet2sj[0].Eta(), 669.)
                  sfsj3down = readerLFdown.eval(sj3flav, jet2sj[0].Eta(), 669.)
		if not isL2:
                 if(jet2sj[1].Pt()<670. and abs(jet2sj[1].Eta())<2.4 ) :
                  sfsj4 = readerHF.eval(sj4flav, jet2sj[1].Eta(), jet2sj[1].Pt())  # jet flavor, eta, pt
                  sfsj4up = readerHFup.eval(sj4flav, jet2sj[1].Eta(), jet2sj[1].Pt())
                  sfsj4down = readerHFdown.eval(sj4flav, jet2sj[1].Eta(), jet2sj[1].Pt())
                 elif abs(jet2sj[1].Eta())>2.4:
                  sfsj4 = readerHF.eval(sj4flav, 2.399, jet2sj[1].Pt())
                  sfsj4up = readerHFup.eval(sj4flav, 2.399, jet2sj[1].Pt())
                  sfsj4down = readerHFdown.eval(sj4flav, 2.399, jet2sj[1].Pt())
                 elif jet2sj[1].Pt()>670.:
                  sfsj4 = readerHF.eval(sj4flav, jet2sj[1].Eta(), 669.)
                  sfsj4up = readerHFup.eval(sj4flav, jet2sj[1].Eta(), 669.)
                  sfsj4down = readerHFdown.eval(sj4flav, jet2sj[1].Eta(), 669.)
		else:
		 if(jet2sj[1].Pt()<670. and abs(jet2sj[1].Eta())<2.4 ) :
                  sfsj4 = readerLF.eval(sj4flav, jet2sj[1].Eta(), jet2sj[1].Pt())  # jet flavor, eta, pt
                  sfsj4up = readerLFup.eval(sj4flav, jet2sj[1].Eta(), jet2sj[1].Pt())
                  sfsj4down = readerLFdown.eval(sj4flav, jet2sj[1].Eta(), jet2sj[1].Pt())
                 elif abs(jet2sj[1].Eta())>2.4:
                  sfsj4 = readerLF.eval(sj4flav, 2.399, jet2sj[1].Pt())
                  sfsj4up = readerLFup.eval(sj4flav, 2.399, jet2sj[1].Pt())
                  sfsj4down = readerLFdown.eval(sj4flav, 2.399, jet2sj[1].Pt())
                 elif jet2sj[1].Pt()>670.:
                  sfsj4 = readerLF.eval(sj4flav, jet2sj[1].Eta(), 669.)
                  sfsj4up = readerLFup.eval(sj4flav, jet2sj[1].Eta(), 669.)
                  sfsj4down = readerLFdown.eval(sj4flav, jet2sj[1].Eta(), 669.)
            
                #sfsj2 = reader.eval(sj2flav, jet1sj[1].Eta(), jet1sj[1].Pt())
#                print "SF " + str(sfsj1) + " for flavor " + str(sj1flav) + " for eta " + str(jet1sj[0].Eta()) + " for pt " + str(jet1sj[0].Pt())

        #for min subjet csv
#	for j in range(len(jet1sjcsv)):
#     	    if jet1sjcsv[j] < jet1mscsv[0]:
#		    jet1mscsv[0] = jet1sjcsv[j]
#	for j in range(len(jet2sjcsv)):
#	    if jet2sjcsv[j] < jet2mscsv[0]:
#		    jet2mscsv[0] = jet2sjcsv[j]
	
	#filling bbtag
	jet1bbtag[0] = jet_bbtag[idxH1] 
	jet2bbtag[0] = jet_bbtag[idxH2]
	
        #writing variables to the tree    
	jet1pt[0] = jets[idxH1].Pt()
	jet2pt[0] = jets[idxH2].Pt()
	jet1eta[0] = jets[idxH1].Eta()
	jet2eta[0] = jets[idxH2].Eta()
	etadiff[0] = abs(jets[idxH1].Eta() - jets[idxH2].Eta())
	dijetmass[0] = (jets[idxH1] + jets[idxH2]).M()
	dijetmass_corr[0] = (jets[idxH1] + jets[idxH2]).M() - (jet1pmass[0]-125)-(jet2pmass[0]-125)
        dijetmass_corr_punc[0] = (jets[idxH1] + jets[idxH2]).M() - (jet1pmassunc[0]-125)-(jet2pmassunc[0]-125)
        if options.isMC:
            puWeights[0]= puweight
            puWeightsUp[0] = puweightUp
            puWeightsDown[0] = puweightDown
            nTrueInt[0] = nTInt 
            xsec[0] = float(options.xsec)
        json[0] = JSON
        evt[0] = EVT
        vtype[0] = vType
        if Data:
            isData[0] = 1
        else:
            isData[0] = 0

        #handling hbb tagger SFs
        sf1 = -1
        sf2 = -1
        sf1change = 1000000
        sf2change = 1000000

        if jet1pt[0] < 400:
            sf1 = 0.929
            sf1change = 0.078
        elif jet1pt[0] >= 400 and jet1pt[0] < 500:
            sf1 = 0.999
            sf1change = 0.126
        elif jet1pt[0] >= 500 and jet1pt[0] < 600:
            sf1 = 0.933
            sf1change = 0.195
        elif jet1pt[0] >= 600:
            sf1 = 1.048
            sf1change = 0.215

        if jet2pt[0] < 400:
            sf2 = 0.929
            sf2change = 0.078
        elif jet2pt[0] >= 400 and jet2pt[0] < 500:
            sf2 = 0.999
            sf2change = 0.126
        elif jet2pt[0] >= 500 and jet2pt[0] < 600:
            sf2 = 0.933
            sf2change = 0.195
        elif jet2pt[0] >= 600:
            sf2 = 1.048
            sf2change = 0.215
        
        SF[0] = sf1*sf2
        SFup[0] = sf1*(1+sf1change)*sf2*(1+sf2change)
        SFdown[0] = sf1*(1-sf1change)*sf2*(1-sf2change)

	SF4sj[0] = -1
	SF4sjUp[0] = -1
	SF4sjDown[0] = -1	
	if n1sj >1 and n2sj>1:
	 SF4sj[0] = sfsj1*sfsj2*sfsj3*sfsj4
	 SF4sjUp[0] = sfsj1up*sfsj2up*sfsj3up*sfsj4up 	
	 SF4sjDown[0] = sfsj1down*sfsj2down*sfsj3down*sfsj4down  
	
	SF3sj[0] =-1.
	SF3sjUp[0] =-1.
	SF3sjDown[0] =-1.
	#3b-tag category
	#SF=[((1-SF1e1))/(1-e1)]*SF2*SF3*SF4
	#e1 estimated in HH signal sample to be 
	'''	
	if (jet1s1csv[0] >0.460 and jet2s1csv[0] >0.460  and jet1s2csv[0] >0.460 and jet2s2csv[0] < 0.460 ) or (jet1s1csv[0] >0.460 and jet2s1csv[0] >0.460  and jet1s2csv[0] <0.460 and jet2s2csv[0] > 0.460 ) or (jet1s1csv[0] >0.460 and jet2s1csv[0] <0.460  and jet1s2csv[0] >0.460 and jet2s2csv[0] > 0.460 ) or (jet1s1csv[0] <0.460 and jet2s1csv[0] > 0.460  and jet1s2csv[0] >0.460 and jet2s2csv[0] > 0.460 ):
	 if n1sj >1 and n2sj>1:
	  if(jet1s1csv[0] >0.460) :
            SF3sj[0] *= sfsj1
            SF3sjUp[0] *= sfsj1up
            SF3sjDown[0] *= sfsj1down
	  else : 
	    SF3sj[0] *= (1-sfsj1*btagging_efficiency_medium(jet1sj[0].Pt()))/(1-btagging_efficiency_medium(jet1sj[0].Pt()))
            SF3sjUp[0] *= (1-sfsj1up*btagging_efficiency_medium(jet1sj[0].Pt()))/(1-btagging_efficiency_medium(jet1sj[0].Pt()))
            SF3sjDown[0] *= (1-sfsj1down*btagging_efficiency_medium(jet1sj[0].Pt()))/(1-btagging_efficiency_medium(jet1sj[0].Pt()))
          if(jet1s2csv[0] >0.460) :
            SF3sj[0] *= sfsj2
            SF3sjUp[0] *= sfsj2up
            SF3sjDown[0] *= sfsj2down
          else :
            SF3sj[0] *= (1-sfsj2*btagging_efficiency_medium(jet1sj[1].Pt()))/(1-btagging_efficiency_medium(jet1sj[1].Pt()))
            SF3sjUp[0] *= (1-sfsj2up*btagging_efficiency_medium(jet1sj[1].Pt()))/(1-btagging_efficiency_medium(jet1sj[1].Pt()))
            SF3sjDown[0] *= (1-sfsj2down*btagging_efficiency_medium(jet1sj[1].Pt()))/(1-btagging_efficiency_medium(jet1sj[1].Pt()))

          if(jet2s1csv[0] >0.460) :
            SF3sj[0] *= sfsj3
            SF3sjUp[0] *= sfsj3up
            SF3sjDown[0] *= sfsj3down
          else:
            SF3sj[0] *= (1-sfsj3*btagging_efficiency_medium(jet2sj[0].Pt()))/(1-btagging_efficiency_medium(jet2sj[0].Pt()))
            SF3sjUp[0] *= (1-sfsj3up*btagging_efficiency_medium(jet2sj[0].Pt()))/(1-btagging_efficiency_medium(jet2sj[0].Pt()))
            SF3sjDown[0] *= (1-sfsj3down*btagging_efficiency_medium(jet2sj[0].Pt()))/(1-btagging_efficiency_medium(jet2sj[0].Pt()))
          if(jet2s2csv[0] >0.460) :
            SF3sj[0] *= sfsj4
            SF3sjUp[0] *= sfsj4up
            SF3sjDown[0] *= sfsj4down
          else:
            SF3sj[0] *= (1-sfsj4*btagging_efficiency_medium(jet2sj[1].Pt()))/(1-btagging_efficiency_medium(jet2sj[1].Pt()))
            SF3sjUp[0] *= (1-sfsj4up*btagging_efficiency_medium(jet2sj[1].Pt()))/(1-btagging_efficiency_medium(jet2sj[1].Pt()))
            SF3sjDown[0] *= (1-sfsj4down*btagging_efficiency_medium(jet2sj[1].Pt()))/(1-btagging_efficiency_medium(jet2sj[1].Pt()))
	
	if SF3sj[0] <0. : SF3sj[0] = -SF3sj[0]
	if SF3sjUp[0] <0. : SF3sjUp[0] = -SF3sjUp[0]
	if SF3sjDown[0] <0. : SF3sjDown[0] = -SF3sjDown[0]	
  	'''
	#filling the tree
        myTree.Fill()

	#filling error values for each object
	jet1pt[0] = -100.0
	jet2pt[0] = -100.0
	jet1eta[0] = -100.0
	jet2eta[0] = -100.0
	etadiff[0] = -100.0
	dijetmass[0] = -100.0
	dijetmass_corr[0]=-100.0
	jet1pmass[0] = -100.0
	jet2pmass[0] = -100.0
	jet1tau21[0] = -100.0
	jet2tau21[0] = -100.0
#	jet1mscsv[0] = -100.0
#	jet2mscsv[0] = -100.0
	jet1bbtag[0] = -100.0
	jet2bbtag[0] = -100.0
	triggerpassbb[0] = -100.0
#	PUWeight[0]= -100.0
	
    
    f1.Close()

print "OK"

f.cd()
f.Write()
f.Close()




