import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import scipy, copy


import Alphabet_Header
from Alphabet_Header import *
import Plotting_Header
from Plotting_Header import *

#gSystem.Load("PU_C.so")
'''
gInterpreter.ProcessLine( 'TFile *f = TFile::Open("trigger_objects.root")')
gInterpreter.ProcessLine( 'TH1F *histo_efficiency = (TH1F*)f->Get("histo_efficiency")')
gInterpreter.ProcessLine( 'TH1F *histo_efficiency_lower = (TH1F*)f->Get("histo_efficiency_lower")')
gInterpreter.ProcessLine( 'TH1F *histo_efficiency_upper = (TH1F*)f->Get("histo_efficiency_upper")')
gInterpreter.ProcessLine(".x trigger_function.cxx")
'''
#background1 = TFile("Hbb_outputup.root")
#qcd_trigger_up=copy.copy(background1.Get("QCD"))
#background2 = TFile("Hbb_outputlow.root")
#qcd_trigger_low=copy.copy(background2.Get("QCD"))
#qcd_trigger_up.SetName("QCD_up")
#qcd_trigger_low.SetName("QCD_low")
#background1.Close()
#background2.Close()

mass=[900,1000,1200,1400,1600,1800,2000,2500,3000]
VAR = "dijetmass_corr"

#variable bin from dijet analysis 788 838 
binBoundaries = [800, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687,
        1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019 ]#, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509,


vartitle = "m_{X} (GeV)"
sigregcutNoTrig = "(dijetmass_corr>800&jet1tau21<0.6&jet2tau21<0.6&jet2pmass<135&jet2pmass>105&&(jet1pmass<135&jet1pmass>105)&jet2bbtag>0.6 &jet1bbtag>0.6&&jet1ID>0.&&jet2ID>0.)"
sigregcut = "(triggerpassbb>0.&dijetmass_corr>800&jet1tau21<0.6&jet2tau21<0.6&jet2pmass<135&jet2pmass>105&&(jet1pmass<135&jet1pmass>105)&jet2bbtag>0.6 &jet1bbtag>0.6&&jet1ID>0.&&jet2ID>0.)"
#(dijetmass_corr>800&jet1tau21<0.6&jet2tau21<0.6&(jet2pmass<130&jet2pmass>90)&(jet1pmass<130&jet1pmass>100)&(jet1bbtag>0.6&jet2bbtag>0.6))"
lumi =2691.
background = TFile("Hbb_Unblind.root")
UD = ['Up','Down']
SF_tau21=0.979
#SF_bbtag=0.95



for m in mass:

	output_file = TFile("datacard/hh_mX_%s_13TeV.root"%(m),"RECREATE")
	hh=output_file.mkdir("hh")
	hh.cd()

	Signal_mX = TH1F("Signal_mX_%s"%(m), "", len(binBoundaries)-1, array('d',binBoundaries))
	Signal_mX_trig_up = TH1F("Signal_mX_%s_CMS_eff_trigUp"%(m), "", len(binBoundaries)-1, array('d',binBoundaries))
	Signal_mX_trig_down = TH1F("Signal_mX_%s_CMS_eff_trigDown"%(m), "", len(binBoundaries)-1, array('d',binBoundaries))
	Signal_mX_btag_up = TH1F("Signal_mX_%s_CMS_eff_btagUp"%(m), "", len(binBoundaries)-1, array('d',binBoundaries))
        Signal_mX_btag_down = TH1F("Signal_mX_%s_CMS_eff_btagDown"%(m), "", len(binBoundaries)-1, array('d',binBoundaries))
	Signal_mX_pu_up = TH1F("Signal_mX_%s_CMS_eff_puUp"%(m), "", len(binBoundaries)-1, array('d',binBoundaries))
 	Signal_mX_pu_down = TH1F("Signal_mX_%s_CMS_eff_puDown"%(m), "", len(binBoundaries)-1, array('d',binBoundaries))
	Signal_mX_FJEC_Up = TH1F("Signal_mX_%s_CMS_eff_JECUp"%(m), "", len(binBoundaries)-1, array('d',binBoundaries))
	Signal_mX_FJEC_Down = TH1F("Signal_mX_%s_CMS_eff_JECDown"%(m), "", len(binBoundaries)-1, array('d',binBoundaries))
	Signal_mX_FJER_Up = TH1F("Signal_mX_%s_CMS_eff_JERUp"%(m), "", len(binBoundaries)-1, array('d',binBoundaries))
        Signal_mX_FJER_Down = TH1F("Signal_mX_%s_CMS_eff_JERDown"%(m), "", len(binBoundaries)-1, array('d',binBoundaries))
	Signal_mX_MJEC_Up = TH1F("Signal_mX_%s_CMS_eff_massJECUp"%(m), "", len(binBoundaries)-1, array('d',binBoundaries))
        Signal_mX_MJEC_Down = TH1F("Signal_mX_%s_CMS_eff_massJECDown"%(m), "", len(binBoundaries)-1, array('d',binBoundaries))

	print(m)

	signal_file= TFile("/uscms_data/d3/cvernier/test-alice/HbbAnalysis/BG_%s_76X_boost_None.root"%(m))
	#signal_file= TFile("/uscms_data/d3/cvernier/test-alice/HbbAnalysis/Radion_%s_76X_boost_None.root"%(m))
	#signal_file= TFile("/uscms_data/d3/asady1/HHcode/CMSSW_7_4_2/src/Analysis/HbbAnalysis/BG_%s_76X_boost_0.root"%(m))
	bbj = signal_file.Get("bbj")
	generatedEvents =bbj.GetEntries()
	print(generatedEvents)
	tree = signal_file.Get("myTree") 
	writeplot(tree, Signal_mX_pu_up, VAR, sigregcut, "puWeightsUp*SF")
	writeplot(tree, Signal_mX_pu_down, VAR, sigregcut, "puWeightsDown*SF")
	writeplot(tree, Signal_mX, VAR, sigregcut, "puWeights*SF")#(trigger_function(int(round(htJet40eta3)))*weight2(nTrueInt))")
        writeplot(tree, Signal_mX_btag_up, VAR, sigregcut, "puWeights*SFup")
	writeplot(tree, Signal_mX_btag_down, VAR, sigregcut, "puWeights*SFdown")
        writeplot(tree, Signal_mX_trig_up, VAR, sigregcutNoTrig, "trigWeightUp*puWeights*SF")
	writeplot(tree, Signal_mX_trig_down, VAR, sigregcutNoTrig, "trigWeightDown*puWeights*SF")


	print(Signal_mX.Integral())

	btaglnN= 1.+ abs(Signal_mX_btag_up.GetSumOfWeights()-Signal_mX_btag_down.GetSumOfWeights())/(2.*Signal_mX_btag_up.GetSumOfWeights())
	pulnN= 1.+ abs(Signal_mX_pu_up.GetSumOfWeights()-Signal_mX_pu_down.GetSumOfWeights())/(2.*Signal_mX.GetSumOfWeights())
	#triglnN= 1.+ abs(Signal_mX_trig_up.GetSumOfWeights()-Signal_mX_trig_down.GetSumOfWeights())/(2.*Signal_mX.GetSumOfWeights())

	Signal_mX.Scale(lumi*SF_tau21*SF_tau21*0.01/generatedEvents)
	
 	Signal_mX_btag_up.Scale(lumi*SF_tau21*SF_tau21*0.01/generatedEvents)
	Signal_mX_btag_down.Scale(lumi*SF_tau21*SF_tau21*0.01/generatedEvents)
	Signal_mX_pu_up.Scale(lumi*SF_tau21*SF_tau21*0.01/generatedEvents)
        Signal_mX_pu_down.Scale(lumi*SF_tau21*SF_tau21*0.01/generatedEvents)
	Signal_mX_trig_up.Scale(0.01*lumi*SF_tau21*SF_tau21/generatedEvents)
	Signal_mX_trig_down.Scale(0.01*lumi*SF_tau21*SF_tau21/generatedEvents)
        
	if not (m==1200 or m==3000 or m==1800):	
	 signal_file_FJEC_Up= TFile("/uscms_data/d3/cvernier/test-alice/HbbAnalysis/BG_%s_76X_boost_FJEC_Up.root"%(m))
	 tree_FJEC_Up = signal_file_FJEC_Up.Get("myTree")
	 writeplot(tree_FJEC_Up, Signal_mX_FJEC_Up, VAR, sigregcut, "puWeights*SF")
	 Signal_mX_FJEC_Up.Scale(lumi*SF_tau21*SF_tau21*0.01/generatedEvents)
		
	
	 signal_file_FJEC_Down= TFile("/uscms_data/d3/cvernier/test-alice/HbbAnalysis/BG_%s_76X_boost_FJEC_Down.root"%(m))
         tree_FJEC_Down = signal_file_FJEC_Down.Get("myTree")
         writeplot(tree_FJEC_Down, Signal_mX_FJEC_Down, VAR, sigregcut, "puWeights*SF")
         Signal_mX_FJEC_Down.Scale(lumi*SF_tau21*SF_tau21*0.01/generatedEvents)
	
	 signal_file_FJER_Up= TFile("/uscms_data/d3/cvernier/test-alice/HbbAnalysis/BG_%s_76X_boost_FJER_Up.root"%(m))
         tree_FJER_Up = signal_file_FJER_Up.Get("myTree")
         writeplot(tree_FJER_Up, Signal_mX_FJER_Up, VAR, sigregcut, "puWeights*SF")
         Signal_mX_FJER_Up.Scale(lumi*SF_tau21*SF_tau21*0.01/generatedEvents)

         signal_file_FJER_Down= TFile("/uscms_data/d3/cvernier/test-alice/HbbAnalysis/BG_%s_76X_boost_FJER_Down.root"%(m))
         tree_FJER_Down = signal_file_FJER_Down.Get("myTree")
         writeplot(tree_FJER_Down, Signal_mX_FJER_Down, VAR, sigregcut, "puWeights*SF")
         Signal_mX_FJER_Down.Scale(lumi*SF_tau21*SF_tau21*0.01/generatedEvents)

	 signal_file_MJEC_Up= TFile("/uscms_data/d3/cvernier/test-alice/HbbAnalysis/BG_%s_76X_boost_MJEC_Up.root"%(m))
         tree_MJEC_Up = signal_file_MJEC_Up.Get("myTree")
         writeplot(tree_MJEC_Up, Signal_mX_MJEC_Up, VAR, sigregcut, "puWeights*SF")
         Signal_mX_MJEC_Up.Scale(lumi*SF_tau21*SF_tau21*0.01/generatedEvents)

         signal_file_MJEC_Down= TFile("/uscms_data/d3/cvernier/test-alice/HbbAnalysis/BG_%s_76X_boost_MJEC_Down.root"%(m))
         tree_MJEC_Down = signal_file_MJEC_Down.Get("myTree")
         writeplot(tree_MJEC_Down, Signal_mX_MJEC_Down, VAR, sigregcut, "puWeights*SF")
         Signal_mX_MJEC_Down.Scale(lumi*SF_tau21*SF_tau21*0.01/generatedEvents)
	
         MJEClnN= 1.+ abs(Signal_mX_MJEC_Up.GetSumOfWeights()-Signal_mX_MJEC_Down.GetSumOfWeights())/(2.*Signal_mX.GetSumOfWeights())
	 FJEClnN= 1.+ abs(Signal_mX_FJEC_Up.GetSumOfWeights()-Signal_mX_FJEC_Down.GetSumOfWeights())/(2.*Signal_mX.GetSumOfWeights())
	 FJERlnN= 1.+ abs(Signal_mX_FJER_Up.GetSumOfWeights()-Signal_mX_FJER_Down.GetSumOfWeights())/(2.*Signal_mX.GetSumOfWeights())
	else :
	

		MJEClnN= 1.02
		FJEClnN= 1.02
		FJERlnN= 1.02
	
	print("%s %s %s")%(MJEClnN,FJEClnN,FJERlnN)

        signal_integral = Signal_mX.Integral()
	print(signal_integral) 
        background.cd() 	
	qcd_integral = EST.Integral()

	qcd =EST
	qcd_antitag = EST_Antitag
	qcd_up = EST_CMS_scale_13TeVUp
	qcd_down = EST_CMS_scale_13TeVDown
	data = data_obs
        data_integral = data_obs.Integral() 
	output_file.cd()
        hh.cd()
	qcd_stat_up =TH1F("qcd_stat_up","",len(binBoundaries)-1, array('d',binBoundaries))
        qcd_stat_down =TH1F("qcd_stat_down","",len(binBoundaries)-1, array('d',binBoundaries))
	
	for bin in range(0,len(binBoundaries)-1):
            for Q in UD:
                qcd_syst =TH1F("%s_bin%s%s"%("EST_CMS_stat_13TeV",bin,Q),"",len(binBoundaries)-1, array('d',binBoundaries))
 		bin_stat = qcd.GetBinContent(bin+1)
		for bin1 in range(0,len(binBoundaries)-1):
			bin_stat1 = qcd.GetBinContent(bin1+1)
			qcd_syst.SetBinContent(bin1+1,bin_stat1)
		#if bin_stat==0 :	
		#	bin_stat = 1.5
		bin_at = qcd_antitag.GetBinContent(bin+1)
		if bin_at < 1 and bin_at >0:  
			bin_at=1.
                if Q == 'Up':
			if bin_at >0 :
			       qcd_stat_up.SetBinContent(bin+1,bin_stat+qcd_antitag.GetBinError(bin+1)/bin_at*bin_stat)	
                               qcd_syst.SetBinContent(bin+1,bin_stat+qcd_antitag.GetBinError(bin+1)/bin_at*bin_stat)
			else : 
				qcd_syst.SetBinContent(bin+1,bin_stat)
				qcd_stat_up.SetBinContent(bin+1,bin_stat)
				
                if Q == 'Down':
			if bin_at >0 :
				if ( bin_stat-qcd_antitag.GetBinError(bin+1)/bin_at*bin_stat >0 ):
                        		qcd_syst.SetBinContent(bin+1,bin_stat-qcd_antitag.GetBinError(bin+1)/bin_at*bin_stat)
					qcd_stat_down.SetBinContent(bin+1,bin_stat-qcd_antitag.GetBinError(bin+1)/bin_at*bin_stat)
				else :
					qcd_syst.SetBinContent(bin+1, 0.1)
					qcd_stat_down.SetBinContent(bin+1, 0.1)
			else : 	
				qcd_syst.SetBinContent(bin+1,bin_stat)
				qcd_stat_down.SetBinContent(bin+1,bin_stat)
		qcd_syst.Write()
		

	#qcd_trigger_up.Write()
	#qcd_trigger_low.Write()
	qcd.Write()
	qcd_up.Write()
	qcd_down.Write()
	qcd_stat_up.Write()
	qcd_stat_down.Write()
	Signal_mX.Write()
	Signal_mX_btag_up.Write()
	Signal_mX_btag_down.Write()
	Signal_mX_pu_up.Write()
	Signal_mX_pu_down.Write()
	Signal_mX_trig_up.Write()
	Signal_mX_trig_down.Write()
	data.Write()
	hh.Write()
	output_file.Write()
	#output_file.Close()

	

	text_file = open("datacard/hh_mX_%s_13TeV.txt"%(m), "w")


	text_file.write("max    1     number of categories\n")
        text_file.write("jmax   1     number of samples minus one\n")
        text_file.write("kmax    *     number of nuisance parameters\n")
        text_file.write("-------------------------------------------------------------------------------\n")
        text_file.write("shapes * * hh_mX_%s_13TeV.root hh/$PROCESS hh/$PROCESS_$SYSTEMATIC\n"%(m))
        text_file.write("-------------------------------------------------------------------------------\n")
        text_file.write("bin                                            hh4b\n")
        text_file.write("observation                                    %f\n"%(data_integral))
        text_file.write("-------------------------------------------------------------------------------\n")
        text_file.write("bin                                             hh4b            hh4b\n")
        text_file.write("process                                          0      1\n")
        text_file.write("process                                         Signal_mX_%s  EST\n"%(m))
        text_file.write("rate                                            %f  %f\n"%(signal_integral,qcd_integral))
        text_file.write("-------------------------------------------------------------------------------\n")
	text_file.write("lumi_13TeV lnN                          1.027       -\n")	
        text_file.write("CMS_eff_tau21_sf lnN                    1.057       -\n") #(0.028/0.979)*(2) 
        text_file.write("CMS_pileup lnN                    %f       -\n"%(pulnN))  
        text_file.write("CMS_eff_Htag_sf lnN                    1.1       -\n")   
        text_file.write("CMS_JEC lnN 		     %f        -\n"%(FJEClnN)) 	
	text_file.write("CMS_massJEC lnN                 %f        -\n"%(MJEClnN))
	text_file.write("CMS_eff_bbtag_sf lnN                    %f       -\n"%(btaglnN))
        text_file.write("CMS_JER lnN                    %f        -\n"%(FJERlnN))
	#text_file.write("CMS_eff_trig lnN           %f   -\n"%(triglnN))
	text_file.write("CMS_eff_trig shapeN2           1.0   -\n")
        text_file.write("CMS_scale_13TeV shapeN2                           -       1.000\n")
	text_file.write("CMS_PDF_Scales lnN   1.02 -       \n")

	for bin in range(0,len(binBoundaries)-1):
		text_file.write("CMS_stat_13TeV_bin%s shapeN2                           -       1.000\n"%(bin))


	text_file.close()


	qcd_up.SetLineColor(kBlack)
        qcd_down.SetLineColor(kBlack)
	qcd_up.SetLineStyle(2)
	qcd_down.SetLineStyle(2)
  	qcd_stat_up.SetLineColor(kAzure+1)
        qcd_stat_down.SetLineColor(kAzure+1)
        qcd_stat_up.SetLineStyle(2)
        qcd_stat_down.SetLineStyle(2)	
	qcd.SetLineColor(kPink+3)
	qcd.SetLineWidth(3)
	qcd.SetFillColor(0)



	data.SetStats(0)
	data.Sumw2()
	data.SetLineColor(1)
	data.SetFillColor(0)
	data.SetMarkerColor(1)
	data.SetMarkerStyle(20)
	#qcd.GetYaxis().SetTitle("events / "+str((bins[2]-bins[1])/bins[0])+" GeV")
	qcd.GetXaxis().SetTitle(vartitle)

	leg2 = TLegend(0.6,0.6,0.89,0.89)
	leg2.SetLineColor(0)
	leg2.SetFillColor(0)
	leg2.AddEntry(data, "QCD in SR", "PL")
	leg2.AddEntry(qcd, "QCD prediction", "LF")
	leg2.AddEntry(qcd_up, "transfer function uncertainty", "F")
	leg2.AddEntry(qcd_stat_up, "statistical uncertainty", "F")

	print("asd")
	continue
	#FindAndSetMax([qcd,qcd_up,qcd_stat_up,data])
	C3 = TCanvas("C3", "", 800, 600)
	C3.cd()
	print("asd")
	data.Draw("same E0")
	qcd_up.Draw("same")
	qcd_down.Draw("same")
	qcd_stat_up.Draw("same")
        qcd_stat_down.Draw("same")
	if m< 1400 :
		C3.Print("split_unc.pdf")


	#output_file.Close()

