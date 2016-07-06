#!/bin/sh

python generalTreeAnalyzer_76X.py --pathIn=/eos/uscms/store/group/lpchbb/HeppyNtuples/V21/JetHT/VHBB_HEPPY_V21_JetHT__Run2015C_25ns-16Dec2015-v1/160318_132855/0000/ --outName=Jet_HT_C --trigger=False --jets=True --deta=True --min=0 --max=31 --file=TxtFiles/76XRunC.txt &

for i in `seq 0 20 780`
do
python generalTreeAnalyzer_76X.py --pathIn=/eos/uscms/store/group/lpchbb/HeppyNtuples/V21/JetHT/VHBB_HEPPY_V21_JetHT__Run2015D-16Dec2015-v1/160317_130618/0000/ --outName=Jet_HT_D --trigger=False --jets=True --deta=True --min=$i --max=$((i+20)) --file=TxtFiles/76XRunD.txt &
done

python generalTreeAnalyzer_76X.py --pathIn=/eos/uscms/store/group/lpchbb/HeppyNtuples/V21/JetHT/VHBB_HEPPY_V21_JetHT__Run2015D-16Dec2015-v1/160317_130618/0000/ --outName=Jet_HT_D --trigger=False --jets=True --deta=True --min=800 --max=804 --file=TxtFiles/76XRunD.txt &

