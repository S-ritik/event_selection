#!/bin/sh
#cd /home/deroy/t3store3/condor_jobs/CMSSW_10_5_0/src/SameWeight/
#cd -
#eval `cmsenv`
export X509_USER_PROXY=/home/rsaxena/x509userproxy/x509up_u56665
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export SCRAM_ARCH=slc6_amd64_gcc700
source $VO_CMS_SW_DIR/cmsset_default.sh 

export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH

cd /home/rsaxena/t3store3/CMSSW_10_2_27/src/emu/dy_jets/
eval `scramv1 runtime -sh`

cp Summer20UL18_DYJetsToLL_M-10to50.log test.log
root runAll.C -l -b -q &> out1.txt
cp test_output.root EMu_Summer20UL18_DYJetsToLL_M-10to50_output.root

cp Summer20UL18_DYJetsToLL_M-50_HT-70to100.log test.log
root runAll.C -l -b -q &> out2.txt
cp test_output.root EMu_Summer20UL18_DYJetsToLL_M-50_HT-70to100_output.root

cp Summer20UL18_DYJetsToLL_M-50_HT-100To200.log test.log
root runAll.C -l -b -q &> out3.txt
cp test_output.root EMu_Summer20UL18_DYJetsToLL_M-50_HT-100To200_output.root

cp Summer20UL18_DYJetsToLL_M-50_HT-200To400.log test.log
root runAll.C -l -b -q &> out4.txt
cp test_output.root EMu_Summer20UL18_DYJetsToLL_M-50_HT-200To400_output.root

cp Summer20UL18_DYJetsToLL_M-50_HT-400To600.log test.log
root runAll.C -l -b -q &> out5.txt
cp test_output.root EMu_Summer20UL18_DYJetsToLL_M-50_HT-400To600_output.root

cp Summer20UL18_DYJetsToLL_M-50_HT-600To800.log test.log
root runAll.C -l -b -q &> out6.txt
cp test_output.root EMu_Summer20UL18_DYJetsToLL_M-50_HT-600To800_output.root

cp Summer20UL18_DYJetsToLL_M-50_HT-800To1200.log test.log
root runAll.C -l -b -q &> out7.txt
cp test_output.root EMu_Summer20UL18_DYJetsToLL_M-50_HT-800To1200_output.root

cp Summer20UL18_DYJetsToLL_M-50_HT-1200To2500.log test.log
root runAll.C -l -b -q &> out8.txt
cp test_output.root EMu_Summer20UL18_DYJetsToLL_M-50_HT-1200To2500_output.root

cp Summer20UL18_DYJetsToLL_M-50_HT-2500toInf.log test.log
root runAll.C -l -b -q &> out9.txt
cp test_output.root EMu_Summer20UL18_DYJetsToLL_M-50_HT-2500toInf_output.root
