#!/bin/sh
#cd /home/deroy/t3store3/condor_jobs/CMSSW_10_5_0/src/SameWeight/
#cd -
#eval `cmsenv`
export X509_USER_PROXY=/home/rsaxena/x509userproxy/x509up_u56665
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export SCRAM_ARCH=slc6_amd64_gcc700
source $VO_CMS_SW_DIR/cmsset_default.sh 

export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH

cd /home/rsaxena/t3store3/CMSSW_10_2_27/src/emu/qcd/
eval `scramv1 runtime -sh`

cp Summer19UL18_QCD_MG_100to200.log test.log
root runAll.C -l -b -q  &> out1.txt
cp test_output.root EMu_Summer19UL18_QCD_MG_100to200_output.root

cp Summer19UL18_QCD_MG_200to300.log test.log
root runAll.C -l -b -q  &> out2.txt
cp test_output.root EMu_Summer19UL18_QCD_MG_200to300_output.root

cp Summer19UL18_QCD_MG_300to500.log test.log
root runAll.C -l -b -q  &> out3.txt
cp test_output.root EMu_Summer19UL18_QCD_MG_300to500_output.root

cp Summer19UL18_QCD_MG_500to700.log test.log
root runAll.C -l -b -q  &> out4.txt
cp test_output.root EMu_Summer19UL18_QCD_MG_500to700_output.root

cp Summer19UL18_QCD_MG_700to1000.log test.log
root runAll.C -l -b -q  &> out5.txt
cp test_output.root EMu_Summer19UL18_QCD_MG_700to1000_output.root

cp Summer19UL18_QCD_MG_1000to1500.log test.log
root runAll.C -l -b -q  &> out6.txt
cp test_output.root EMu_Summer19UL18_QCD_MG_1000to1500_output.root

cp Summer19UL18_QCD_MG_1500to2000.log test.log
root runAll.C -l -b -q  &> out7.txt
cp test_output.root EMu_Summer19UL18_QCD_MG_1500to2000_output.root

cp Summer19UL18_QCD_MG_2000toInf.log test.log
root runAll.C -l -b -q  &> out8.txt
cp test_output.root EMu_Summer19UL18_QCD_MG_2000toInf_output.root
