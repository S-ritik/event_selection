#!/bin/sh
#cd /home/deroy/t3store3/condor_jobs/CMSSW_10_5_0/src/SameWeight/
#cd -
#eval `cmsenv`
export X509_USER_PROXY=/home/rsaxena/x509userproxy/x509up_u56665
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export SCRAM_ARCH=slc6_amd64_gcc700
source $VO_CMS_SW_DIR/cmsset_default.sh 

export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH

cd /home/rsaxena/t3store3/CMSSW_10_2_27/src/emu/ttbar_signal/
eval `scramv1 runtime -sh`

cp Summer20UL18_TTBar_DiLeptonic.log test.log
root runAll.C -l -b -q &> out1.txt
cp test_output.root EMu_Summer20UL18_TTBar_DiLeptonic_signal_output.root

cp Summer20UL18_TTBar_Hadronic.log test.log
root runAll.C -l -b -q &> out2.txt
cp test_output.root EMu_Summer20UL18_TTBar_Hadronic_signal_output.root

cp Summer20UL18_TTBar_SemiLeptonic.log test.log
root runAll.C -l -b -q &> out3.txt
cp test_output.root EMu_Summer20UL18_TTBar_SemiLeptonic_signal_output.root
