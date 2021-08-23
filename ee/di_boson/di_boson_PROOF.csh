#!/bin/sh
#cd /home/deroy/t3store3/condor_jobs/CMSSW_10_5_0/src/SameWeight/
#cd -
#eval `cmsenv`
export X509_USER_PROXY=/home/rsaxena/x509userproxy/x509up_u56665
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export SCRAM_ARCH=slc6_amd64_gcc700
source $VO_CMS_SW_DIR/cmsset_default.sh 

export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH

cd /home/rsaxena/t3store3/CMSSW_10_2_27/src/ee/di_boson/
eval `scramv1 runtime -sh`

cp Summer20UL18_WW.log test.log
root runAll.C -l -b -q &> out1.txt
cp test_output.root EE_Summer20UL18_WW_output.root

cp Summer20UL18_WZ.log test.log
root runAll.C -l -b -q &> out2.txt
cp test_output.root EE_Summer20UL18_WZ_output.root

cp Summer20UL18_ZZ.log test.log
root runAll.C -l -b -q &> out3.txt
cp test_output.root EE_Summer20UL18_ZZ_output.root
