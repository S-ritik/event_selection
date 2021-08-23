#!/bin/sh
#cd /home/deroy/t3store3/condor_jobs/CMSSW_10_5_0/src/SameWeight/
#cd -
#eval `cmsenv`
export X509_USER_PROXY=/home/rsaxena/x509userproxy/x509up_u56665
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
export SCRAM_ARCH=slc6_amd64_gcc700
source $VO_CMS_SW_DIR/cmsset_default.sh 

export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH

cd /home/rsaxena/t3store3/CMSSW_10_2_27/src//ee/single_t/
eval `scramv1 runtime -sh`

cp Summer19UL18_ST_t4f_antitop.log test.log
root runAll.C -l -b -q &> out1.txt
cp test_output.root EE_Summer19UL18_ST_t4f_antitop_output.root

cp Summer19UL18_ST_t4f_top.log test.log
root runAll.C -l -b -q &> out2.txt
cp test_output.root EE_Summer19UL18_ST_t4f_top_output.root

cp Summer20UL18_ST_s_lep.log test.log
root runAll.C -l -b -q &> out3.txt
cp test_output.root EE_Summer20UL18_ST_s_lep_output.root

cp Summer20UL18_ST_tW_antitop.log test.log
root runAll.C -l -b -q &> out4.txt
cp test_output.root EE_Summer20UL18_ST_tW_antitop_output.root

cp Summer20UL18_ST_tW_top.log test.log
root runAll.C -l -b -q &> out5.txt
cp test_output.root EE_Summer20UL18_ST_tW_top_output.root
