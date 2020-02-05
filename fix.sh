#!/bin/tcsh
setenv HOME /afs/cern.ch/user/c/ctuysuz
setenv Offline /afs/cern.ch/ams/Offline
##setenv AMSINCLUDE /afs/cern.ch/ams/Offline/vdev

setenv STAGE_HOST castorpublic
setenv RFIO_USE_CASTOR_V2 YES
setenv STAGE_SVCCLASS amsuser
setenv CASTOR_INSTANCE castorpublic

source /afs/cern.ch/exp/ams/Offline/root/Linux/root-v5-34-9-gcc64-slc6/amsvar 

#source /afs/cern.ch/project/eos/installation/ams/etc/setup.csh
#eosforceumount /afs/cern.ch/user/c/ctuysuz/eos
#eosmount /afs/cern.ch/user/c/ctuysuz/eos

setenv AMSWD /afs/cern.ch/ams/Offline/vdev
setenv AMSSRC /afs/cern.ch/ams/Offline/vdev

##setenv AMSWD /afs/cern.ch/user/c/ckonak/vdev
##setenv AMSSRC /afs/cern.ch/user/c/ckonak/vdev
setenv ROOTSYS /cvmfs/ams.cern.ch/Offline/AMSsoft/linux_slc6_gcc64/root534/

setenv LANG C
setenv LC_ALL C
