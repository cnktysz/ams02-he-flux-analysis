#!/bin/tcsh
setenv Offline    /afs/cern.ch/ams/Offline
setenv Offlinea   /afs/cern.ch/ams/Offline
setenv Offline    /cvmfs/ams.cern.ch/Offline

#setenv INTEL      $Offline/intel/Compiler/11.1/073
#setenv INTEL      /opt/intel/Compiler/11.1/073
#setenv INTEL      /afs/cern.ch/ams/Offline/intel/Compiler/11.1/073
setenv INTEL /cvmfs/projects.cern.ch/intelsw/psxe/linux/x86_64/2018/compilers_and_libraries_2018.2.199/linux/

setenv ROOTSYSI   $Offline/root/Linux/root-v5-34-9-icc64-slc5
setenv ROOTSYSG   $Offline/root/Linux/root-v5-34-9-gcc64-slc5
setenv ROOTSYS6   $Offline/root/Linux/root-v5-34-9-icc64-slc6
setenv ROOTSYS64  $Offline/root/Linux/root-v5-34-9-icc64.14-slc6
setenv ROOTSYS604 $Offline/root/Linux/root6-04-08-icc16
setenv ROOTSYS    $ROOTSYS6

setenv XRDLIB /afs/cern.ch/ams/local2/opt/xrootd-icc64-12
setenv EOSAMS root://eosams.cern.ch///eos/ams/Data/AMS02/2014/ISS.B950/pass6

setenv ARCH linuxx8664icc5.34
setenv AMSSRC $Offline/vdev/
setenv AMSWD  $Offline/vdev/lib/linuxx8664icc5.34/

if `awk '{print $6}' /etc/redhat-release | cut -c 1-1` == "6" then
  setenv INTEL  /afs/cern.ch/ams/opt/intel/composer_xe_2013_sp1.3.174/compiler
  setenv XRDLIB /afs/cern.ch/ams/local2/opt/xrootd-icc64-12
  setenv ANTLIB libntuple_slc6_PG.a
  setenv NTUPLE ntuple_slc6_PG.so
  setenv SLC    6
  #source /afs/cern.ch/project/eos/installation/ams/etc/setup.csh
endif

setenv AMSDataDir $Offline/AMSDataDir

setenv AMSDataDirFC /fc02dat0/Offline/AMSDataDir
setenv RunsDir    $Offline/RunsDir
setenv AMSGeoDir  $Offline/vdev/display/
setenv amsedPG    $Offline/vdev/exe/linuxx8664icc5.34/amsedcPG
setenv amsedcPG   $Offline/vdev/exe/linuxx8664icc5.34/amsedcPG
setenv amsed      $AMSWD/exe/linuxx8664icc5.34/amsedcPG

setenv CASTORSTATIC 1

setenv AMSEDCPG 1
setenv PGTRACK 1
setenv ECALBDT 1
setenv MINUIT2 1
setenv AMSICC  1
setenv AMSP    1

#setenv CLHEP_BASE_DIR /afs/cern.ch/exp/ams/Offline/CLHEP.2.0.4.7.icc64
#setenv G4INSTALL /afs/cern.ch/exp/ams/Offline/geant4.9.3.p02

#setenv CLHEP_BASE_DIR /cvmfs/ams.cern.ch/Offline/CLHEP.2.1.0.1.icc64
#setenv G4INSTALL /cvmfs/ams.cern.ch/Offline/geant4.9.4.p04/

setenv CLHEP_BASE_DIR $Offline/CLHEP.2.1.0.1.icc64
#setenv G4INSTALL /cvmfs/ams.cern.ch/Offline/geant4.9.6.p02/
setenv G4INSTALL /Offline/geant4.9.6.p02

setenv GEANT4NEW 1
#setenv G4AMS   1
setenv G4SYSTEM Linux-icc.slc6
setenv G4USE_STL 1
setenv OGLHOME /usr
setenv G4WORKDIR /tmp
setenv G4LIB $G4INSTALL/lib
setenv G4_NO_VERBOSE 1

alias ag4 "setenv G4AMS 1 ; setenv G4SYSTEM Linux-icc.slc6"


set path = ( $ROOTSYS/bin $INTEL/bin/intel64 $path )


## IMPORT FROM GCC
set gcc_config_version = 5.2.0
set LCGPLAT = x86_64-slc6-gcc52-opt
set LCG_lib_name = lib64
set LCG_arch = x86_64

set LCG_contdir = /cvmfs/sft.cern.ch/lcg/contrib
set LCG_gcc_home = ${LCG_contdir}/gcc/${gcc_config_version}/${LCGPLAT}

setenv PATH ${LCG_gcc_home}/bin:${PATH}
setenv COMPILER_PATH ${LCG_gcc_home}/lib/gcc/${LCG_arch}-unknown-linux-gnu/${gcc_config_version}

if ($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH ${LCG_gcc_home}/${LCG_lib_name}:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH ${LCG_gcc_home}/${LCG_lib_name}
endif
setenv FC `which gfortran`
setenv CXX `which g++`
setenv CC `which gcc`
## END OF IMPORT GCC


## IMPORT FROM ICC CSH
set arch = intel64
set archdir = bin64

setenv INTEL_LOCAL_DIR /cvmfs/projects.cern.ch/intelsw/psxe/linux
source $INTEL_LOCAL_DIR/setup.csh

source $INTEL_LOCAL_DIR/x86_64/2018/parallel_studio_xe_2018.2.046/psxevars.csh $arch
## END OF IMPORT ICC

#setenv LD_LIBRARY_PATH ".:${INTEL}/idb/lib/intel64:${INTEL}/lib/intel64"
setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:${ROOTSYS}/lib" 
setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:$XRDLIB/lib64"
setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:/usr/lib64"
setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:/usr/lib"

#if `awk '{print $6}' /etc/redhat-release` == "6.4" then
  #  setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:/opt/intel/compiler80/lib"
  #endif

limit core 0
limit data 700000
limit stack 32000
