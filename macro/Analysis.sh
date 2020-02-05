#!/bin/tcsh
@ run = $2
echo $run
source /afs/cern.ch/work/c/ctuysuz/private/he_flux/fix.sh
source /afs/cern.ch/work/c/ctuysuz/private/he_flux/amsvar.csh
source /afs/cern.ch/work/c/ctuysuz/private/he_flux/macro/Run.sh $run
