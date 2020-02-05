#!/bin/tcsh
@ run1 = $2
@ run2 = $run1 + 10000
@ run3 = $run2 + 10000
@ run4 = $run3 + 10000
source /afs/cern.ch/work/c/ctuysuz/private/he_flux/fix.sh
source /afs/cern.ch/work/c/ctuysuz/private/he_flux/amsvar.csh
echo $run1
source /afs/cern.ch/work/c/ctuysuz/private/he_flux/macro/Run.sh $run1
echo $run2
source /afs/cern.ch/work/c/ctuysuz/private/he_flux/macro/Run.sh $run2
echo $run3
source /afs/cern.ch/work/c/ctuysuz/private/he_flux/macro/Run.sh $run3
echo $run4
source /afs/cern.ch/work/c/ctuysuz/private/he_flux/macro/Run.sh $run4
