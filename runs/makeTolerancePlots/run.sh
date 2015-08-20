#!/bin/sh

TOTAL_SAMPLES=1000

TOLERANCE_VALUES="1e-5 2.5e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 2.5e-3 5e-3 1e-2"
PATIENT_LIST="12 33 35 36 37 39 41 42 43 44"

OCTAVE="octave -q --path ~/bin/PTPB_mfiles"

ALPHA=0.25
BETA=0.033
RBE_MIN=1.25
RBE_MAX=6

if [ "$1" = "plot" ] ; then
    for PATIENT in $PATIENT_LIST ; do
        FILES=""
        for TOLERANCE in $TOLERANCE_VALUES ; do
            OUTFILE="data_patient${PATIENT}_alpha${ALPHA}_beta${BETA}_RBEmin${RBE_MIN}_RBEmax${RBE_MAX}_tol${TOLERANCE}.mat"
            FILES="$FILES, '$OUTFILE'"
        done
        EPSFILE="data_patient${PATIENT}_alpha${ALPHA}_beta${BETA}_RBEmin${RBE_MIN}_RBEmax${RBE_MAX}.eps"
        $OCTAVE --eval "makeToleranceScanPlot([$TOLERANCE_VALUES] $FILES); print('-deps', '$EPSFILE');"
    done
    exit
fi

CORE_COUNT=`python -c "import multiprocessing; print multiprocessing.cpu_count();"`
N_SAMPLES=`python -c "import math; print int(math.ceil(float($TOTAL_SAMPLES) / float($CORE_COUNT)));"`

for PATIENT in $PATIENT_LIST ; do
    for TOLERANCE in $TOLERANCE_VALUES ; do
        OUTFILE="data_patient${PATIENT}_alpha${ALPHA}_beta${BETA}_RBEmin${RBE_MIN}_RBEmax${RBE_MAX}_tol${TOLERANCE}.mat"
        make -j "$CORE_COUNT" PATIENT="$PATIENT" ALPHA="$ALPHA" BETA="$BETA" RBE_MIN="$RBE_MIN" \
            RBE_MAX="$RBE_MAX" TOLERANCE="$TOLERANCE" ITERATIONS="$CORE_COUNT" N_SAMPLES="$N_SAMPLES" \
            OUTFILE="$OUTFILE" OCTAVE="$OCTAVE"
    done
done
