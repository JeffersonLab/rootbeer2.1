#!/bin/bash
#setup rootbeer environment
#if (!($?ROOTBEER)) then
#        echo Error: '$ROOTBEER' must be defined for rootbeer
#        exit
#fi

export ROOTBEER_SCRIPTS="${ROOTBEER}/scripts"

#workout the osname
#export OSCLAS=`${ROOTBEER_SCRIPTS}/uname_clas`
#export OS_NAME=`${ROOTBEER_SCRIPTS}/uname_clas`

export ROOTBEER_LIB="${ROOTBEER}/lib/Linux"
export ROOTBEER_SLIB="${ROOTBEER}/slib/Linux"
export ROOTBEER_BIN="${ROOTBEER}/bin/Linux"
export ROOTBEER_OBJ="${ROOTBEER}/obj/Linux"

#make directories if they don't exist
#if (! -d $ROOTBEER_LIB) mkdir -p $ROOTBEER_LIB
#if (! -d $ROOTBEER_SLIB) mkdir -p $ROOTBEER_SLIB
#if (! -d $ROOTBEER_BIN) mkdir -p $ROOTBEER_BIN
#if (! -d $ROOTBEER_OBJ) mkdir -p $ROOTBEER_OBJ

export PATH="${ROOTBEER_BIN}:${PATH}"
  #if ($?LD_LIBRARY_PATH) then
        export LD_LIBRARY_PATH="${ROOTBEER_SLIB}:${LD_LIBRARY_PATH}"
  #else
  #      export LD_LIBRARY_PATH="${ROOTBEER_SLIB}"
  #fi

alias rootbeer="cd $ROOTBEER;root -l RootBeerSetup.cxx"
