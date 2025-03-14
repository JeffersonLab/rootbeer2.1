#!/bin/cshrc
#setup rootbeer environment
if (!($?ROOTBEER)) then
        echo Error: '$ROOTBEER' must be defined for rootbeer
        exit
endif

setenv ROOTBEER_SCRIPTS     $ROOTBEER/scripts

#workout the osname
setenv OSCLAS       `${ROOTBEER_SCRIPTS}/uname_clas`
setenv OS_NAME      `${ROOTBEER_SCRIPTS}/uname_clas`

setenv ROOTBEER_LIB     $ROOTBEER/lib/${OS_NAME}
setenv ROOTBEER_SLIB    $ROOTBEER/slib/${OS_NAME}
setenv ROOTBEER_BIN     $ROOTBEER/bin/${OS_NAME}
setenv ROOTBEER_OBJ     $ROOTBEER/obj/${OS_NAME}

#make directories if they don't exist
if (! -d $ROOTBEER_LIB) mkdir -p $ROOTBEER_LIB
if (! -d $ROOTBEER_SLIB) mkdir -p $ROOTBEER_SLIB
if (! -d $ROOTBEER_BIN) mkdir -p $ROOTBEER_BIN
if (! -d $ROOTBEER_OBJ) mkdir -p $ROOTBEER_OBJ

setenv PATH "${ROOTBEER_BIN}:${PATH}"
  if ($?LD_LIBRARY_PATH) then
        setenv LD_LIBRARY_PATH "${ROOTBEER_SLIB}:${LD_LIBRARY_PATH}"
  else
        setenv LD_LIBRARY_PATH "${ROOTBEER_SLIB}"
  endif

alias rootbeer "cd $ROOTBEER;root -l RootBeerSetup.cxx"
