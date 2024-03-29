#!/bin/bash -f

echo ">>> Setting Environment Variables for EXTERNALS <<<"

export ARCH="64bit"
export ARCHT="64bit"

export EXTERNDIR=${HOME}/EXTERNALS
export SOFT=/user/a/alaoui/software # external, ahmed dir

echo ""
echo " > EXTERNDIR is set to $EXTERNDIR"
echo " > SOFT is set to      $SOFT"
echo " > ARCH is set to      $ARCH"
echo " > ARCHT is set to     $ARCHT"

################
# cern settings
################

echo ""
echo ">> Setting Environment Variables for CERN"
echo ""

export CERNVER=2005 # hardcoded
export CERN=${SOFT}/cern # ahmed dir
export CERN_LEVEL=$CERNVER
export CERN_ROOT=${CERN}/${CERN_LEVEL}
export CERN_LIB=${CERN_ROOT}/lib
export CERN_BIN=${CERN_ROOT}/bin
export CERN_INC=${CERN_ROOT}/include
export CERNLIBDIR=${CERN_LIB}
export CERN_INCLUDEDIR=${CERN_INC}

if [ ! -d $CERN_LIB ]; then
  echo "CERN Error: $CERN_LIB Not Found"
  return -1
fi

if [ ! -d $CERN_BIN ]; then
  echo "CERN Error: $CERN_BIN Not Found"
  return -1
fi

if [ ! -d $CERN_INC ]; then
  echo "CERN Error: $CERN_INC Not Found"
  return -1
fi

echo "CERN_LEVEL           is set to ${CERN_LEVEL}"
echo "CERN_ROOT            is set to ${CERN_ROOT}"
echo "CERN_LIB             is set to ${CERN_LIB}"
echo "CERN_BIN             is set to ${CERN_BIN}"
echo "CERN_INC             is set to ${CERN_INC}"

if [ -z "${PATH}" ]; then
  PATH=${CERN_BIN}
else
  PATH=${CERN_BIN}:${PATH}
fi

if [ -z "$LIBRARY_PATH" ]; then
  export LIBRARY_PATH="${CERN_LIB}"
else
  export LIBRARY_PATH="${CERN_LIB}:$LIBRARY_PATH"
fi

if [ -z "${LD_LIBRARY_PATH}" ]; then
  LD_LIBRARY_PATH=${CERN_LIB}
else
  LD_LIBRARY_PATH=${CERN_LIB}:${LD_LIBRARY_PATH}
fi

################
# root settings
################

echo ""
echo ">> Setting Environment Variables for ROOT"
echo ""

export ROOTVER=6.10.02
export ROOTSYS=/opt/software/root/${ROOTVER}
export ROOTLIB=${ROOTSYS}/lib
export ROOTBIN=${ROOTSYS}/bin
export ROOTINC=${ROOTSYS}/include

if [ ! -d $ROOTLIB ]; then
  echo "ROOT Error: $ROOTLIB Not Found"
  return -1
fi

if [ ! -d $ROOTBIN ]; then
  echo "ROOT Error: $ROOTBIN Not Found"
  return -1
fi

if [ ! -d $ROOTINC ]; then
  echo "ROOT Error: $ROOTINC Not Found"
  return -1
fi

echo "ROOTVER              is set to ${ROOTVER}"
echo "ROOTSYS              is set to ${ROOTSYS}"
echo "ROOTLIB              is set to ${ROOTLIB}"
echo "ROOTBIN              is set to ${ROOTBIN}"
echo "ROOTINC              is set to ${ROOTINC}"

if [ -z "$PATH" ]; then
  export PATH="$ROOTBIN"
else
  export PATH="$ROOTBIN:$PATH"
fi

if [ -z "$LIBRARY_PATH" ]; then
  export LIBRARY_PATH="$ROOTLIB"
else
  export LIBRARY_PATH="$ROOTLIB:$LIBRARY_PATH"
fi

if [ -z "$LD_LIBRARY_PATH" ]; then
  export LD_LIBRARY_PATH="$ROOTLIB"
else
  export LD_LIBRARY_PATH="$ROOTLIB:$LD_LIBRARY_PATH"
fi

################
# CLAS settings
################

echo ""
echo ">> Setting Environment variables for CLAS"
echo ""

export CLAS_ROOT=${SOFT}/clas_software_ver1 # ahmed dir

export CVS_RSH="/usr/bin/ssh"
export CVSROOT="ahmede@cvs.jlab.org:/group/clas/clas_cvs"
export TOP_DIR=${CLAS_ROOT}
export CLAS_BUILD=${CLAS_ROOT}
export CLAS_PACK=${CLAS_ROOT}/packages
export CLAS_CMS=${CLAS_PACK}/cms
export OSNAME=`${CLAS_CMS}/uname_clas` # hardcoded, idk
export OS_NAME=${OSNAME}
export CLAS_TOOL=${CLAS_PACK}/ClasTool
export CLASTOOL=${CLAS_TOOL}
export CLAS_TOOL_BIN=${CLAS_TOOL}/bin/${OS_NAME}
export CLAS_SCRIPTS=${CLAS_PACK}/scripts
export CLAS_LIB=${CLAS_ROOT}/lib/${OS_NAME}
export CLAS_BIN=${CLAS_ROOT}/bin/${OS_NAME}
export HV_LOCATION=${CLAS_PACK}/Hv
export RECSIS=${CLAS_PACK}
export COBRASYS=${CLAS_PACK}/utilities/cobra
export CLAS_SLIB=${CLAS_ROOT}/slib/${OS_NAME}

export CLAS_PARMS=${HOME}/simulations/parms # this dir
export CLAS_CALDB_HOST=atlasusr.fis.utfsm.cl
export CLAS_CALDB_USER=alaoui
export CLAS_CALDB_DBNAME=clas_calib
export CLAS_CALDB_RUNINDEX="clas_calib.RunIndex"
export RECSIS_RUNTIME=${HOME}/simulations/clsrc/recsis/runtime # this dir

echo "CVS_RSH              is set to ${CVS_RSH}"
echo "CVSROOT              is set to ${CVSROOT}"
echo "CLAS_ROOT            is set to ${CLAS_ROOT}"
echo "TOP_DIR              is set to ${TOP_DIR}"
echo "CLAS_BUILD           is set to ${CLAS_BUILD}"
echo "CLAS_PACK            is set to ${CLAS_PACK}"
echo "CLAS_CMS             is set to ${CLAS_CMS}"
echo "CLASTOOL             is set to ${CLASTOOL}"
echo "CLAS_TOOL            is set to ${CLAS_TOOL}"
echo "CLAS_SCRIPTS         is set to ${CLAS_SCRIPTS}"
echo "OSNAME               is set to ${OSNAME}"
echo "OS_NAME              is set to ${OS_NAME}"
echo "CLAS_LIB             is set to ${CLAS_LIB}"
echo "CLAS_BIN             is set to ${CLAS_BIN}"
echo "HV_LOCATION          is set to ${HV_LOCATION}"
echo "RECSIS               is set to ${RECSIS}"
echo "COBRASYS             is set to ${COBRASYS}"
echo "CLAS_SLIB            is set to ${CLAS_SLIB}"
echo "CLAS_PARMS           is set to ${CLAS_PARMS}"
echo "CLAS_CALDB_HOST      is set to ${CLAS_CALDB_HOST}"
echo "CLAS_CALDB_USER      is set to ${CLAS_CALDB_USER}"
echo "CLAS_CALDB_DBNAME    is set to ${CLAS_CALDB_DBNAME}"
echo "CLAS_CALDB_RUNINDEX  is set to ${CLAS_CALDB_RUNINDEX}"
echo "RECSIS_RUNTIME       is set to ${RECSIS_RUNTIME}"

if [ -z "${PATH}" ]; then
  PATH=${CLAS_BIN}:${CLAS_TOOL_BIN}:${CLAS_SCRIPTS}:${COBRASYS}/bin
else
  PATH=${CLAS_BIN}:${CLAS_TOOL_BIN}:${CLAS_SCRIPTS}:${COBRASYS}/bin:${PATH}
fi

if [ -z "$LIBRARY_PATH" ]; then
  export LIBRARY_PATH=${CLAS_LIB}:${COBRASYS}/lib:${CLAS_SLIB}
else
  export LIBRARY_PATH=${CLAS_LIB}:${COBRASYS}/lib:${CLAS_SLIB}:$LIBRARY_PATH
fi

if [ -z "${LD_LIBRARY_PATH}" ]; then
  LD_LIBRARY_PATH=${CLAS_LIB}:${COBRASYS}/lib:${CLAS_SLIB}
else
  LD_LIBRARY_PATH=${CLAS_LIB}:${COBRASYS}/lib:${CLAS_SLIB}:${LD_LIBRARY_PATH}
fi
