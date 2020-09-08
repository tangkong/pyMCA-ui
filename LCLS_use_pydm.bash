#!/bin/bash
LCLS_ROOT=${LCLS_ROOT:-/afs/slac/g/lcls}
PACKAGE_TOP=${PACKAGE_TOP:-${PACKAGE_SITE_TOP:-${LCLS_ROOT}/package}}

# Use a reasonable version of Qt 5
source $PACKAGE_TOP/Qt/Qt-5.9.5/linux-x86_64/use_qt_5.9.5.bash

# Use a 64-bit version of EPICS
source /afs/slac/g/lcls/epics/setup/epicsenv-cur.bash

#Override the python path to use 2.7.9 x86-64:
export PATH=$PACKAGE_TOP/python/python2.7.9_bak/linux-x86_64/bin/:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PACKAGE_TOP/python/python2.7.9_bak/linux-x86_64/lib

export PYDM_EPICS_LIB=pyepics
export PYDM_TOOLS_PATH=$TOOLS/pydm/tools

export PYQTDESIGNERPATH=$PACKAGE_TOP/pydm/current:$PYQTDESIGNERPATH

#PyDM EDM Button Widget setup
export PYQTDESIGNERPATH=$TOOLS/pydm/widgets/edmbutton/edmbutton:$PYQTDESIGNERPATH
export PATH=$PATH:$PACKAGE_TOP/wmctrl/current/linux-x86_64/bin


export PYDM_DEFAULT_PROTOCOL=ca

#Remove old QT definitions that are added by the default environment.
unset QT_VER
unset QTDIR
unset QTINC
unset QTLIB
unset QT_PLUGIN_PATH
