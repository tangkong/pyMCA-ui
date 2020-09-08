#!/bin/bash

conda activate pydm-env

export PYQTDESIGNERPATH=$PACKAGE_TOP/pydm/pydm-master:$PYQTDESIGNERPATH
export PYQTDESIGNERPATH=$TOOLS/pydm/widgets/edmbutton/edmbutton:$PYQTDESIGNERPATH

export PYDM_EPICS_LIB=pyepics
export PYDM_TOOLS_PATH=$TOOLS/pydm/tools

