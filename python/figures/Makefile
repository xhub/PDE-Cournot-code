PYTHON   := python
ROOT_DIR := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
HDF5_FILES := $(wildcard *.hdf5)
# options for figures
OPTS := 
OPTS_CMP ?= --label-reduced-cost

.PHONY: econ_fig physical_qty_loop all clean 

all: econ_fig physical_qty_loop

econ_fig:
	$(PYTHON) $(ROOT_DIR)/fig-econ-vals.py $(OPTS) $(HDF5_FILES)

physical_qty_loop:
	$(PYTHON) $(ROOT_DIR)/figures-cmp.py $(OPTS) $(OPTS_CMP) $(HDF5_FILES)

