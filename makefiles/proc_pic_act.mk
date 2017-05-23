#!/bin/bash
subjDir=$(shell pwd)
subjects=$(wildcard $(subjDir)/109[0-9][0-9][0-9])
.PHONY: all

all: $(subjects)

%: FORCE

$(subjDir)/%: FORCE
	proc_pic_act2 $*
	
FORCE:
	
