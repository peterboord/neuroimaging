#!/bin/bash
subjDir=$(shell pwd)
subjPaths=$(wildcard $(subjDir)/RC4[1-2][0-9][0-9]-2)
.PHONY: all

all: $(subjPaths)
	
$(subjDir)/%: FORCE
	proc_pic $*

FORCE:
	
