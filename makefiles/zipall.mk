outFile=bigNiiPaths
niiList=$(shell cat $(outFile))
niiGzList=$(niiList:%=%.gz)
.PHONY: all
all: $(niiGzList)

%.gz:
	cmd="gzip $*" ;\
	echo $$cmd ;\
	eval $$cmd
