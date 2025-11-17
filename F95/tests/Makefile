MODELS = InnerHotModel ReadModelFromFile
# MODELS = ReadModelFromFile
LANGUAGE = F95

%: %.sh

test: clean $(MODELS) HST

$(MODELS): clean
	$(MAKE) -s -C ../ Driver GIMEOBJ=$@ 
	@./test.sh ../$@ $(LANGUAGE) 
	$(MAKE) -s -C ../ clean

HST: clean
	$(MAKE) -s -C ../ Driver GIMEOBJ=ReadModelFromFile
	@./test.sh ../$@ $(LANGUAGE) 
	$(MAKE) -s -C ../ clean

.PHONY: $(MODELS)

.DEFAULT:

.PHONY: clean

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG test
