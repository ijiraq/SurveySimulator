MODELS = InnerHotModel ReadModelFromFile
# MODELS = ReadModelFromFile
LANGUAGE = F95

%: %.sh

test: clean $(MODELS) HST

$(MODELS): clean
	$(MAKE) -s -C ../$(LANGUAGE) Driver GIMEOBJ=$@ 
	@./test.sh ../$(LANGUAGE)/$@ $(LANGUAGE) 
	$(MAKE) -s -C ../$(LANGUAGE) clean

HST: clean
	$(MAKE) -s -C ../$(LANGUAGE) Driver GIMEOBJ=ReadModelFromFile
	@./test.sh ../$(LANGUAGE)/$@ $(LANGUAGE) 
	$(MAKE) -s -C ../$(LANGUAGE) clean

.PHONY: $(MODELS)

.DEFAULT:

.PHONY: clean

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG test
