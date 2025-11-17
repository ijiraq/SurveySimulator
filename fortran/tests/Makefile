MODELS = InnerHotModel ReadModelFromFile
# MODELS = ReadModelFromFile
LANGUAGE = F95

%: %.sh

test: clean $(MODELS) HST

$(MODELS): clean
	$(MAKE) -s -C ../src/$(LANGUAGE) Driver GIMEOBJ=$@ 
	@./test.sh ../src/$(LANGUAGE)/$@ $(LANGUAGE) 
	$(MAKE) -s -C ../src/$(LANGUAGE) clean

HST: clean
	$(MAKE) -s -C ../src/$(LANGUAGE) Driver GIMEOBJ=ReadModelFromFile
	@./test.sh ../src/$(LANGUAGE)/$@ $(LANGUAGE) 
	$(MAKE) -s -C ../src/$(LANGUAGE) clean

.PHONY: $(MODELS)

.DEFAULT:

.PHONY: clean

clean:
	\rm -f ModelUsed.dat SimulDetect.dat SimulTrack.dat LOG test
