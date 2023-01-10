.PHONY: all condainstall

CC  ?= gcc
CXX ?= g++

CXXFLAGS += -std=c++11 -I./ -I./utils
LD_FLAGS += -lz -lpthread -lhts

all : realign_cm merge_result read_extractor
	mkdir -p circlemapcppbin
	cp realign_cm merge_result read_extractor circlemapcppbin/

edlib.o: utils/edlib.cpp utils/edlib.h
	${CXX} -c ${CXXFLAGS} utils/edlib.cpp -o edlib.o

realign_cm: realign_cm.cpp  utils/incr_array.h utils/MultiThread.h  edlib.o utils/edlib.cpp utils/edlib.h
	${CXX} ${CXXFLAGS} realign_cm.cpp -o realign_cm ${LD_FLAGS} edlib.o

merge_result: merge_result.cpp
	${CXX} ${CXXFLAGS} merge_result.cpp  -o merge_result ${LD_FLAGS}

read_extractor: read_extractor.cpp
	${CXX} ${CXXFLAGS} read_extractor.cpp  -o read_extractor ${LD_FLAGS}

clean:
	rm -rf *.o realign_cm merge_result read_extractor
	rm -rf circlemapcppbin

condainstall: all
	echo  "Installed into ${PREFIX}"
	mkdir -p  ${PREFIX}/bin
	cp circle_map++ ${PREFIX}/bin
	cp -r circlemapcppbin ${PREFIX}/bin/
	chmod a+x ${PREFIX}/bin/circle_map++
	chmod a+x ${PREFIX}/bin/circlemapcppbin/*
	mkdir -p  ${PREFIX}/info
	cp LICENSE ${PREFIX}/info/
