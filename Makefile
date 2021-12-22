
.PHONY: all



all : realign_cm check_1bp_cov


edlib.o: utils/edlib.cpp utils/edlib.h
	g++ -std=c++11 -I./ -I./utils -g -c -Wall utils/edlib.cpp -o edlib.o

realign_cm: realign_cm.cpp  utils/incr_array.h utils/MultiThread.h  edlib.o utils/edlib.cpp utils/edlib.h
	g++ -std=c++11 -I./   -g -c -Wall realign_cm.cpp -o realign_cm.o
	g++ -std=c++11 -lhts -lz -lpthread edlib.o realign_cm.o -o realign_cm

check_1bp_cov.o: check_1bp_cov.cpp
	g++ -std=c++11 -g -c -Wall check_1bp_cov.cpp -o check_1bp_cov.o

check_1bp_cov: check_1bp_cov.o
	g++ -std=c++11 -lhts -lz -lpthread check_1bp_cov.o -o check_1bp_cov


clean:
	rm -rf *.o realign_cm check_1bp_cov
