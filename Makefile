
.PHONY: all

all : realign_cm merge_result

edlib.o: utils/edlib.cpp utils/edlib.h
	g++ -std=c++11 -I./ -I./utils -g -c -Wall utils/edlib.cpp -o edlib.o

realign_cm: realign_cm.cpp  utils/incr_array.h utils/MultiThread.h  edlib.o utils/edlib.cpp utils/edlib.h
	g++ -std=c++11 -I./ -g -c -Wall realign_cm.cpp -o realign_cm.o
	g++ -std=c++11 -lhts -lz -lpthread edlib.o realign_cm.o -o realign_cm

merge_result: merge_result.cpp
	g++ -std=c++11 -g -c -Wall merge_result.cpp -o merge_result.o
	g++ -std=c++11 merge_result.o  -lpthread  -o merge_result 
        

clean:
	rm -rf *.o realign_cm merge_result 
