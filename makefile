CC = g++
INC_DIR = include
CFLAGS = -c -Wall -std=c++11 -fno-exceptions -O2 -I$(INC_DIR)
DEPS = $(INC_DIR)/curve.h $(INC_DIR)/Number.h $(INC_DIR)/xmatrix.h $(INC_DIR)/betti.h
objects = main.o src/curve.o src/Number.o src/xmatrix.o src/betti.o


syzygies:	$(objects)
	$(CC) -o syzygies $(objects) -s -lpthread

%.o: 		%.cpp $(DEPS)
	$(CC) -o $@ $< $(CFLAGS)

clean:		
	rm -f $(objects) syzygies

