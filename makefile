CC=g++
INC_DIR = include
CFLAGS = -c -Wall -fexceptions -I$(INC_DIR)
DEPS = $(INC_DIR)/curve.h $(INC_DIR)/Number.h $(INC_DIR)/xmatrix.h


syzygies:	main.o src/curve.o src/Number.o src/xmatrix.o
	$(CC) -o syzygies main.o src/curve.o src/Number.o src/xmatrix.o -s

%.o: 		%.cpp $(DEPS)
	$(CC) -o $@ $< $(CFLAGS)

clean:		
	rm -rf *.o src/*.o syzygies

