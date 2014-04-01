#CXXFLAGS=-ggdb -O0 -Wall
CXXFLAGS=-O2
all: main
main: main.cpp
	(cd third; make tabix)
	(cd base; make)
	g++ $(CXXFLAGS) -I./base -I./third/tabix main.cpp ./base/lib-base.a ./third/tabix/libtabix.a -o main  -lz -lbz2
test:
	./main
clean:
	rm -f main
