
CXX = g++
CFLAGS = -std=c++11 -g -Wall
IDFLAGS = 
LDFLAGS = 

# Detect the operating system
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Linux)
	IDFLAGS += -I/usr/local/include/paraview-5.6
	LDFLAGS += -lvtkCommonCore-pv5.6 -lvtkCommonExecutionModel-pv5.6 -lvtkIOXML-pv5.6 -lvtkCommonDataModel-pv5.6
endif
ifeq ($(UNAME_S), Darwin)
	IDFLAGS += -I/usr/local/opt/vtk/include/vtk-8.2
	LDFLAGS += -lvtkCommonCore-8.2 -lvtkCommonExecutionModel-8.2 -lvtkIOXML-8.2 -lvtkCommonDataModel-8.2
endif

main: ReadFile.o MergeTree.o
	${CXX} ${CFLAGS} ReadFile.o MergeTree.o ${IDFLAGS} ${LDFLAGS} -o bin/main

ReadFile.o: MergeTree.o ReadFile.cpp
	${CXX} ${CFLAGS} ReadFile.cpp MergeTree.o ${IDFLAGS} ${LDFLAGS}

MergeTree.o: MergeTree.cpp MergeTree.h
	${CXX} ${CFLAGS} MergeTree.cpp ${IDFLAGS} ${LDFLAGS}

clean:
	rm -rf *.o bin/* 
