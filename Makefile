
CXX = g++
CFLAGS = -std=c++11 -g -Wall -fopenmp

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

all: serial parallel

serial: SerialMain.cpp MergeTree.cpp Utils.cpp
	${CXX} ${CFLAGS} $^ ${IDFLAGS} ${LDFLAGS} -o bin/$@

parallel: ParallelMain.cpp MergeTree.cpp Utils.cpp
	${CXX} ${CFLAGS} $^ ${IDFLAGS} ${LDFLAGS} -o bin/$@

clean:
	rm -rf *.o bin/* 
