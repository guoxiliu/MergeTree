CXX = g++
CFLAGS = -std=c++11 -g -Wall
IDFLAGS = -I/usr/local/opt/vtk/include/vtk-8.2
LDFLAGS = -lvtkCommonCore-8.2 -lvtkCommonExecutionModel-8.2 -lvtkIOXML-8.2 

PROJECT = ReadUnstructuredGrid

${PROJECT}: ${PROJECT}.cpp
	${CXX} ${CFLAGS} -o bin/${PROJECT} ${PROJECT}.cpp ${IDFLAGS} ${LDFLAGS}

clean:
	rm -rf ${PROJECT}