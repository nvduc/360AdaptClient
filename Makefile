CXX=g++
FLAGS=-std=c++11 -g
DEBUG1=-DDEBUG1
DEBUG2=-DDEBUG2
DEBUG3=-DDEBUG3
DEBUG4=-DDEBUG4

main_fix_vp:
	${CXX} ${FLAGS} ${DEBUG3} ${FLAGS} -o main_fix_vp main.cpp AdaptLogic.cpp common.cpp Metadata.cpp Projection.cpp

