# Complier
CC=g++
# FLAG
CFLAGS=-std=c++11 -g -Wall
DEBUG1=-DDEBUG1
DEBUG2=-DDEBUG2
DEBUG3=-DDEBUG3
DEBUG4=-DDEBUG4
ARGV1 = data/Viewport.cfg
ARGV2 = data/run.cfg
# Sources files
SRCS = main.cpp AdaptLogic.cpp common.cpp Metadata.cpp Projection.cpp
# Object files
OBJS = $(SRCS:.c=.o)
# Build target
MAIN = main_fix_vp
all: $(MAIN)
	@echo  File has been sucessfully compiled!
# Build
$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(DEBUG3) -o $(MAIN) $(OBJS) 
# 
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@
#
run:
	./$(MAIN) $(ARGV1) $(ARGV2)
# 
clean:
	$(RM) *.o *~ $(MAIN)