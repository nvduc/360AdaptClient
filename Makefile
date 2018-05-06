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
MAIN_FIX_VP = main_fix_vp
MAIN_REAL_TRACE = main_real_trace
all: $(MAIN_FIX_VP) $(MAIN_REAL_TRACE)
	@echo  File has been sucessfully compiled!
# Build
$(MAIN_FIX_VP): $(OBJS)
	$(CC) $(CFLAGS) $(DEBUG3) -o $(MAIN_FIX_VP) $(OBJS) 
# 
$(MAIN_REAL_TRACE): $(OBJS)
	$(CC) $(CFLAGS) $(DEBUG2) -o $(MAIN_REAL_TRACE) $(OBJS)
# 
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@
#
run_fix_vp:
	./$(MAIN_FIX_VP) $(ARGV1) $(ARGV2)
# 
run_real_trace:
	./$(MAIN_REAL_TRACE) $(ARGV1) $(ARGV2)
# 
clean:
	$(RM) *.o *~ $(MAIN_FIX_VP) ${MAIN_REAL_TRACE}