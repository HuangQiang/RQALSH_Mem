SRCS=util.cc random.cc pri_queue.cc rqalsh.cc rqalsh_star.cc \
	afn.cc main.cc
OBJS=$(SRCS:.cc=.o)

CXX=g++ -std=c++11
CPPFLAGS=-w -O3

.PHONY: clean

all: ${OBJS}
	${CXX} ${CPPFLAGS} -o rqalsh ${OBJS}

util.o: util.h

random.o: random.h

pri_queue.o: pri_queue.h

rqalsh.o: rqalsh.h

rqalsh_star.o: rqalsh_star.h

afn.o: afn.h

main.o:

clean:
	-rm $(OBJS) rqalsh
