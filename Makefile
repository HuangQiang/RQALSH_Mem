SRCS=util.cc random.cc pri_queue.cc qdafn.cc drusilla_select.cc \
	rqalsh.cc rqalsh_star.cc ml_rqalsh.cc afn.cc main.cc
OBJS=${SRCS:.cc=.o}

CXX=g++ -std=c++11
CPPFLAGS=-w -O3

.PHONY: clean

all: ${OBJS}
	${CXX} ${CPPFLAGS} -o rqalsh ${OBJS}

util.o: util.h

random.o: random.h

pri_queue.o: pri_queue.h

qdafn.o: qdafn.h

drusilla_select.o: drusilla_select.h

rqalsh.o: rqalsh.h

rqalsh_star.o: rqalsh_star.h

ml_rqalsh.o: ml_rqalsh.h

afn.o: afn.h

main.o:

clean:
	-rm ${OBJS}
