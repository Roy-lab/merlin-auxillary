LFLAG = -lgsl -lgslcblas 
SRC = estimateEdgeConf.C

CC=g++
CFLAGS = -g

estimateEdgeConf: $(SRC)
	$(CC) $(SRC) $(LFLAG) $(CFLAGS) -o estimateEdgeConf

clean:
	rm -f estimateEdgeConf *~
