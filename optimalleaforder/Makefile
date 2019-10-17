LFLAG = -lgsl -lgslcblas 
SRC = Distance.C  Framework.C  Heap.C  HierarchicalCluster.C  HierarchicalClusterNode.C  OptimalLeafOrder.C

CC=g++
CFLAGS = -g

reorder: $(SRC)
	$(CC) $(SRC) $(LFLAG) $(CFLAGS) -o reorder

clean:
	rm -f reorder *~
