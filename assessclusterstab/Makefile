LFLAG = -lgsl -lgslcblas 
SRC = Framework.C

CC=g++
CFLAGS = -g

assessClusterStab: $(SRC)
	$(CC) $(SRC) $(LFLAG) $(CFLAGS) -o assessClusterStab

clean:
	rm -f assessClusterStab *~
