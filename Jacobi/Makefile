CC       = g++
CCFLAGS  = -O3

BIN =  laplace2d

all: $(BIN)

laplace2d: laplace2d.c
	$(CC) $(CCFLAGS) $(ACCFLAGS) -o $@ $<

clean:
	$(RM) $(BIN)
