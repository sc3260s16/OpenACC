CC=pgcc
CFLAGS=-acc -ta=nvidia:maxwell,time -Minfo=acc

BIN = jacobi_kernels

all: $(BIN)

jacobi_kernels : jacobi_kernels.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	$(RM) $(BIN)
