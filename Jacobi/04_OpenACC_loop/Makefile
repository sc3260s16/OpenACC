CC=pgcc
CFLAGS=-acc -ta=nvidia:maxwell,time -Minfo=acc

BIN = jacobi_loop

all: $(BIN)

jacobi_data : jacobi_loop.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	$(RM) $(BIN)
