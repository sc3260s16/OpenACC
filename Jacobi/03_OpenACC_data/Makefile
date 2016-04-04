CC=pgcc
CFLAGS=-acc -ta=nvidia:maxwell,time -Minfo=acc

BIN = jacobi_data

all: $(BIN)

jacobi_data : jacobi_data.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	$(RM) $(BIN)
