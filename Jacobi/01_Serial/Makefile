CC=pgcc
CFLAGS=-Mipa=fast -fast

BIN = jacobi_serial

all: $(BIN)

jacobi_serial : jacobi_serial.c
	$(CC) $(CFLAGS) -o $@ $<
	$(RM) jacobi_serial_ipa13_jacobi_serial.oo  jacobi_serial.o 

clean:
	$(RM) $(BIN)
