CC=pgcc
CFLAGS=-acc -Minfo=acc -ta=nvidia:maxwell -fast

mat_multiply : mat_multiply.c
	$(CC) -o $@ $< $(CFLAGS)
	rm -f mat_multiply_ipa12_mat_multiply.oo mat_multiply.o

.PHONY: clean

clean:
	rm mat_multiply