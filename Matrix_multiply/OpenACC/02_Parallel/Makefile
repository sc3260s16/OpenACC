CC=pgcc
GPUFLAGS=-acc -Minfo -ta=nvidia:maxwell
MCFLAGS=-acc -Minfo -ta=multicore -fast

BIN = mat_multiply_gpu mat_multiply_multi

all: $(BIN)

mat_multiply_gpu : mat_multiply.c
	$(CC) $(GPUFLAGS) -o $@ $<

mat_multiply_multi : mat_multiply.c
	$(CC) $(MCFLAGS) -o $@ $<

clean:
	$(RM) $(BIN)