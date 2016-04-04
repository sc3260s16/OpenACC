CC=nvcc
CFLAGS=-arch sm_52

mat_multiply : mat_multiply.cu
	$(CC)  $(CFLAGS) -o $@ $<

.PHONY: clean

clean:
	rm mat_multiply