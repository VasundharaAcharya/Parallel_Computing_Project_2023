ALL:exp.c par.cu
	mpixlc -O3 exp.c -c -o exp.o -lm -D DEBUG_MODE
	nvcc -O3 par.cu -c -o par.o -D DEBUG_CUDA
	mpixlc -O3 exp.o par.o -o exp-exe \
-L/usr/local/cuda-11.2/lib64/ -lcudadevrt -lcudart -lstdc++
