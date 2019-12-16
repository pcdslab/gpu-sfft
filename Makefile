 CC=gcc
#For GNU compiler
ifeq ($(CC),gcc)
	CC=gcc
	CXX=g++
	CFLAGS=  -O3 -g -std=gnu99 -Wall  -I/usr/local/cuda-10.1/include -I/usr/local/cuda-10.1/samples/common/inc
	CXXFLAGS= -O3 -g -Wall
	LDFLAGS= -O3 -g
	LDLIBS=  -lm -lstdc++  -L/usr/local/cuda-10.1/lib64 -lcufft -lcudart
endif


#For NVCC Compiler
CUDA_DIR=/usr/local/cuda-10.1
NVCC= $(CUDA_DIR)/bin/nvcc
NVCCFLAGS = -O3 -I$(CUDA_DIR)/include  -I/$(CUDA_DIR)/samples/common/inc -lineinfo  
NVCCLIBS = -L$(CUDA_DIR)/lib64 -lcufft -lcudart -lcudadevrt  

ifeq ($(CC),nvcc)
	CC=$(NVCC)
	CFLAGS = $(NVCCFLAGS) -Xcompiler -std=gnu99  --use_fast_math --compiler-options
	CXX = $(NVCC)
	CXXFLAGS = -O3 --use_fast_math
        LDLIBS = $(NVCCLIBS) -lgomp
endif

C_SRCS       = $(wildcard *.c)
C_OBJS       = $(C_SRCS:.c=.c.o)

CXX_SRCS     = $(wildcard *.cc)
CXX_OBJS     = $(CXX_SRCS:.cc=.cc.o)

CU_SRCS		 = $(wildcard *.cu)
CU_OBJS		 = $(CU_SRCS:.cu=.cu.o)

SRCS         = $(C_SRCS) $(CXX_SRCS) $(CU_SRCS)
OBJS         = $(C_OBJS) $(CXX_OBJS) $(CU_OBJS)

.PHONY: all	clean


.SUFFIXES: .c .cc .cu .o

%.c.o: %.c
	 $(CC) $(CFLAGS) -c $< -o $@

%.cu.o: %.cu
	 $(NVCC) $(NVCCFLAGS) -c $< -o $@

%.cc.o: %.cc
	 $(CXX) $(CXXFLAGS) -c $< -o $@

exec: $(OBJS)
	 $(CC) -o $@ $(OBJS) $(NVCCLIBS) $(LDFLAGS) $(LDLIBS)

all:	exec


clean:
	$(RM) $(OBJS) *.gpu ./exec *.co 



