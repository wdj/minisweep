
CFLAGS = -O3 -acc -Minfo=accel -ta=tesla,cc60 

#all: atomics locks locks_llvm locks_xlc atomics_xlc
TARGETS = gauss_seidel gauss_seidel_acc

OBJECTS = $(addsuffix .o, $(TARGETS))

all: $(TARGETS)

ifdef SIZE
    MYSIZE=$(SIZE)
else
    MYSIZE=10
endif
CFLAGS += -DSIZE=$(MYSIZE) 

ifdef SIZE_U
    MYSIZE_U=$(SIZE_U)
else
    MYSIZE_U=10
endif
CFLAGS += -DSIZE_U=$(MYSIZE_U)

ifdef BLOCK_GRID_W
    BLOCK_GRID_W=$(BLOCK_GRID_W)
else
    BLOCK_GRID_W=1
endif
CFLAGS += -DBLOCK_GRID_W=$(BLOCK_GRID_W)

ifdef BLOCK_GRID_H
    BLOCK_GRID_H=$(BLOCK_GRID_H)
else
    BLOCK_GRID_H=1
endif
CFLAGS += -DBLOCK_GRID_H=$(BLOCK_GRID_H)

CFLAGS += -DCOUNTEXTENT=$(COUNTEXTENT)

%:%.c
	pgcc $(CFLAGS) $< -o $@ -lm

clean:
	rm -f $(TARGETS) $(OBJECTS)

