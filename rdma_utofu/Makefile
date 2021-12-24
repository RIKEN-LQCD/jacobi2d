CC = mpifccpx
#CFLAGS = -Kfast
AR      = ar
RANLIB  = ranlib

PROG = librdma_comlib.a 
SRC = rdma_utofu_comlib.c
OBJ = $(SRC:%.c=%.o)
HDR = rdma_utofu_comlib.h utofu_tool.h

all : $(PROG)

$(OBJ) : $(SRC) $(HDR)
	$(CC) $(CFLAGS) -c $<

$(PROG) : $(OBJ)
	rm -f  $@
	$(AR) -r $@  $(OBJ)
	$(RANLIB) $@

clean :
	rm -f *.o $(PROG)



