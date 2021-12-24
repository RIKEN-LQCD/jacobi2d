CC = mpifccpx
CFLAGS = -Kfast -D_UTOFU_RDMA -D_USE_RANKMAP -D_UTOFU_RANKMAP $(TNI_FLAGS)
AR      = ar
RANLIB  = ranlib

PROG = librankmap_utofu.a
SRC = get_tni_2d.c get_tofu_coord_2d.c  rankmap_lib_utofu.c
OBJ = $(SRC:%.c=%.o)
HDR = 

all : $(PROG)

.o : .c $(HDR)
	$(CC) $(CFLAGS) -c $<

$(PROG) : $(OBJ)
	rm -f  $@
	$(AR) -r $@  $(OBJ)
	$(RANLIB) $@

clean :
	rm -f *.o $(PROG)




