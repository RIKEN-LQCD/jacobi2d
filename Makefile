##############################################
# choose compiler
##############################################

#compiler = intel
#compiler = gnu
#compiler = fujitsu
compiler = fugaku


##############################################
# rankmap and TNI options
#  -D_USE_TNI_PATTERNx : uses a predefined TNI assignment (x=1,2,3,4)
##############################################
use_rankmap_utofu = yes
TNI_FLAGS=-D_USE_TNI_PATTERN3


#flag=-DDEBUG -g
#flag=-D_USE_RANKMAP
#flag=-D_USE_RANKMAP -D_USE_TNI_PATTERN4
#flag=-D_USE_TNI_PATTERN2




##############################################################
ifeq ($(use_rankmap_utofu),yes)
  flags+=-D_USE_RANKMAP
endif
flags += $(TNI_FLAGS)


TAGET=EMPTY
ifeq ($(compiler),intel)
  TARGET=INTEL
  CXX       = mpicxx
  CXX_NOMPI = icc
  CXXFLAGS  = -O3 -D_MPI_ $(flag)
endif

ifeq ($(compiler),gnu)
  TARGET=GNU
  CXX       = mpic++
  CXX_NOMPI = g++
  CXXFLAGS  = -O3 -D_MPI_ $(flag)
endif

ifeq ($(compiler),fujitsu)
  TARGET=FUJITSU
  CXX       = mpiFCCpx
  CXX_NOMPI = FCCpx
  CXXFLAGS  = -D_MPI_ -Xg -Kfast,-restp=all,optmsg=2,ocl,preex,unroll=100 -Nline,lst=t  $(flags)
endif

ifeq ($(compiler),fugaku)
  TARGET=FUGAKU
  CC        = mpifccpx
  CXX       = mpiFCCpx
  CXX_NOMPI = FCCpx
  CFLAGS  = -D_MPI_  -D_UTOFU_RDMA -Kfast,openmp -I$(RDMA_LIB) -Nline,lst=t $(flags)
  CXXFLAGS = $(CFLAGS)
  TOFUFLAGS = -ltofucom
endif

ifeq ($(TARGET),EMPTY)
  @echo unknown comiler
  exit
endif


PROGRAM_SERIAL    = jacobi2d
PROGRAMS_PARALLEL = mpi_single_buf_jacobi2d mpi_double_buf_jacobi2d
PROGRAM_UTOFU     = utofu_double_buf_jacobi2d
#utofu_single_buf_jacobi2d

ifeq ($(use_rankmap_utofu),yes)
  RANKMAP_LIB = rankmap_utofu_2d
  RANKMAP_A = $(RANKMAP_LIB)/librankmap_utofu.a
  RANKMAP_LDFLAGS = -L$(RANKMAP_LIB) -lrankmap_utofu
else
  OBJS_RANKMAP = rankmap_lib.o
  RANKMAP_LIB =
endif

SRC_COMMON   = field_util.cpp comm_helper.cpp

SRC_MPI_DB  := mult_mpi_double_buf.cpp mpi_double_buf_jacobi2d.cpp $(SRC_COMMON)
OBJS_MPI_DB := $(SRC_MPI_DB:%.cpp=%.o) $(OBJS_RANKMAP)
DEPS_MPI_DB := $(OBJS_MPI_DB:%.o=%.d)

SRC_MPI_SB  := mult_mpi_single_buf.cpp mpi_single_buf_jacobi2d.cpp $(SRC_COMMON)
OBJS_MPI_SB := $(SRC_MPI_SB:%.cpp=%.o)  $(OBJS_RANKMAP)
DEPS_MPI_SB := $(OBJS_MPI_SB:%.o=%.d)


RDMA_LIB = rdma_utofu
SRC_UTOFU_DB  := utofu_double_buf_jacobi2d.cpp mult_utofu_double_buf.cpp rdma_comlib_2buf.cpp $(SRC_COMMON)
OBJS_UTOFU_DB := $(SRC_UTOFU_DB:%.cpp=%.o) $(OBJS_RANKMAP)
DEPS_UTOFU_DB := $(OBJS_UTOFU_DB:%.o=%.d)



all : $(PROGRAMS_PARALLEL) $(PROGRAM_SERIAL) $(PROGRAM_UTOFU)
#all : $(PROGRAMS_PARALLEL) $(PROGRAM_SERIAL)
#all : $(PROGRAM_UTOFU)

#export CC
#export CFLAGS

$(RDMA_LIB)/librdma_comlib.a:
	$(MAKE) -C $(RDMA_LIB)

$(RANKMAP_LIB)/librankmap_utofu.a:
	$(MAKE) -C $(RANKMAP_LIB) RANKMAP_FLAGS="$(RANKMAP_FLAGS)"

jacobi2d: jacobi2d.cpp
	$(CXX_NOMPI) $(CXXFLAGS) $< -o jacobi2d

mpi_single_buf_jacobi2d: $(OBJS_MPI_SB) $(RANKMAP_A)
	$(CXX) $(CXXFLAGS) $(TOFUFLAGS) -o $@ $(OBJS_MPI_SB) $(RANKMAP_LDFLAGS) 

mpi_double_buf_jacobi2d: $(OBJS_MPI_DB) $(RANKMAP_A)
	$(CXX) $(CXXFLAGS) $(TOFUFLAGS) $(RANKMAP_LDFLAGS) -o $@ $(OBJS_MPI_DB) $(RANKMAP_LDFLAGS) 

utofu_double_buf_jacobi2d: $(OBJS_UTOFU_DB) $(RANKMAP_A) $(RDMA_LIB)/librdma_comlib.a
	$(CXX) $(CXXFLAGS) $(TOFUFLAGS) -L$(RDMA_LIB) -o $@ $(OBJS_UTOFU_DB) -lrdma_comlib $(RANKMAP_LDFLAGS)

clean:
	rm -f $(OBJS_MPI_SB) $(DEPS_MPI_SB) \
	      $(OBJS_MPI_DB) $(DEPS_MPI_DB) \
	      $(OBJS_UTOFU_DB) $(DEPS_UTOFU_DB) \
	      $(OBJS_UTOFU_SB) $(DEPS_UTOFU_SB) \
	      $(PROGRAMS_PARALLEL) $(PROGRAM_UTOFU) $(PROGRAM_SERIAL) *~ *.bak *.lst

dist-clean: clean
	$(MAKE) -C $(RDMA_LIB) clean
	$(MAKE) -C $(RANKMAP_LIB) clean

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -MMD -MP $<

%.o: %.c
	$(CC) $(CFLAGS) -c -MMD -MP $<

-include $(DEPS_MPI)
-include $(DEPS_MPI_SB)
-include $(DEPS_MPI_DB)
-include $(DEPS_UTPFI_DB)
