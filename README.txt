========================================================================================== 
Test program for Double Buffering for stensil applications.

  2018 Dec. 18, 2019 Feb. 7
  Issaku Kanamori (Hiroshima U., kanamori@hiroshima-u.ac.jp)

  2019 Sep. 10 for public version
  Ken-Ichi Ishikawa (Hiroshima U.,ishikawa@thoe.phys.sci.hiroshima-u.ac.jp)

  2021 Dec.   version with uTofu interface
  Issaku Kanamori (RIKEN R-CCS, kanamori-i@riken.jp)


==========================================================================================

This program solves the linear equation Mx = b of Poisson equation with Jacobi method.

The coefficient matrix M for 2dim Poisson eq. (= Laplacian + const.) is

 (M x)(i,j)  = (4+m2)x(i,j) - x(i+1,j) - x(i-1,j) - x(i,j+1) -x(i,j-1)
               ~~~~~~~~~~~  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  = D x              = H x
  M = D + (L+U) = D + H

Functions hop(out, in) in jacobi2d.cpp and mult_*.cpp gives out += H in.
The size of 2dim. lattice and number of MPI ranks should be defined in "config.h"

The iteration number of Jacobi method is fixed to a constant defiend in "config.h",
and the solver does not monitor the residual to converge.

This program benchmarks the stencil communication overlap in parallelized cases.
This program involves two versions for the stencil communication, (1) single buffering
and (2) double buffering, versions. Buffers are used to exchange the stencil halos.

In order to check the overlap between communcation and computation, and check
the network bandwidth saturatoin, the programs repeats the same solver with
different send/recieve buffer sizes.

* If the buffer size is small enough, all the communication should be hidden
   by the computation, and the elapsed time for communication (and the solver)
   should not depend on the buffer size.

* If the buffer size is large enough, communcation cannot be hidden by
   the computation and the communication time should be propotional
   to the buffer size (+ some latency).  From the communication time, one
   can check if the network bandwidth is saturated.

------------------------------------------------------------------------------------------

Compilation:

  Edit "Makefile" for compiler setting and "config.h" for 2-D lattice size and
  parallelism. Then type make.

Excutables:

  jacobi2d:
     serial version

  mpi_single_buf_jacobi2d:
     single buffering, standard MPI version with persistent communication

  mpi_double_buf_jacobi2d:
     double buffering, standard MPI version with persistent communication

  utofu_double_buf_jacobi2d:
     double buffering, uTofu version

Usage:

  For the serial version, simply run the executable.
  For the parallel MPI versions, run them properly with the MPI launcher.

Outputs:

  To see the consistency, we provide the output logs in ./LOGS directory 
  for each executable. Please compare your results of "norm2" with them.

Files:

  README.txt:
    This file

  LICENCE:
    Licence file
  
  Makefile:
    Edit the first lines to choose the compiler
  
  config.h:
    Defines lattice size, mpi rank size, number of fixed iterations, max buffer size
    (and some constants)
  
  comm_helper.cpp, comm_helper.h:
    Communication related helper functions
  
  field_util.cpp, field_util.h:
    Some utilities: linear algebra related functions
  
  simple_timer.h:
  simple_timer_mpi.h:
  simple_timer_nompi.h:
    Timer with clock_gettime in CLOCK_MONOTONIC mode (nompi)
    or MPI_Wtime (mpi)
    
  solver.h:
    Solver with Jacobi method
  
  mult_*.cpp,  mult_*.h:
    Implementations of the hopping
  
  *jacobi2d.cpp:
    main function

========================================================================================== 
