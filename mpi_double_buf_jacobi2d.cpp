//****************************************************************************************
//
//  Copyright (c) 2019-2021, Ken-Ichi Ishikawa (ishikawa@theo.phys.sci.hirosima-u.ac.jp)
//  Copyright (c) 2019, Issaku Kanamori   (kanamori@hiroshima-u.ac.jp)
//  Copyright (c) 2019-2021, Issaku Kanamori   (kanamori-i@riken.jp)
//  Copyright (c) 2021, Yoshifumi Nakamura <nakamura@riken.jp>
//  Copyright (c) 2021, Yuta Mukai         <mukai.yuta@fujitsu.com>
//
//
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are
//  met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer. 
//
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer listed
//    in this license in the documentation and/or other materials
//    provided with the distribution.
//
//  * Neither the name of the copyright holders nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
//  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
//  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
//
//----------------------------------------------------------------------------------------
//  ACKNOWLEDGMENT
//  
//  This software has been developed in a co-design working group for the lattice QCD
//  application on Post K computer and supported by Priority Issue 9 to be tackled
//  by Using Post K Computer.
//
//****************************************************************************************
/*
  solving 2dim Laplace eq. with Jacobi method: double buffers with MPI functions
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2018 Dec.14, 2019 Feb. 6
*/


#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "comm_helper.h"
#include "simple_timer.h"
#include "field_util.h"
#include "mult_mpi_double_buf.h"
#include "solver.h"
#ifdef _USE_RANKMAP
extern "C" {
#include "rankmap_lib.h"
}
#endif

int main(int argc, char** argv){

  if(argc<3){
    fprintf(stderr, "usage: %s Px Py\n", argv[0]);
    return -1;
  }
  int Px=0;
  int Py=0;

  MPI_Init(&argc, &argv);

  simple_timer timer;
  timer.start();


  // set some global variables: rankid etc.
  communicator::init(false);
  if(argc>2){
    Px=atoi(argv[1]);
    Py=atoi(argv[2]);
    if(communicator::rankid==0){
      printf("from command line: Px=%d, Py=%d\n", Px, Py);
    }
  }


  // fields
  double source[Nx*Ny];
  double solution[Nx*Ny];

  int proc_size[2]={Px,Py};
  int proc_coord[2]={-1,-1};
  int proc_neighbors[4];


  {
#ifdef _USE_RANKMAP
    int map_id=rankmap_lib_set_rankmap2d_size(proc_size);
#else
    int map_id=-1;
#endif
    int err=0;
    if(map_id<0){
      err=map_id;
      int myrank=communicator::rankid;
      if(myrank==0){
        printf("bad map_id: %d\n", map_id);
        printf("using the default lexical rank map\n");
        fflush(0);
      }
      int npx=myrank % Px;
      int npy=myrank / Px;
      int np_xp=(npx+1) % Px;
      int np_xm=(npx-1+Px) % Px;
      int np_yp=(npy+1) % Py;
      int np_ym=(npy-1+Py) % Py;

      int pxp = np_xp + Px*npy;
      int pxm = np_xm + Px*npy;
      int pyp = npx + Px*np_yp;
      int pym = npx + Px*np_ym;
      proc_coord[0]=npx;
      proc_coord[1]=npy;
      proc_neighbors[0]=pxp;
      proc_neighbors[1]=pxm;
      proc_neighbors[2]=pyp;
      proc_neighbors[3]=pym;
    } else {
#ifdef _USE_RANKMAP
      rankmap_lib_get_rankmap(proc_coord, proc_neighbors, proc_size);
      Px=proc_size[0];
      Py=proc_size[1];
#endif
    }
    NNx=Nx*Px;
    NNy=Ny*Py;
    if(communicator::rankid==0){
      printf("Px, Py = %d, %d; NNx, NNy = %d, %d\n", Px, Py, NNx, NNy);
      fflush(0);
    }
  }
  if(communicator::rankid==0){
    printf("before set_field: Px, Py = %d, %d; NNx, NNy = %d, %d\n", Px, Py, NNx, NNy);
    fflush(0);
  }
  int np;
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  if(np != Px*Py){
    char msg[256];
    sprintf(msg,"wrong number of nodes, np=%d: must be %d\n", np, Px*Py);
    communicator::error_stop(msg);
  }


  set_field(source, proc_coord[0], proc_coord[1], Px, Py);
  MPI_Barrier(MPI_COMM_WORLD);
  int factor=1;
  for(int i=0; i<FactIterMax; i++)
  {
    using namespace MultMpiDoubleBuf;
    if(communicator::rankid==0){
      printf("factor for extra buffer: %d\n", factor);
    }
    communicator::barrier();
    Laplace2d  mat(factor, proc_neighbors);
    communicator::barrier();
    solve(solution, source, mat);
    factor*=2;
  }

  // done
  MPI_Barrier(MPI_COMM_WORLD);
  timer.stop();
  MPI_Finalize();

  double etime=timer.get_elapsed_msec();
  if(communicator::rankid==0){
    printf("total elapsed time:  %f [msec]\n", etime);
  }

  return 0;

}

