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
  communication realted helpler funcionts
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
    2018 Dec.14, 2019 Feb. 6
*/


#include "config.h"
#include "comm_helper.h"
#include <stdlib.h>
#include <stdio.h>
#ifdef _MPI_
#include <mpi.h>
#endif

int Px=0;
int Py=0;
int NNx;
int NNy;

namespace communicator{
  // this node
  int rankid;

  simple_timer timer[Timer::NUM]; // global timer

  // safe error exit
  void error_stop(const char* err_msg){
    if(communicator::rankid==0){
      printf("%s", err_msg);
      printf("aborting...\n");
      fflush(0);
    }
#ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    exit(0);
  }

  // print message to stdout
  void debug_print(const char* msg){
#ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
    if(rankid==0){
      printf("%s\n", msg);
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#else
    printf("%s\n", msg);
    fflush(stdout);
#endif
  }

  // set rankid and sanity check
  void init(const bool err_stop){
    int np, rankid;
#ifdef _MPI_
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankid);
#else
    np=1;
    rankid=1;
#endif
    communicator::rankid=rankid;
    if(np != _Pnum){
      if(err_stop){
        char msg[256];
        sprintf(msg,"wrong number of nodes, np=%d: must be %d\n", np, _Pnum);
        error_stop(msg);
      }
    }
    if(communicator::rankid==0){
      printf("NNx, NNy= %d, %d\n",NNx,NNy);
      printf("Nx,  Ny = %d, %d\n",Nx,Ny);
      //      printf("Px,  Py = %d, %d\n",Px,Py);
    }
    for(int i=0; i<communicator::Timer::NUM; i++){
      timer[i].reset();
    }
  }

  // barrier
  void barrier(){
#ifdef _MPI_
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  double reduction(const double x){
    double sum=x;
#ifdef _MPI_
    MPI_Allreduce(&x, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    return sum;
  }

} //communicator

