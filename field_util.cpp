//****************************************************************************************
//
//  Copyright (c) 2019-21, Ken-Ichi Ishikawa (ishikawa@theo.phys.sci.hirosima-u.ac.jp)
//  Copyright (c) 2019, Issaku Kanamori   (kanamori@hiroshima-u.ac.jp)
//  Copyright (c) 2019-21, Issaku Kanamori   (kanamori-i@riken.jp)
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
  some utilities: linear algebra related functions
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2018 Dec.14, 2019 Feb. 6
*/

#include "config.h"
#include "comm_helper.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef _MPI_
#include <mpi.h>
#endif

namespace{
  // thread global variables
  double sum2;
  double result;
}

// squared norm
double norm2(const double *in){
#pragma omp barrier
#pragma omp master
  {
    sum2=0.0;
    result=0.0;
  }
#pragma omp barrier
#pragma omp for reduction(+:sum2)
  for(int i=0; i<Nx*Ny; i++){
    sum2+=in[i]*in[i];
  }
#pragma omp master
  {
#ifdef _MPI_
  communicator::timer[communicator::Timer::ALL_REDUCE].start();
  MPI_Allreduce(&sum2, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  communicator::timer[communicator::Timer::ALL_REDUCE].stop();
#else
  result=sum2;
#endif
  }
#pragma omp barrier
  return result;
}

// initialize with random numbers
//   for a given (NNx, NNy), this routine gives a unique value
//   independet from (Px, Py) and gives the the same on in the serial version.
void set_field(double *v, int npx, int npy, int Px, int Py){
  int idx=0;
  int rankid=communicator::rankid;
  if(npx<0){
    npx=rankid % Px;
  }
  if(npy<0){
    npy=rankid / Px;
  }
  double n2=0.0;
  //  printf("[npx=%d, npy=%d]\n",npx, npy);
  //printf("rankid=%d, Px=%d, Py=%d\n", rankid, Px, Py);
  NNx=Nx*Px;
  NNy=Ny*Py;
  srand(1);
  for(int y=0; y<NNy; y++){
    for(int x=0; x<NNx; x++){
      double val=( (double) rand())/RAND_MAX;
      n2+=val*val;
      if(x / Nx != npx) { continue; }
      if(y / Ny != npy) { continue; }
      v[idx]=val;
      idx++;
    }
  }
  if(idx!=Nx*Ny){
    printf("cannot happen [npx=%d, npy=%d], something is wrong: idx=%d  (must be %d)\n",npx, npy, idx, Nx*Ny);
  }
  if(communicator::rankid==0){
    printf("set_field: n2 before normalizing=%e\n", n2);
  }

  double factor=sqrt(1.0/n2);
  for(int i=0; i<Nx*Ny; i++){
    v[i]*=factor;
  }
}
