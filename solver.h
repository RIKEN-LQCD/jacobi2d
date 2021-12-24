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
  solver with Jacobi method
    Issaku Kanamori kanamori@hiroshima-u.ac.jp
    2018 Dec.14, 2019 Feb. 6

  requires:
    template class MAT has
      void hop(*double, const *double)
        hop(inout, in) gives
        inout(x):= inout(x) + sum_mu( in(x+mu) + in(x-mu) )
      void print_etime()
      static const char *name

    global constants Nx, Ny, mass2, NiterForSolver are given (config.h)

 */

#ifndef SOLVER_H
#define SOLVER_H
#include "config.h"
#include <omp.h>
#include "comm_helper.h"
#include "simple_timer.h"
#include "field_util.h"


template<typename MAT>
void solve(double* solution, const double* source,
           MAT &mat){
  int n_iter=IterForSolver;
  double v[Nx*Ny]; // work vector
  double factor=-1.0/(4.0+mass2);
  simple_timer timer_solver;
  simple_timer timer_norm2;
  simple_timer &timer_allreduce=communicator::timer[communicator::Timer::ALL_REDUCE];
  timer_solver.reset();
  int norm2_count=0;

  if(communicator::rankid==0){
    printf("---------------------------------------------\n");
    printf("   solving with %s\n", mat.name);
    printf("---------------------------------------------\n");
    printf("Nx=%d, Ny=%d, n_iter=%d\n", Nx, Ny, n_iter);
    fflush(stdout);
  }

#pragma omp parallel
  {
    double n2_src=norm2(source);
#pragma omp master
    {
      communicator::barrier();
      timer_norm2.reset();
      timer_allreduce.reset();
      timer_solver.start();
    }

  // initial guess
#pragma omp for
    for(int xy=0; xy<Nx*Ny; xy++){
      solution[xy]=0.0;
    }

    for(int i=0; i<n_iter; i++){

#pragma omp for
      for(int xy=0; xy<Nx*Ny; xy++){
        v[xy]=-source[xy];
      }

      mat.hop(v,solution);

#pragma omp barrier
#pragma omp for
      for(int xy=0; xy<Nx*Ny; xy++){
        solution[xy]=factor*v[xy];
      }

      if(i%norm2_freq == 0){
#pragma omp master
        {
          timer_norm2.start();
        }
        double n2=norm2(solution);
#pragma omp master
        {
          timer_norm2.stop();
          norm2_count++;
          if(communicator::rankid==0){
            printf("iteration=%d, norm2(solution): %g\n",i, n2);
          }
        }
      }
#pragma omp barrier
    } // i

#pragma omp master
    {
      communicator::barrier();
      timer_solver.stop();
      double etime_solver=timer_solver.get_elapsed_msec();
      double etime_norm2=timer_norm2.get_elapsed_msec();
      double etime_allreduce=timer_allreduce.get_elapsed_msec();
      if(communicator::rankid==0){
        printf("etime_solver:  %f [msec]  (%d iterations)\n", etime_solver,n_iter);
        printf("etime_norm2:       %f [msec]   (%d times)\n",etime_norm2, norm2_count);
        printf("etime_allreduce:   %f [msec]   (%d times)\n",etime_allreduce, norm2_count);
        mat.print_etime();
      }
    }
#pragma omp barrier
    mat.mult(v,solution);
#pragma omp barrier
#pragma omp for
    for(int xy=0; xy<Nx*Ny; xy++){
      v[xy]-=source[xy];
    }
    double n2_r=norm2(v);
    double n2  =norm2(solution);
#pragma omp master
    {
      if(communicator::rankid==0){
        printf("done: %s\n", mat.name);
        printf("iteration=%d, relative diff2: %g\n",n_iter, n2_r/n2_src);
        printf("norm2 of the solution=%24.15e\n", n2);
        printf("\n");
      }
    }
  } // omp parallel

}

#endif
