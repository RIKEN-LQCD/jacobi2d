//****************************************************************************************
//
//  Copyright (c) 2019, Ken-Ichi Ishikawa (ishikawa@theo.phys.sci.hirosima-u.ac.jp)
//  Copyright (c) 2019, Issaku Kanamori   (kanamori@hiroshima-u.ac.jp)
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
  simple timer class with MPI_Wtime
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2018 Dec.29, 2019 Feb. 6

  member functions

  void reset()
    reset the elapsed time

  void start()
    start the timer

  void stop()
    stop the timer and accumlate the elapsed time

  double get_elapsed_msec()
    return the elapsed tiem in unit of msec

  void print_resolution()
    print the resutlion of the timer

*/

#ifndef SIMPLE_TIMER_MPI_INC_H
#define SIMPLE_TIMER_MPI_INC_H
#include <mpi.h>
#include <stdio.h>
class simple_timer_mpi{

 public:
 simple_timer_mpi(): diff_sec(0) {}

  void reset(){
    diff_sec=0.0;
  }

  void start(){
    start_time=MPI_Wtime();
  }
  void stop(){
    end_time=MPI_Wtime();
    diff_sec+=(end_time-start_time);
  }

  double get_elapsed_msec() const {
    double diff_total_msec=diff_sec*1.0e+3;
    return diff_total_msec;
  }

  void print_resolution(){
    double ticks=MPI_Wtick();
    double msec=1.0e3*ticks;
    printf("time resolution: %g msec.  (MPI version)\n",msec);
  }

 private:
  double start_time, end_time;
  double diff_sec;
};
#endif
