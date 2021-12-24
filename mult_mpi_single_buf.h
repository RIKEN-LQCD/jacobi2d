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
  Hopping term with single buffer (pesrsistent MPI communication)
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2018 Dec.14
  2018 Dec.29: change the timer to use MPI_Wtime
, 2019 Feb. 6
*/

#ifndef MULT_MPI_SINGLE_BUF_H
#define MULT_MPI_SINGLE_BUF_H

#include "config.h"
#include <mpi.h>
#include "simple_timer.h"

namespace MultMpiSingleBuf {

  class Laplace2d{
public:

  Laplace2d(){ init(); }
  Laplace2d(int buffer_extra){ init(buffer_extra); }
  Laplace2d(int buffer_extra, const int *neighbors){
    init(buffer_extra, neighbors);
  }
  ~Laplace2d(){ finalize(); }
  void init(const int buffer_extra=1, const int *neighbors=0);
  void hop(double *out, const double *in);
  void mult(double *out, const double *in);
  void print_etime();
  void init_timer();
  void set_buffer_extra(const int factor){ m_buffer_extra=factor; }
  static const char *name;
  
private:

  double *rbuf[4], *sbuf[4];

  // send buffer
  double *sbuf_xp;
  double *sbuf_xm;
  double *sbuf_yp;
  double *sbuf_ym;
  
  // receive buffer
  double *rbuf_xp;
  double *rbuf_xm;
  double *rbuf_yp;
  double *rbuf_ym;

  // buffer size
  int size[4];

  // rank for each direction
  int rnode[4];  // for receiveing
  int snode[4];  // for sending

  MPI_Request rrequest[4];
  MPI_Request srequest[4];

  // timers
  simple_timer timer_hop;
  simple_timer timer_communicate;
  simple_timer timer_overlap;
  simple_timer timer_boundary;
  simple_timer timer_pack;
  simple_timer timer_bulk;
  simple_timer timer_irecv;
  simple_timer timer_swait;
  
  void finalize();

  int m_buffer_extra;
};

} //MultMpiSingleBuf

#endif
