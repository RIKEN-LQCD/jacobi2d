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
  Hopping term with double buffer (RDMA version)
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2018 Dec.17, 2019 Feb. 6
*/

#ifndef MULT_UTOFU_DOUBLE_BUF_H
#define MULT_UTFU_DOUBLE_BUF_H

#include "config.h"
#include <stdint.h>
#include "simple_timer.h"
#include "rdma_comlib_2buf.h"

namespace MultUtofuDoubleBuf{
class Laplace2d{
public:

  Laplace2d(){ init(); }
  Laplace2d(int buffer_extra, const int *proc_neighbors){
    init(buffer_extra, proc_neighbors);
  }
  ~Laplace2d(){ finalize(); }
  void init(const int buffer_extra=1, const int *proc_neighbors=0);
  void hop(double *out, const double *in);
  void mult(double *out, const double *in);
  void print_etime();
  void init_timer();
  static const char *name;

private:

  // buffer + communicator
  rdma_comlib_2buf buff_rdma[4];

  // buffer size
  int size[4];

  // rank for each direction
  int rnode[4];  // for receiveing
  int snode[4];  // for sending

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

}//MultUtofuDoubleBuf

#endif
