//****************************************************************************************
//
//  Copyright (c) 2015-2021, Yoshifumi Nakamura <nakamura@riken.jp>
//  Copyright (c) 2015-2021, Yuta Mukai         <mukai.yuta@fujitsu.com>
//  Copyright (c) 2018-2021, Ken-Ichi Ishikawa  <ishikawa@theo.phys.sci.hirosima-u.ac.jp>
//  Copyright (c) 2019-2021, Issaku Kanamori    <kanamori-i@riken.jp>
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
//  supported by MEXT's programs for the Development and Improvement for the Next
//  Generation Ultra High-Speed Computer System, under its Subsidies for Operating the
//  Specific Advanced Large Research Facilities, and Priority Issue 9 
//  (Elucidation of the Fundamental Laws and Evolution of the Universe) to be tackled by
//  using the Supercomputer Fugaku.
//
//****************************************************************************************
/*
  rdma commucation with double buffering:
      wrapper functions for rdma communication
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2018 Dec. 17, 2019 Feb. 6

  Copying and distribution of this file, with or without modification,
  are permitted in any medium without royalty provided the copyright
  notice and this notice are preserved.  This file is offered as-is,
  without any warranty.


  Usage:


void rdma_comlib_2buf::comlib_init()
void rdma_comlib_2buf::comlib_finalize()
  initialize/finalize the rdma communication
  NOTE: Re-initialization does not work in this FJMPI RDMA implementation.

member functions:

void init(const int tni_id, const int dst_rank, const int rcv_rank, const size_t size)
  allocate the buffer
  tni_id:  specifies the Tofu Netowork Interface (TNI), [0...5]
  dst_rank: "target" rank for the 1-to-1 communication
  rcv_rank: rank of which "target" is this rank
  size:    buffer size in byte

void finalize()
  deallocate the buffer

void* sbuff()
  retunrs pointer to the send buffer

void* rbuff()
  recieve buffer with the current parity

void irecv()
  asynchronous RDMA get the data from reciev buffer with the current parity

void isend()
   asynchronous RDMA send/receive with the current parity

void send_wait()
   check and wait for RDMA isend finishes

void recv_wait()
   check and wait for RDMA irecieve finishes

void irecv_ok()
   set ok status to the receive buffer with the current parity

void reset_comm()
   clear the has_started flag set by isend(), in order to supress *_wait()


// interface to rankmap
static int set_rankmap(int dim)
  prepare rankmap and its internal id.
  (in) dim: dimension of the system, [2,4]

static int get_rank_map(int *rank_coord, int *neighbor_ranks, int *rank_size)
  obtain the prepared the rankmap
  (out) rank_coord: 2 or 4 dim array to specify the logical rank cooridanate
  (out) neighbor_ranks: 4 or 8 dim array to spedcify the logical neghbors.  Directions are [+x, -x, +y, -y,...]
  (out) rank_size: logical size of 2 or 4 dim ranks.

 */
#ifndef RDMA_COMLIB_2BUF_H
#define RDMA_COMLIB_2BUF_H
#include <stdlib.h>
#include <mpi.h>


#ifdef _UTOFU_RDMA

#include "rdma_utofu_comlib.h"
extern "C" {
  /*
  struct rdma_comlib_data;

  void rdma_comlib_init(void);
  void rdma_comlib_finalize(void);

  void rdma_comlib_new(rdma_comlib_data *id, const int *tni_id, const int *dst_rank, const int *rcv_rank, const size_t *size);
  void rdma_comlib_delete(rdma_comlib_data *id);

  void rdma_comlib_isend(rdma_comlib_data *id);
  void rdma_comlib_send_wait(rdma_comlib_data *id);

  void rdma_comlib_irecv(rdma_comlib_data *id);
  void rdma_comlib_recv_wait(rdma_comlib_data *id);

  void *get_sbuff_ptr_rdma_comlib_data(rdma_comlib_data *id);
  void *get_rbuff_ptr_rdma_comlib_data(rdma_comlib_data *id);

  void destroy_rdma_comlib_data(rdma_comlib_data *id);
  void *create_rdma_comlib_data();

  int rdma_comlib_get_ssize(const rdma_comlib_data *id);
  void rdma_comlib_swap_vcq_for_sending(rdma_comlib_data *id1, rdma_comlib_data *id2);
  */
  #include "rankmap_lib.h"
}

#endif

class rdma_comlib_2buf{

 public:

  rdma_comlib_2buf(): m_is_allocated(false), m_has_started(false) { }
  ~rdma_comlib_2buf() {
    finalize();  // delete the allocated buffers
  }

  // wrapper to rdma_comlib: init/finalize
  static void comlib_init() {
    rdma_comlib_init();
  }

  static void comlib_finalize(void) {
    rdma_comlib_finalize();
  }

  // todo: separate rank map/tni assignment to indepedent class
  static int set_rankmap(int dim, const int *proc_size=0){
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if(dim==2){
      if(proc_size==0){
        m_rankmap_id=rankmap_lib_set_rankmap2d();
      } else {
        m_rankmap_id=rankmap_lib_set_rankmap2d_size(proc_size);
      }
    } else if(dim==4) {
      if(proc_size==0){
        m_rankmap_id=rankmap_lib_set_rankmap4d();
      } else {
        m_rankmap_id=rankmap_lib_set_rankmap4d_size(proc_size);
      }
    } else {
      fprintf(stderr, "err in making rankmap: unkown dim = %d\n",dim);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if(m_rankmap_id<0){
      if(myrank==0){
        fprintf(stderr, "WARNING: err in making rankmap, using the default map.\n");
      }
    }
    return m_rankmap_id;
  }

  static int get_rank_map(int *rank_coord, int *neighbor_ranks,
                          int *rank_size){
    if(m_rankmap_id>=0) {
      rankmap_lib_get_rankmap(rank_coord, neighbor_ranks, rank_size);
    }
    return m_rankmap_id;
  }


  static void get_tni_ids(int *tni_ids, int m_rankmap_id){
    rankmap_lib_get_tni_list(tni_ids, &m_rankmap_id);
  }

  // wapper to rdma_comlib:  new/delete
  void init(const int tni_id, const int dst_rank, const int rcv_rank, const size_t size);
  void finalize();

  // accesser to the buffers
  void* sbuff() const
  {
    return get_sbuff_ptr_rdma_comlib_data( m_rdma_data );
  }
  void* rbuff() const
  {
    return get_rbuff_ptr_rdma_comlib_data( m_rdma_data );
  }

  // wrapper to rdma_comlib, cont'd
  void irecv()
  {
    rdma_comlib_irecv( m_rdma_data );
  }

  void isend()
  {
    rdma_comlib_isend( m_rdma_data );
    m_has_started = true;
  }

  void send_wait()
  {
    if(m_has_started){
      rdma_comlib_send_wait( m_rdma_data );
    }
    m_has_started = false;
  }

  void recv_wait()
  {
    // wait for the recv buffer is ready
    rdma_comlib_recv_wait( m_rdma_data );
  }

  void irecv_ok()
  {
    // do nothing
  }

  void reset_comm()
  {
    m_has_started=false;
  }

  static void swap_vcq_for_sending(rdma_comlib_2buf &buf1, rdma_comlib_2buf &buf2){
    rdma_comlib_swap_vcq_for_sending(buf1.m_rdma_data, buf2.m_rdma_data);
  };

 private:
  volatile bool m_has_started;
  bool m_is_allocated;
  rdma_comlib_data *m_rdma_data;

 private:
  volatile static int m_rankmap_id;

};


#endif


