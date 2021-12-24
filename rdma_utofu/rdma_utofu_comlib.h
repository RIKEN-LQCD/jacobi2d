//****************************************************************************************
//
//  Copyright (c) 2015-2021, Yoshifumi Nakamura <nakamura@riken.jp>
//  Copyright (c) 2015-2021, Yuta Mukai         <mukai.yuta@fujitsu.com>
//  Copyright (c) 2018-2021, Ken-Ichi Ishikawa  <ishikawa@theo.phys.sci.hirosima-u.ac.jp>
//  Copyright (c) 2019-2021, Issaku Kanamori    <kanamori-i@riken.jp>
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
#ifndef _RDMA_UTOFU_COMLIB_H
#define _RDMA_UTOFU_COMLIB_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <utofu.h>

typedef struct {
  utofu_vcq_hdl_t local_vcq_hdl; // handle to the local virtual control queue
  utofu_vcq_id_t  local_vcq_id;  // id of the local virtual control queue     (recv-put target)
  utofu_vcq_id_t target_vcq_id;  // id of the hadle to the remote virtual control queue (send-put target)
  int local_rank;                // local MPI rank
  int target_rank;               // target MPI rank (send to)
  int from_rank;                 // MPI rank from which data comes
  int tni_id;                    // TNI ID
  int dummy0;                    // padding
  int dummy1;
  int dummy2;
} rdma_comlib_communicator;

typedef struct {
  volatile int   *rbuff;
  utofu_stadd_t  target_stadd_rbuff;
  int            last_component;
} rbuff_t;

typedef struct {
  volatile int *sbuff;           // pointer to int array (additional +1 component is required to recv check)
  volatile int *rbuff;           // pointer to int array
  volatile int *rbuff_all;       // pointer to int double buffering array
  size_t length;                 // size of array (in byte data part)
  size_t data_size;              // size in byte (data + additilanl commponent),
  size_t data_len_int;           // the same as data size but in sizeof(int) // not used
  size_t data_size_all;          // size in byte (data + additilanl commponent)*2, 256 byte aligned (fugaku cache line) size of rbuff_all.
  size_t rbuff_ptr_disp;         // pointer displacement in rbuff_all with parity (byte), 256 byte aligned
  int    parity;                 // double buffering parity
  rbuff_t rbuff_wp[2];           // double receive buffer 

  rdma_comlib_communicator comm;
  utofu_stadd_t  sending_stadd;         // RDMA address (stadd) of local send buffer
  utofu_stadd_t   local_stadd_rbuff;    // RDMA address (stadd) of local reciev buffer
  utofu_stadd_t  target_stadd_rbuff;    // RDMA address (stadd) of target receive buffer
  utofu_vcq_hdl_t sending_vcq_hdl;      // vcq handle used for data put == local_vcd_hdl
  utofu_stadd_t   local_stadd_sbuff;    // stadd of send buffer  == sending_stadd
  int tag;                              // message tag
  int last_component;                   // checking the data arrival // not used (instead rbuff_wp[0:1].last_component is used)
  FILE *fp;                             // debugging output
} rdma_comlib_data;

void rdma_comlib_init(void);
void rdma_comlib_finalize(void);

void rdma_comlib_new(rdma_comlib_data *id, const int *tni_id, const int *dst_rank, const int *rcv_rank, const size_t *size);
void rdma_comlib_delete(rdma_comlib_data *id);

void rdma_comlib_irecv(rdma_comlib_data *id);
void rdma_comlib_isend(rdma_comlib_data *id);
void rdma_comlib_send_wait(rdma_comlib_data *id);
void rdma_comlib_recv_wait(rdma_comlib_data *id);

void *create_rdma_comlib_data();
void destroy_rdma_comlib_data(rdma_comlib_data *id);

void *get_sbuff_ptr_rdma_comlib_data(rdma_comlib_data *id);
void *get_rbuff_ptr_rdma_comlib_data(rdma_comlib_data *id);

void rdma_comlib_swap_vcq_for_sending(rdma_comlib_data *id1, rdma_comlib_data *id2);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
