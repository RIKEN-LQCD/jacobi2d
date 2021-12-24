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
////////////////////////////////////////////////////////////////////////////////
//
//  The original version is by Ken-Ichi Ishikawa <ishikawa@theo.phys.sci.hiroshima-u.ac.jp>
//
//  2018 Dec. 17:
//    modification of the original version
//    allows the user allocate the buffers in the outside
//    Issaku Kanamori <kanamori@hiroshima-u.ac.jp>
//
//  2019 Nov. 28:
//    UTOFU version.  Removed interfaces which are not used anymore.
//
//  2020 Jan. 9:
//    added rdma_comlib_new_duplicate
//
//  2021 Apr. 2:
//    lddhmc version:
//       sbuff/rbuff are allocated in rdma_comlib_new. No user allocation is permitted.
//       vcq pool is used. swap, duplicate functions are removed.
//       double buffering is involved.
//
//  2021 Sep. 6:
//    recovered rdma_comlib_swap_vcq_for_sending
//
//  Copying and distribution of this file, with or without modification,
//  are permitted in any medium without royalty provided the copyright
//  notice and this notice are preserved.  This file is offered as-is,
//  without any warranty.
//

#define _UTOFU_RDMA
#ifdef _UTOFU_RDMA

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <utofu.h>

//#define _RDMA_DEBUG
#define RDMA_NO_REMOTE_MRQ_POLLING
#define RDMA_USE_CACHE_INJECTION


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef __cplusplus
const int MAX_VCQ = 7;                   // Maximum number of VCQ for each TNI
const int MAXNUM_TNIID = 6;              // Maximum number of TNI ID
#else
#define MAX_VCQ       (7)
#define MAXNUM_TNIID  (6)
#endif

const int FUGAKU_ALIGN = 256;            // memory arrignment for send/recv buffer
const int MAX_RDMA_DATASIZE = 16777212;  // RDMA put/get MAX data size in byte
const int MAXNUM_TAG   = 255;            // RDMA put/get MAX number of message tag
const int TAG_OFFSET   = 128;            // offset for the message tag
const int MAX_LAST_COMPONENT  = 256;     // maximum value of the watchdog data for receive data polling. (watch_dog value is in [0...255])

static int m_rdma_comlib_is_initialized = 0;
static int m_rdma_comlib_myrank = 0;     // MPI LOCAL RANK is stored.
static int m_rdma_comlib_tag = 128;      // counts RDMA message tag. [0..127] + TAG_OFFSET
static size_t m_num_tnis = 0;

#include "utofu_tool.h"
#include "rdma_utofu_comlib.h"

typedef struct {
  utofu_vcq_hdl_t vcq_hdl;
  utofu_vcq_id_t  vcq_id;
} local_vcq_t;

typedef struct {
  local_vcq_t local_vcq[MAX_VCQ];
  int         ref_count[MAX_VCQ];
  int         tni_id;
  int         top_id;
  int         num_allocated;
} local_vcq_pool_t;

static local_vcq_pool_t m_local_vcq_pool[MAXNUM_TNIID];

void *create_rdma_comlib_data()
{
  rdma_comlib_data* id = (rdma_comlib_data *)malloc(sizeof(rdma_comlib_data));
  id->sbuff = NULL;
  id->rbuff = NULL;
  id->length = 0;
  id->data_size = 0;
  id->data_len_int = 0;
  return (void*)id;
}

void destroy_rdma_comlib_data(rdma_comlib_data *id)
{
  rdma_comlib_delete(id);
  free(id);
  id = NULL;
}

void *get_sbuff_ptr_rdma_comlib_data(rdma_comlib_data *id)
{
  return (void *)(id->sbuff);
}

void *get_rbuff_ptr_rdma_comlib_data(rdma_comlib_data *id)
{
  return (void *)(id->rbuff);
}


int rdma_comlib_get_ssize(const rdma_comlib_data *id)
{
   return (*id).length;
}

void init_local_vcq_pool()
{
  //
  // initialize vcq pool for each tni
  //
  for (int tni_id = 0; tni_id < m_num_tnis; ++tni_id) {

    m_local_vcq_pool[tni_id].tni_id = tni_id;
    m_local_vcq_pool[tni_id].top_id = 0;
    m_local_vcq_pool[tni_id].num_allocated = 0;
  }
}

void extend_local_vcq_pool(int tni_id){
  const unsigned long int vcq_flag = UTOFU_VCQ_FLAG_THREAD_SAFE;
  int vcq_id = m_local_vcq_pool[tni_id].num_allocated;
  m_local_vcq_pool[tni_id].ref_count[vcq_id] = 0;

  int rc = utofu_create_vcq(tni_id, vcq_flag, &(m_local_vcq_pool[tni_id].local_vcq[vcq_id].vcq_hdl));
  UTOFU_CHECK_ERROR(rc);

  rc = utofu_query_vcq_id(m_local_vcq_pool[tni_id].local_vcq[vcq_id].vcq_hdl,&(m_local_vcq_pool[tni_id].local_vcq[vcq_id].vcq_id));
  UTOFU_CHECK_ERROR(rc);

  ++m_local_vcq_pool[tni_id].num_allocated;
}

void delete_local_vcq_pool()
{
  for (int tni_id = 0; tni_id < m_num_tnis; ++tni_id) {

     m_local_vcq_pool[tni_id].tni_id = tni_id;
     m_local_vcq_pool[tni_id].top_id = 0;

     for (int vcq_id = 0; vcq_id < m_local_vcq_pool[tni_id].num_allocated; ++vcq_id) {

       int rc = utofu_free_vcq(m_local_vcq_pool[tni_id].local_vcq[vcq_id].vcq_hdl);
       UTOFU_CHECK_ERROR(rc);
       m_local_vcq_pool[tni_id].ref_count[vcq_id] = 0;

     }
     m_local_vcq_pool[tni_id].num_allocated=0;
  }
}

void get_local_vcq_hdl_and_id(const int tni_id, utofu_vcq_hdl_t *vcq_hdl, utofu_vcq_id_t *vcq_id)
{
  //
  // get local_vcq_hdl and local_vcq_id from vcq_pool
  //
  //   tni_id : TNI id number
  //  vcq_hdl : vcq_hdl
  //  vcq_id  : vcq_id
  //
  if(m_local_vcq_pool[tni_id].num_allocated < MAX_VCQ){
    extend_local_vcq_pool(tni_id); // extend the pool up to MAX_VCQ
  }

  int top_id = m_local_vcq_pool[tni_id].top_id;
  *vcq_hdl = m_local_vcq_pool[tni_id].local_vcq[top_id].vcq_hdl;
  *vcq_id  = m_local_vcq_pool[tni_id].local_vcq[top_id].vcq_id;

  m_local_vcq_pool[tni_id].ref_count[top_id]++;  // increment reference counter of the vcq
  m_local_vcq_pool[tni_id].top_id = (m_local_vcq_pool[tni_id].top_id + 1) % MAX_VCQ; // move top pointer (cyclic)
}


void rdma_comlib_init(void)
{
//
// Initialize UTOFU communication
//
// assumes MPI is already initialized
  if (0 == m_rdma_comlib_is_initialized) {

    MPI_Comm_rank(MPI_COMM_WORLD,&m_rdma_comlib_myrank);

    if ( 0 == m_rdma_comlib_myrank ) {
      printf("%%%% UTOFU interface for one sided communcation  (RDMA) is used.\n");
      printf("sizeof(rdma_comlib_data) = %d, sizeof(int*) = %d\n", sizeof(rdma_comlib_data), sizeof(int*));
#ifdef RDMA_NO_REMOTE_MRQ_POLLING
      printf("RDMA_NO_REMOTE_MRQ_POLLING is defined.\n");
#endif
      fflush(stdout);
    }

    // check if we can use full TNIs
    utofu_tni_id_t *tni_ids;  // array of TNI IDs
    if (UTOFU_SUCCESS != utofu_get_onesided_tnis(&tni_ids, &m_num_tnis)) {
      fprintf(stderr, "rank %d: Failed at utofu_get_onesided_tnis()!\n", m_rdma_comlib_myrank);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (int id_tni = 0; id_tni < m_num_tnis; ++id_tni) {
      // obtain available features and dump them
      struct utofu_onesided_caps *onesided_caps;
      int rc = utofu_query_onesided_caps(tni_ids[id_tni], &onesided_caps);
      UTOFU_CHECK_ERROR(rc);

#ifdef _RDMA_DEBUG
      // dump caps...
      if (m_rdma_comlib_myrank==0) {
        int rank=m_rdma_comlib_myrank;
        fprintf(stderr, "rank %d: tni_ids[%d] = %lu\n", rank, id_tni, tni_ids[id_tni]);
        fprintf(stderr, "rank %d: utofu_onesided_caps: flags = %lu\n", rank, onesided_caps->flags);
        fprintf(stderr, "rank %d: utofu_onesided_caps: arwm_ops  = %lu\n", rank, onesided_caps->armw_ops);
        fprintf(stderr, "rank %d: utofu_onesided_caps: num_cmp_ids  = %u\n", rank, onesided_caps->num_cmp_ids);
        fprintf(stderr, "rank %d: utofu_onesided_caps: num_reserved_stags  = %u\n", rank, onesided_caps->num_reserved_stags);
        fprintf(stderr, "rank %d: utofu_onesided_caps: cache_line_size  = %zu\n", rank, onesided_caps->cache_line_size);
        fprintf(stderr, "rank %d: utofu_onesided_caps: stag_address_alignment  = %zu\n", rank, onesided_caps->stag_address_alignment);
        fprintf(stderr, "rank %d: utofu_onesided_caps: max_toq_desc_size  = %zu\n", rank, onesided_caps->max_toq_desc_size);
        fprintf(stderr, "rank %d: utofu_onesided_caps: max_putget_size  = %zu\n", rank, onesided_caps->max_putget_size);
        fprintf(stderr, "rank %d: utofu_onesided_caps: max_piggyback_size  = %zu\n", rank, onesided_caps->max_piggyback_size);
        fprintf(stderr, "rank %d: utofu_onesided_caps: max_edata_size  = %zu\n", rank, onesided_caps->max_edata_size);
        fprintf(stderr, "rank %d: utofu_onesided_caps: max_mtu  = %zu\n", rank, onesided_caps->max_mtu);
        fprintf(stderr, "rank %d: utofu_onesided_caps: max_gap  = %zu\n", rank, onesided_caps->max_gap);
      }
#endif
    }

    if (m_num_tnis < MAXNUM_TNIID) {
      fprintf(stderr, "rank %d: only %d TNIs are available, aborting...\n", m_rdma_comlib_myrank, m_num_tnis);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    free(tni_ids);

    //
    // initialize local vcq pool
    //
    init_local_vcq_pool();

    m_rdma_comlib_is_initialized = 1;

  }
}


void rdma_comlib_finalize(void)
{
  if (1 == m_rdma_comlib_is_initialized) {
    delete_local_vcq_pool();
    m_rdma_comlib_is_initialized = 0;
  }
}

void rdma_comlib_communicator_new(rdma_comlib_communicator *comm, const int *tni_id, const int *dst_rank, const int *rcv_rank)
{
//
// Assign local vcq handle and id to communicator form vcq pool and set remote target vcq id
//

  int err = 0;
  if ( (MAXNUM_TNIID <= *tni_id) || ( *tni_id < 0 ) ) {
    fprintf(stderr,"TNI ID should be greater than -1 and less than MAXNUM_TNIID. 0 <= tni_id = %d < %d = MAXNUM_TNIID\n",tni_id,MAXNUM_TNIID);
    err = 1;
  }
  if (err) {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  //
  // set values
  //
  comm->tni_id      = *tni_id;
  comm->target_rank = *dst_rank;
  comm->from_rank   = *rcv_rank;
  comm->local_rank  = m_rdma_comlib_myrank;

  //
  // get local vcq handle and id from vcq pool
  //
  get_local_vcq_hdl_and_id(comm->tni_id, &(comm->local_vcq_hdl), &(comm->local_vcq_id));

  //
  // send local_vcq_id to from_rank and recev target_vcq_id from target_rank
  //
  // from_rank  ->   local_rank   -> target_rank
  //            <-  local_vcq_id  <- target_vcq_id
  //
  uint64_t tmp_sbuff;
  uint64_t tmp_rbuff;
  int tag = ((m_rdma_comlib_tag-TAG_OFFSET) % (MAXNUM_TAG-TAG_OFFSET)) + TAG_OFFSET;
  tmp_sbuff = comm->local_vcq_id;
#ifdef _RDMA_DEBUG
  fprintf(stderr, "rank %d: calling MPI_Sendrecv: tmp_sbuff=%lu\n", m_rdma_comlib_myrank, tmp_sbuff);
#endif

  MPI_Sendrecv(&tmp_sbuff, 1, MPI_UINT64_T, comm->from_rank,   tag,
               &tmp_rbuff, 1, MPI_UINT64_T, comm->target_rank, tag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  comm->target_vcq_id = tmp_rbuff;

  //
  // set the default path to the target rank vcq
  //
  int rc = utofu_set_vcq_id_path(&(comm->target_vcq_id), NULL);
  UTOFU_CHECK_ERROR(rc);

}

void rdma_comlib_new_impl(rdma_comlib_data *id, const size_t *size)
{
//
// Initialize comlib ID, local-send/recv buffer, data size (in byte)
//
//   *id         : comlib_id
//                 user put/get data in id->sbuff : send buffef, id->rbuff : recv buffef
//   size        : data size to be send/recv'd in byte unit.
//

#ifdef _RDMA_DEBUG
  id->fp = stderr;
  fprintf(id->fp,"rank %d: rdma_comlib_new_impl start.\n",m_rdma_comlib_myrank);
#endif

  //
  // compute send/receive data size including polling buffer.
  //
  id->length = *size; // data length in byte

  id->data_size = (id->length + sizeof(int));  // additional last 4-byte component is used for polling data receive.

  if ( MAX_RDMA_DATASIZE < id->data_size ) {
    fprintf(stderr,"RDMA data size is too large.  data size should be data_size = %d < %d = MAX_RDMA_DATASIZE\n",id->data_size,MAX_RDMA_DATASIZE);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if ( (id->data_size % (sizeof(int))) != 0 ) {
    fprintf(stderr,"RDMA data size shoule be a multiple of sizeof(int)=%d bytes.  data_size = %d\n",sizeof(int),id->data_size);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  id->data_len_int = (id->data_size) / (sizeof(int));

  //
  // allocate local send and local receive buffers (cache aligned)
  //
  id->rbuff_ptr_disp = (((id->data_size)/FUGAKU_ALIGN) + 1)*FUGAKU_ALIGN;
  id->data_size_all  = id->rbuff_ptr_disp*2;
  posix_memalign((void **)(&(id->sbuff)),     FUGAKU_ALIGN, id->data_size);
  posix_memalign((void **)(&(id->rbuff_all)), FUGAKU_ALIGN, id->data_size_all); // doubel buffering (two rbuffs)

  //
  // set RDMA message tag
  //
  id->tag = ((m_rdma_comlib_tag-TAG_OFFSET + MAXNUM_TAG) % (MAXNUM_TAG-TAG_OFFSET)) + TAG_OFFSET;
  m_rdma_comlib_tag++;

  //
  // register local send buffer
  //
  int rc;
  rc = utofu_reg_mem(id->comm.local_vcq_hdl, (void *)(id->sbuff), id->data_size, 0, &id->local_stadd_sbuff);
  UTOFU_CHECK_ERROR(rc);
  id->sending_vcq_hdl = id->comm.local_vcq_hdl;
  id->sending_stadd   = id->local_stadd_sbuff;

  //
  // register local receive buffer
  //
  rc = utofu_reg_mem(id->comm.local_vcq_hdl, (void *)(id->rbuff_all), id->data_size_all, 0, &id->local_stadd_rbuff);
  UTOFU_CHECK_ERROR(rc);

  //
  // send local_stadd(rbuff) to from_rank and recv target_stadd(rbuff) from target_rank
  //
  //  from_rank  ->        local_rank  ->  target_rank
  //             <- local_stadd(rbuff) <-  target_stadd(rbuff)
  //
  uint64_t tmp_sbuff;
  uint64_t tmp_rbuff;
  tmp_sbuff = id->local_stadd_rbuff;
  MPI_Sendrecv(&tmp_sbuff, 1, MPI_UINT64_T, id->comm.from_rank,   id->tag,
               &tmp_rbuff, 1, MPI_UINT64_T, id->comm.target_rank, id->tag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  id->target_stadd_rbuff = tmp_rbuff;
  //
  // clear watchdog data in the additional last component
  //
  memset((void *)(id->sbuff)    ,0,id->data_size);
  memset((void *)(id->rbuff_all),0,id->data_size_all);
  id->sbuff[id->data_len_int-1] = 0;
  for (int ip = 0; ip < 2; ++ip) {
    id->rbuff_wp[ip].rbuff              = id->rbuff_all          + (id->rbuff_ptr_disp/sizeof(int))*ip;
    id->rbuff_wp[ip].target_stadd_rbuff = id->target_stadd_rbuff +  id->rbuff_ptr_disp*ip;
    id->rbuff_wp[ip].rbuff[id->data_len_int-1] = MAX_LAST_COMPONENT-1;
    id->rbuff_wp[ip].last_component            = MAX_LAST_COMPONENT-1;
  }

  //
  // set parity and rbuff
  //
  id->parity = 0;
  id->rbuff = id->rbuff_wp[id->parity].rbuff;

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"message length=%ld  size=%ld  len_int=%ld\n", id->length, id->data_size, id->data_len_int); fflush(NULL);
  fprintf(id->fp,"local_rank = %3d, target_rank = %3d ",id->comm.local_rank, id->comm.target_rank);
  fprintf(id->fp,"local_stadd_sbuff = %ld, target_stadd_rbuff = %ld\n",
                 id->local_stadd_sbuff,id->target_stadd_rbuff);
  fprintf(id->fp,"rank %d: rdma_comlib_new OK.\n",m_rdma_comlib_myrank); fflush(NULL);
  fflush(NULL);
#endif

}

void rdma_comlib_new(rdma_comlib_data *id, const int *tni_id, const int *dst_rank, const int *rcv_rank, const size_t *size)
{
//
// Initialize comlib ID, set destination send-to-rank and receive recv-from-rank, local-send/recv buffer, and data size (in byte)
// this version does not malloc,  the resource must be managed by the user
//
//   *id         : comlib_id
//                 user put/get data in id->sbuff : send buffef, id->rbuff : recv buffef
//   *nic_id     : TNI ID [0...5]
//   *dst_rank   : destination MPI rank
//   *rcv_rank   : MPI rank from which data comes
//   *size       : data size to be send/recv'd in byte unit.
//
#ifdef _RDMA_DEBUG
  id->fp = stderr;
  fprintf(id->fp,"rank %d: rdma_comlib_new start.\n",m_rdma_comlib_myrank);
#endif

  id->length = 0;
  if ( 0 == *size ) {
    id->length = *size;
    return;
  }

  //
  // check buffer association
  //
  {
    int err = 0;
    if (id->sbuff != NULL) {
      fprintf(stderr, "sbuff is allocated\n");
      err = 1;
    };
    if (id->rbuff != NULL) {
      fprintf(stderr, "rbuff is allocated\n");
      err = 1;
    };
    if (err) {
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }

  //
  // create vcq
  //
  rdma_comlib_communicator_new(&id->comm, tni_id, dst_rank, rcv_rank);

  //
  // register send/recv buffers and data size
  //
  rdma_comlib_new_impl(id, size);

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_new_ext, done.: id->comm=%p\n",m_rdma_comlib_myrank, id->comm);
#endif

}

void rdma_comlib_delete(rdma_comlib_data *id)
{
//
// Delete comlib ID
//

  if ( 0 == id->length ) return;

  utofu_dereg_mem(id->comm.local_vcq_hdl, id->local_stadd_sbuff, 0);
  utofu_dereg_mem(id->comm.local_vcq_hdl, id->local_stadd_rbuff, 0);

  free((void *)(id->sbuff));
  free((void *)(id->rbuff_all));
  for (int ip = 0; ip < 2; ++ip) {
    id->rbuff_wp[ip].rbuff = NULL;
    id->rbuff_wp[ip].target_stadd_rbuff = 0;
    id->rbuff_wp[ip].last_component = 0;
  }
  id->sbuff = NULL;
  id->rbuff = NULL;
  id->rbuff_all = NULL;

  // cleaning
  id->length = 0;
  id->data_size = 0;
  id->data_size_all = 0;
  id->data_len_int = 0;
  id->tag = 0;
  id->parity = 0;

  id->sending_stadd = 0;
  id->local_stadd_sbuff = 0;
  id->local_stadd_rbuff = 0;
  id->target_stadd_rbuff = 0;
  id->last_component = 0;

  id->sending_vcq_hdl = 0;
}

void rdma_comlib_clear_mrq(rdma_comlib_data *id)
{
//
// MRQ polling to clear MRQ avoiding MRQ overflow.
//
  struct utofu_mrq_notice notice;
  int rc = 0;
  do {
    rc = utofu_poll_mrq(id->comm.local_vcq_hdl, 0UL, &notice);
  } while ( rc != UTOFU_ERR_NOT_FOUND ) ;
#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: rdma_comlib_clear_mrq done.\n",m_rdma_comlib_myrank);
  fflush(id->fp);
#endif
}

void rdma_comlib_irecv(rdma_comlib_data *id)
{
  if ( 0 == id->length ) return;
  //  rdma_comlib_clear_mrq(id);
}

void rdma_comlib_isend(rdma_comlib_data *id)
{
//
// asynchronous RDMA send/receive start with comlib ID
//

  if ( 0 == id->length ) return;

  //
  // reset watchdog data for send.
  //
  id->sbuff[id->data_len_int-1] = (id->sbuff[id->data_len_int-1]+1) % MAX_LAST_COMPONENT;

  //
  // Start data send
  //
#ifdef RDMA_NO_REMOTE_MRQ_POLLING
  const unsigned long int mrq_flag = 0UL;
#else
  const unsigned long int mrq_flag = UTOFU_ONESIDED_FLAG_REMOTE_MRQ_NOTICE;
#endif
#ifdef RDMA_USE_CACHE_INJECTION
  const unsigned long int cache_injection_flag = UTOFU_ONESIDED_FLAG_CACHE_INJECTION;
#else
  const unsigned long int cache_injection_flag = 0UL;
#endif
  const unsigned long int send_flags = 
    UTOFU_ONESIDED_FLAG_TCQ_NOTICE | UTOFU_ONESIDED_FLAG_STRONG_ORDER | mrq_flag | cache_injection_flag;

  uint64_t  edata   = 0; // for mrq polling; the value is not used
  uintptr_t cbvalue = 0; // for tcq polling; the value is not used

  int rc = 0;
  const utofu_stadd_t rmt_stadd = id->rbuff_wp[id->parity].target_stadd_rbuff;
  rc = utofu_put(id->sending_vcq_hdl, id->comm.target_vcq_id, id->sending_stadd,
                 rmt_stadd, id->data_size, edata, send_flags, (void *)cbvalue);
  UTOFU_CHECK_ERROR(rc);

}

void rdma_comlib_recv_wait(rdma_comlib_data *id)
{
//
// check/wait for RDMA irecv finish.
//

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: %s start.\n",m_rdma_comlib_myrank,__func__); fflush(NULL);
#endif

  if ( 0 == id->length ) return;

  rdma_comlib_clear_mrq(id);

  //
  // polling watch dog data in the last element of the receive buffer (strong order put must be used)
  //
  const int watch_dog = id->rbuff_wp[id->parity].last_component;
  while( id->rbuff_wp[id->parity].rbuff[id->data_len_int-1] == watch_dog ){ };

  id->rbuff_wp[id->parity].last_component = id->rbuff_wp[id->parity].rbuff[id->data_len_int-1];  // update the watch dog data

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank: %d, local_vcq_hdl=%p: new id->last_compoient_id=%d\n",
          m_rdma_comlib_myrank, id->comm.local_vcq_hdl, id->last_component, id->rbuff[id->data_len_int-1]);
  fprintf(id->fp,"rank %d: %s OK.\n",m_rdma_comlib_myrank,__func__); fflush(NULL);
#endif

}

void rdma_comlib_send_wait(rdma_comlib_data *id)
{
//
// check/wait for RDMA isend finish.
//

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: %s.\n",m_rdma_comlib_myrank,__func__); fflush(NULL);
#endif

  if ( 0 == id->length ) return;

  //
  // tcq polling and clear tcq avoiding TCQ busy
  //
  void *cbdata;
  int rc;
  do {
    rc = utofu_poll_tcq(id->sending_vcq_hdl, 0, &cbdata);
  } while (rc == UTOFU_ERR_NOT_FOUND);

  //
  // check if the polling is successfully finished
  //
  if (rc != UTOFU_SUCCESS) {
    fprintf(stderr,"rank %d: %s at %d : %s, ERROR: %d, target_rank=%d\n",
            m_rdma_comlib_myrank, __FILE__, __LINE__,__func__, rc, id->comm.target_rank);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  //
  // flip parity and receive buffer
  //
  id->parity = (id->parity + 1) % 2;
  id->rbuff = id->rbuff_wp[id->parity].rbuff;

#ifdef _RDMA_DEBUG
  fprintf(id->fp,"rank %d: %s OK.\n",m_rdma_comlib_myrank,__func__); fflush(NULL);
#endif
}


void rdma_comlib_swap_vcq_for_sending(rdma_comlib_data *id1, rdma_comlib_data *id2){

  if(id1->length != id2->length){
    fprintf(stderr, "wrong lengths in rdma_comlib_swap_vcq_for_sending\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // tell actuall vcq and stadd
  id1->sending_vcq_hdl=id2->comm.local_vcq_hdl;
  id2->sending_vcq_hdl=id1->comm.local_vcq_hdl;

  id1->sending_stadd=id2->local_stadd_sbuff;
  id2->sending_stadd=id1->local_stadd_sbuff;

  // Note that sbuffs are also swapped.
  volatile int *tmp_sbuff=id1->sbuff;
  id1->sbuff=id2->sbuff;
  id2->sbuff=tmp_sbuff;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif
