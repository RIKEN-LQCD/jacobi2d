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
      wrapper for functions rdma communication
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2018 Dec. 17, 2019 Feb. 6

  UTOFU version
  Issaku Kanamori kanamori-i@riken.jp
  2019 Nov. 29

  Copying and distribution of this file, with or without modification,
  are permitted in any medium without royalty provided the copyright
  notice and this notice are preserved.  This file is offered as-is,
  without any warranty.

 */


#ifdef _UTOFU_RDMA
#include "rdma_comlib_2buf.h"


volatile int rdma_comlib_2buf::m_rankmap_id=-999;

// initializatoin of the buffer
void rdma_comlib_2buf::init(const int tni_id, const int dst_rank, const int rcv_rank, const size_t size)
{

  if(m_is_allocated){  // free the buffer if allocated
    finalize();
  }
  // create and initialize the data for th rdma comm.
  m_rdma_data=(rdma_comlib_data*)create_rdma_comlib_data();
  rdma_comlib_new(m_rdma_data, &tni_id, &dst_rank, &rcv_rank, &size);

  // flag to check if the commucation has started
  m_has_started=false;

  // set the flag
  m_is_allocated=true;


}

// cleaning up the buffer
void rdma_comlib_2buf::finalize()
{
  if(!m_is_allocated){ // check if the buffer is allocated
    return;
  }

  destroy_rdma_comlib_data(m_rdma_data);
  // clear the flag
  m_is_allocated=false;
}


#endif
