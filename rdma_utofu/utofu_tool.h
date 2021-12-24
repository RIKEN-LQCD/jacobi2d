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
#ifndef _UTOFU_TOOL_H
#define _UTOFU_TOOL_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

const char *get_utofu_error_string(int rc)
{
  switch (rc) {
    case UTOFU_SUCCESS :
      return "UTOFU_SUCCESS";

    case UTOFU_ERR_NOT_FOUND :
      return "UTOFU_ERR_NOT_FOUND";
      
    case UTOFU_ERR_NOT_COMPLETED :
      return "UTOFU_ERR_NOT_COMPLETED";
      
    case UTOFU_ERR_NOT_PROCESSED :
      return "UTOFU_ERR_NOT_PROCESSED";
      
    case UTOFU_ERR_BUSY :
      return "UTOFU_ERR_BUSY";
      
    case UTOFU_ERR_USED :
      return "UTOFU_ERR_USED";
      
    case UTOFU_ERR_FULL :
      return "UTOFU_ERR_FULL";
      
    case UTOFU_ERR_NOT_AVAILABLE :
      return "UTOFU_ERR_NOT_AVAILABLE";
      
    case UTOFU_ERR_NOT_SUPPORTED :
      return "UTOFU_ERR_NOT_SUPPORTED";
      

    case UTOFU_ERR_TCQ_OTHER :
      return "UTOFU_ERR_TCQ_OTHER";
      
    case UTOFU_ERR_TCQ_DESC :
      return "UTOFU_ERR_TCQ_DESC";
      
    case UTOFU_ERR_TCQ_MEMORY :
      return "UTOFU_ERR_TCQ_MEMORY";
      
    case UTOFU_ERR_TCQ_STADD :
      return "UTOFU_ERR_TCQ_STADD";
      
    case UTOFU_ERR_TCQ_LENGTH :
      return "UTOFU_ERR_TCQ_LENGTH";
      

    case UTOFU_ERR_MRQ_OTHER :
      return "UTOFU_ERR_MRQ_OTHER";
      
    case UTOFU_ERR_MRQ_PEER :
      return "UTOFU_ERR_MRQ_PEER";
      
    case UTOFU_ERR_MRQ_LCL_MEMORY :
      return "UTOFU_ERR_MRQ_LCL_MEMORY";
      
    case UTOFU_ERR_MRQ_RMT_MEMORY :
      return "UTOFU_ERR_MRQ_RMT_MEMORY";
      
    case UTOFU_ERR_MRQ_LCL_STADD :
      return "UTOFU_ERR_MRQ_LCL_STADD";
      
    case UTOFU_ERR_MRQ_RMT_STADD :
      return "UTOFU_ERR_MRQ_RMT_STADD";
      
    case UTOFU_ERR_MRQ_LCL_LENGTH :
      return "UTOFU_ERR_MRQ_LCL_LENGTH";
      
    case UTOFU_ERR_MRQ_RMT_LENGTH :
      return "UTOFU_ERR_MRQ_RMT_LENGTH";
      

    case UTOFU_ERR_BARRIER_OTHER :
      return "UTOFU_ERR_BARRIER_OTHER";
      
    case UTOFU_ERR_BARRIER_MISMATCH :
      return "UTOFU_ERR_BARRIER_MISMATCH";
      

    case UTOFU_ERR_INVALID_ARG :
      return "UTOFU_ERR_INVALID_ARG";
      
    case UTOFU_ERR_INVALID_POINTER :
      return "UTOFU_ERR_INVALID_POINTER";
      
    case UTOFU_ERR_INVALID_FLAGS :
      return "UTOFU_ERR_INVALID_FLAGS";
      
    case UTOFU_ERR_INVALID_COORDS :
      return "UTOFU_ERR_INVALID_COORDS";
      
    case UTOFU_ERR_INVALID_PATH :
      return "UTOFU_ERR_INVALID_PATH";
      
    case UTOFU_ERR_INVALID_TNI_ID :
      return "UTOFU_ERR_INVALID_TNI_ID";
      
    case UTOFU_ERR_INVALID_CQ_ID :
      return "UTOFU_ERR_INVALID_CQ_ID";
      
    case UTOFU_ERR_INVALID_BG_ID :
      return "UTOFU_ERR_INVALID_BG_ID";
      
    case UTOFU_ERR_INVALID_CMP_ID :
      return "UTOFU_ERR_INVALID_CMP_ID";
      
    case UTOFU_ERR_INVALID_VCQ_HDL :
      return "UTOFU_ERR_INVALID_VCQ_HDL";
      
    case UTOFU_ERR_INVALID_VCQ_ID :
      return "UTOFU_ERR_INVALID_VCQ_ID";
      
    case UTOFU_ERR_INVALID_VBG_ID :
      return "UTOFU_ERR_INVALID_VBG_ID";
      
    case UTOFU_ERR_INVALID_PATH_ID :
      return "UTOFU_ERR_INVALID_PATH_ID";
      
    case UTOFU_ERR_INVALID_STADD :
      return "UTOFU_ERR_INVALID_STADD";
      
    case UTOFU_ERR_INVALID_ADDRESS :
      return "UTOFU_ERR_INVALID_ADDRESS";
      
    case UTOFU_ERR_INVALID_SIZE :
      return "UTOFU_ERR_INVALID_SIZE";
      
    case UTOFU_ERR_INVALID_STAG :
      return "UTOFU_ERR_INVALID_STAG";
      
    case UTOFU_ERR_INVALID_EDATA :
      return "UTOFU_ERR_INVALID_EDATA";
      
    case UTOFU_ERR_INVALID_NUMBER :
      return "UTOFU_ERR_INVALID_NUMBER";
      
    case UTOFU_ERR_INVALID_OP :
      return "UTOFU_ERR_INVALID_OP";
      
    case UTOFU_ERR_INVALID_DESC :
      return "UTOFU_ERR_INVALID_DESC";
      
    case UTOFU_ERR_INVALID_DATA :
      return "UTOFU_ERR_INVALID_DATA";
      

    case UTOFU_ERR_OUT_OF_RESOURCE :
      return "UTOFU_ERR_OUT_OF_RESOURCE";
      
    case UTOFU_ERR_OUT_OF_MEMORY :
      return "UTOFU_ERR_OUT_OF_MEMORY";
      
    case UTOFU_ERR_FATAL :
      return "UTOFU_ERR_FATAL";
      
  }

  return "UTOFU_NO_ERROR_CODE";
}

#define UTOFU_CHECK_ERROR(err)    __utofuCheckError(err, __FILE__, __LINE__ )

#ifdef __cplusplus
inline
#endif
void __utofuCheckError(const int err, const char *file, const int line )
{
  if (err != UTOFU_SUCCESS) {
    fprintf(stderr, "rank %d: Failed at %s : %d : %s\n", m_rdma_comlib_myrank,file,line,get_utofu_error_string(err));
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  return;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif
