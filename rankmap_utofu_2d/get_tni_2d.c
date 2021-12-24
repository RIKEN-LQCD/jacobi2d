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

#include <stdint.h>
#include <stdio.h>

#define TNI_PATTERN 1

#ifdef _USE_TNI_PATTERN2
#undef  TNI_PATTERN
#define TNI_PATTERN 2
#endif

#ifdef _USE_TNI_PATTERN3
#undef  TNI_PATTERN
#define TNI_PATTERN 3
#endif

#ifdef _USE_TNI_PATTERN4
#undef  TNI_PATTERN
#define TNI_PATTERN 4
#endif

int get_tni_list(int *tni,
                 const int myrank, const int *nranks,
                 const int dim, const int flag){
  int rank_in_node = myrank % 4;

#if TNI_PATTERN == 1
  if(myrank==0){
    printf("using TNI_PATTERN 1 (default)\n");
  }

  // from +x direction  [i.e., to -x]
  int dir=0;
  if(rank_in_node == 0 || rank_in_node == 2){
    tni[dir]=1;
  } else {
    tni[dir]=4;
  }

  // from -x direction  [i.e., to +x]
  dir=1;
  if(rank_in_node == 0 || rank_in_node == 2){
    tni[dir]=4;
  } else {
    tni[dir]=0;
  }

  // from +y direction [i.e., to -y]
  dir=2;
  if(rank_in_node == 0 || rank_in_node == 1){
    tni[dir]=3;
  } else {
    tni[dir]=5;
  }

  // from -y direction [i.e., to +y]
  dir=3;
  if(rank_in_node == 0 || rank_in_node == 1){
    tni[dir]=5;
    } else {
    tni[dir]=2;
  }
#elif  TNI_PATTERN == 2
{
  if(myrank==0){
    printf("using NTI_PATTERN 2 (defined _USE_TNI_PATTERN2\n");
  }
  switch(rank_in_node){
  case 0:  // (0,0)
    tni[0] = 0; // to -x
    tni[1] = 2; // to +x
    tni[2] = 1; // to -y
    tni[3] = 4; // to +y
    break;
  case 1:  // (1,0)
    tni[0] = 2; // to -x
    tni[1] = 0; // to +x
    tni[2] = 1; // to -y
    tni[3] = 5; // to +y
    break;
  case 2:  // (0,1)
    tni[0] = 0; // to -x
    tni[1] = 3; // to +x
    tni[2] = 4; // to -y
    tni[3] = 1; // to +y
    break;
  case 3:  // (1,1)
    tni[0] = 3; // to -x
    tni[1] = 0; // to +x
    tni[2] = 5; // to -y
    tni[3] = 1; // to +y
    break;
  default:
    fprintf(stderr, "rank %d: cannot happen in get_tni_2d [for PATTERN 2]: rank_in_node=%d\n",
	    myrank, rank_in_node);
    return -1;
  }
}

#elif  TNI_PATTERN == 3
{
  if(myrank==0){
    printf("using NTI_PATTERN 3 (defined _USE_TNI_PATTERN3\n");
  }
  switch(rank_in_node){
  case 0:  // (0,0)
    tni[0] = 0; // to -x
    tni[1] = 1; // to +x
    tni[2] = 2; // to -y
    tni[3] = 3; // to +y
    break;
  case 1:  // (1,0)
    tni[0] = 4; // to -x
    tni[1] = 5; // to +x
    tni[2] = 0; // to -y
    tni[3] = 1; // to +y
    break;
  case 2:  // (0,1)
    tni[0] = 2; // to -x
    tni[1] = 3; // to +x
    tni[2] = 4; // to -y
    tni[3] = 5; // to +y
    break;
  case 3:  // (1,1)
    tni[0] = 0; // to -x
    tni[1] = 1; // to +x
    tni[2] = 2; // to -y
    tni[3] = 3; // to +y
    break;
  default:
    fprintf(stderr, "rank %d: cannot happen in get_tni_2d [for PATTERN 3]: rank_in_node=%d\n",
	    myrank, rank_in_node);
    return -1;
  }
}
#elif TNI_PATTERN == 4
  if(myrank==0){
    printf("using TNI_PATTERN 4\n");
  }

  // from +x direction  [i.e., to -x]
  int dir=0;
  if(rank_in_node == 0 || rank_in_node == 2){
    tni[dir]=4;
  } else {
    tni[dir]=0;
  }

  // from -x direction  [i.e., to +x]
  dir=1;
  if(rank_in_node == 0 || rank_in_node == 2){
    tni[dir]=1;
  } else {
    tni[dir]=4;
  }

  // from +y direction [i.e., to -y]
  dir=2;
  if(rank_in_node == 0 || rank_in_node == 1){
    tni[dir]=5;
  } else {
    tni[dir]=2;
  }

  // from -y direction [i.e., to +y]
  dir=3;
  if(rank_in_node == 0 || rank_in_node == 1){
    tni[dir]=3;
    } else {
    tni[dir]=5;
  }

#else
  #error unknown TNI_PATTERN
#endif
  return 0;
}


int get_tni_list_default(int *tni,
                         const int myrank, const int *nranks,
                         const int dim, const int flag){

  return get_tni_list(tni, myrank, nranks, dim, flag);
}
