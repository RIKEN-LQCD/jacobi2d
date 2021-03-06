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

#ifndef RANKMAP_LIB_H
#define RANKMAP_LIB_H

int rankmap_lib_set_rankmap2d();
int rankmap_lib_set_rankmap2d_size(const int *proc_size);
int rankmap_lib_set_rankmap4d();
int rankmap_lib_set_rankmap4d_size(const int *proc_size);
int rankmap_lib_get_rankmap(int *myrand_coord, int *neighbors, int *sizes);
int rankmap_lib_get_tni_list(int *tni_list, const int *flag);

/*
int rankmap_lib_set_rankmap2d();
int rankmap_lib_set_rankmap4d();
  prepare rankmap for 2 or 4 dimensional lattice.
  return value:
     gives the internal id for obtained map.
     < 0 if a proper map is not found.

int rankmap_lib_set_rankmap2d_size(const int *proc_size);
int rankmap_lib_set_rankmap4d_size(const int *proc_size);
  the same as rankmap_lib_set_rankmap2d() and rankmap_lib_set_rankmap4d(),
  but uses the given proc_size



int rankmap_lib_get_rankmap(int *myrand_coord, int *neighbors, int *sizes);
  obtain the rankmap.
  (out) myrank_coord: 2 or 4 dim array to specify the logical rank cooridanate
  (out) neighbors: 4 or 8 dim array to spedcify the logical neghbors.  Directions are [+x, -x, +y, -y,...]
  (out) sizes: logical size of 2 or 4 dim ranks.

int rankmap_lib_get_tni_list(int *tni_list, const int *flag);
  obtain the tni list
  (out) tni_list:  tni for [+x, -x, +y, -y,...] directions
  (in)  flag:      (not used) intended to specify the rank map
*/

#endif
