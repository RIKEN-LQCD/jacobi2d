//****************************************************************************************
//
//  Copyright (c) 2015-2020, Yoshifumi Nakamura <nakamura@riken.jp>
//  Copyright (c) 2015-2020, Yuta Mukai         <mukai.yuta@fujitsu.com>
//  Copyright (c) 2018-2020, Ken-Ichi Ishikawa  <ishikawa@theo.phys.sci.hirosima-u.ac.jp>
//  Copyright (c) 2019-2020, Issaku Kanamori    <kanamori-i@riken.jp>
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

#ifdef _USE_RANKMAP

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <assert.h>

//#define DEBUG_RANKMAP

#ifdef DEBUG_RANKMAP
#include <mpi-ext.h>
#include <utofu.h>
#endif

// rankmapt without using tofu coordiante:
//    lexical rankmap and round robin TNI

// tofu coordinate etc.: up to 4-dim
volatile int m_rankmap_lib_neighbor_ranks[8];
volatile int m_rankmap_lib_process_size[4];
volatile int m_rankmap_lib_process_coord[4];
volatile int m_rankmap_lib_dim;


int rankmap_lib_set_rankmap4d(){
  // not implemented
  return -1;
};

int rankmap_lib_set_rankmap4d_size(const int *proc_size ){
  // not implemented
  return -1;
};

int rankmap_lib_set_rankmap2d(){
  // not implemented
  return -1;
};

int rankmap_lib_set_rankmap2d_size(const int *proc_size){
  // assumes rank-map-bychip
  int rc;
  int myrank;
  int np;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  if(myrank==0){
    printf("setting 2-dim rank map\n");
    fflush(stdout);
  }
  if(np != proc_size[0] * proc_size[1]){
    if(myrank==0){
      fprintf(stderr, "np = %d != %d x %d\n", np, proc_size[0], proc_size[1]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  m_rankmap_lib_dim=2;
  // assumes 4 process per node
  const int node_size[2]={proc_size[0]/2, proc_size[1]/2};
  const int nodeid=myrank/4;

  const int node_x=nodeid % node_size[0];
  const int node_y=nodeid / node_size[0];

  const int rank_in_node = myrank % 4;
  const int rank_x_in_node = rank_in_node % 2;
  const int rank_y_in_node = rank_in_node / 2;

  const int Px = 2*node_x + rank_x_in_node;
  const int Py = 2*node_y + rank_y_in_node;

  m_rankmap_lib_process_size[0]=proc_size[0];
  m_rankmap_lib_process_size[1]=proc_size[1];
  m_rankmap_lib_process_coord[0]=Px;
  m_rankmap_lib_process_coord[1]=Py;

  // +x
  int node_xp = nodeid;
  int xp_in_rank=(rank_x_in_node+1)%2 + 2*rank_y_in_node;
  if(rank_x_in_node == 1){
    node_xp = (node_x+1)%node_size[0] + node_y*node_size[0];
  }
  m_rankmap_lib_neighbor_ranks[0]=4*node_xp + xp_in_rank;

  // -x
  int node_xm = nodeid;
  int xm_in_rank=(rank_x_in_node+1)%2 + 2*rank_y_in_node;
  if(rank_x_in_node == 0){
    node_xm = (node_x-1+node_size[0])%node_size[0] + node_y*node_size[0];
  }
  m_rankmap_lib_neighbor_ranks[1]=4*node_xm + xm_in_rank;


  // +y
  int node_yp = nodeid;
  int yp_in_rank=rank_x_in_node + 2*( (rank_y_in_node+1)%2 );
  if(rank_y_in_node == 1){
    node_yp = node_x + ( (node_y+1)%node_size[1] )*node_size[0];
  }
  m_rankmap_lib_neighbor_ranks[2]=4*node_yp + yp_in_rank;

  // -y
  int node_ym = nodeid;
  int ym_in_rank=rank_x_in_node + 2*( (rank_y_in_node+1)%2 );
  if(rank_y_in_node == 0){
    node_ym = node_x + ( (node_y-1+node_size[1])%node_size[1] )*node_size[0];
  }
  m_rankmap_lib_neighbor_ranks[3]=4*node_ym + ym_in_rank;


#ifdef DEBUG_RANKMAP
  {
    int dim=0;
    int rc = FJMPI_Topology_get_dimension(&dim);
    if (FJMPI_SUCCESS != rc) {
      fprintf(stderr, "rank %d: error from FJMPI_Topology_get_dimension\n", myrank);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (2 != dim) {
      fprintf(stderr, "rank %d: bad dim from FJMPI_Toplogy_get_dimension: dim=%d (must be 3)\n", myrank, dim);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    int shape_x, shape_y, shape_z;
    rc = FJMPI_Topology_get_shape(&shape_x, &shape_y, &shape_z);
    if (FJMPI_SUCCESS != rc) {
      fprintf(stderr, "rank %d: error from FJMPI_Topology_get_shape\n", myrank);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // get logical 3-dim coordiate
    int coords[2];
    rc = FJMPI_Topology_get_coords(MPI_COMM_WORLD, myrank, FJMPI_LOGICAL, 2, coords);
    if (FJMPI_SUCCESS != rc) {
      fprintf(stderr, "rank %d: error from FJMPI_Topology_get_coords\n", myrank);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    int this_x=coords[0];
    int this_y=coords[1];

    for(int i=0; i<np; i++){
      if(i==myrank){
        fprintf(stderr, "hoge: myrank=%d: Px=%d Py=%d; this_x=%d, this_y=%d; xp=%d, xm=%d, yp=%d, ym=%d\n",
                myrank, Px, Py, this_x, this_y,
                m_rankmap_lib_neighbor_ranks[0],
                m_rankmap_lib_neighbor_ranks[1],
                m_rankmap_lib_neighbor_ranks[2],
                m_rankmap_lib_neighbor_ranks[3]);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
#endif

  return 1;
};


void rankmap_lib_get_rankmap(int *myrank_coord, int *neighbors, int *process_size) {
  for(int i=0; i < m_rankmap_lib_dim; i++){
    myrank_coord[i]=m_rankmap_lib_process_coord[i];
    neighbors[2*i]=m_rankmap_lib_neighbor_ranks[2*i];
    neighbors[2*i+1]=m_rankmap_lib_neighbor_ranks[2*i+1];
    process_size[i]=m_rankmap_lib_process_size[i];
  }
}

//
// tni assignemnt: round robin
//
int rankmap_lib_get_tni_list(int *tni_list, const int *flag){

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if(m_rankmap_lib_dim!=2 && m_rankmap_lib_dim!=4){
    fprintf(stderr, "rank %d: cannot happen, m_rankmap_lib_dim=%d  (did you call set_rankmap* function?)\n",myrank, m_rankmap_lib_dim);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  if(myrank==0){
    fprintf(stdout, "rankmap_lib: using round robin TNI assignment");
  }
  for(int i=0; i<2*m_rankmap_lib_dim; i++){
    tni_list[i]=(2*m_rankmap_lib_dim*myrank+i) % 6;  // 6 is max of TIN id
  }
  return 0;

}

//

#endif
