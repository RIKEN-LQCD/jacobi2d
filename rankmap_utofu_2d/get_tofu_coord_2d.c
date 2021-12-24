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

/*
   rank map for 2-dim system

 */
#ifdef _UTOFU_RDMA

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <utofu.h>


enum {DirX=0, DirY, DirZ, DirA, DirB, DirC};

// 2x2 in each node
const int proc_per_node=4;
const int pnx_in_node=2;
const int pny_in_node=2;


#define TOFU_MAX_IN_1AXIS 32
int check_tofu_volume(const uint8_t *my_coords, uint8_t *coords_org, uint8_t *coords_size, uint8_t *coords_min, uint8_t *coords_max, int *np){

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, np);

  MPI_Allreduce(my_coords, coords_min, 6, MPI_INTEGER1, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(my_coords, coords_max, 6, MPI_INTEGER1, MPI_MAX, MPI_COMM_WORLD);

  int my_occupation[6][TOFU_MAX_IN_1AXIS]={0};
  int occupied[6][TOFU_MAX_IN_1AXIS]={0};
  for(int i=0; i<6; i++){
    my_occupation[i][my_coords[i]]=1;
  }
  MPI_Allreduce(my_occupation, occupied, 6*TOFU_MAX_IN_1AXIS,
                MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // count the number of nodes in each direction
  for(int i=0; i<6; i++) {
    int size=0;
    for(int n=0; n<TOFU_MAX_IN_1AXIS; n++){
      if(occupied[i][n]>0) { size++; }
    }
    coords_size[i] = size;
  }
  // look for the origin
  for(int i=0; i<6; i++){
    int org=0;
    int prev=0;
    for(int n=TOFU_MAX_IN_1AXIS; n>0; n--){
      if(occupied[i][n-1] < prev) {
        org=n;
        break;
      }
      prev=occupied[i][n-1];
    }
    coords_org[i]=org;
  }

  //
  // periodic condition:
  //   if coord[i] > coord_max[i]
  //     coord[i] = coords_min[i] + coord[i] % coord_max[i]
  //
  //   if coord[i] < coord_org[i]
  //     coord[i] = coord_max[i] + coord[i] - coord_min[i]
  //
  //  if tofu in i-th dir is torus and
  //     sytem         7 8 9 10 11 12 13 14 15
  //     *:using/-not  * * - -  -  -  -  *  *
  //  coords_min[i]=7
  //  coords_max[i]=15
  //  coords_org[i]=14
  //  coords_size[i]=4
  //    (N.B. in the above case, naive min/max of tofu coordinate
  //     becomes 7 and 15, respectively.)

  // check the FJMPI interface as well
  //int rel_coords[6];
  //  FJMPI_Topology_get_coords(MPI_COMM_WORLD, myrank, FJMPI_TOFU_REL, 6, rel_coords);
  //  printf("rank %d: tofu coords sys = %d %d %d %d %d %d; rel = %d %d %d %d %d %d\n",
  //         myrank,
  //         my_coords[0],my_coords[1],my_coords[2],my_coords[3],my_coords[4],my_coords[5],
  //         rel_coords[0],rel_coords[1],rel_coords[2],rel_coords[3],rel_coords[4],rel_coords[5]);

  int tofu_vol=1;
  for(int i=0; i<6; i++){
    tofu_vol*=coords_size[i];
  }
  if(tofu_vol*4 != *np){
    if(myrank==0){
      fprintf(stderr, "warning: allocated size is not (hyper-)rectangluer:\n");
      fprintf(stderr, "       np = %d, tofu_vol * 4 = %d\n", *np, tofu_vol*4);
      fprintf(stderr, "       dir: min  max  size origin\n");
      for(int i=0; i<6; i++){
        fprintf(stderr, "        %d: %3d  %3d  %3d  %3d\n", i, coords_min[i], coords_max[i], coords_size[i], coords_org[i]);
      }
    }
    return -1;
    //MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    //    return RANKMAP_TOPOLOGY_Y;

    // or find a largest hyper rectangluar which fit in the given nodes
    //  np = ...
  }

  if(coords_size[DirA] !=2 || coords_size[DirB] !=3 || coords_size[DirC] !=2 ){
    if(myrank==0){
      fprintf(stderr, "bad mpi size: must be [X,Y,Z,A,B,C]=[*,*,*,2,3,2]  (but [%d,%d,%d,%d,%d,%d])\n",
              coords_size[DirX], coords_size[DirY], coords_size[DirZ], coords_size[DirA], coords_size[DirB], coords_size[DirC]);
    }
    return -1;
  }

  if(myrank==0){
    fprintf(stderr, "       np = %d, tofu_vol * 4 = %d\n",*np, tofu_vol*4);
    fprintf(stderr, "       dir: min  max  size origin\n");
    for(int i=0; i<6; i++){
      fprintf(stderr, "        %d: %3d  %3d  %3d  %3d\n", i, coords_min[i], coords_max[i], coords_size[i], coords_org[i]);
    }
  }

  return 0;
}



int get_tofu_coord_2d_eval1(const int myrank, const uint8_t *my_coords,
                            const uint8_t *coords_org, const uint8_t *coords_size, const uint8_t *coords_min, const uint8_t *coords_max,
                            int *rank_coord, int *rank_size,
                            uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                            uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node){
  if(myrank==0){
    fprintf(stderr, "rankmap for eval1\n");
    fprintf(stderr, "  coords_org: %d %d %d %d %d %d\n", coords_org[0], coords_org[1], coords_org[2], coords_org[3], coords_org[4], coords_org[5]);
    fprintf(stderr, "  coords_size: %d %d %d %d %d %d\n", coords_size[0], coords_size[1], coords_size[2], coords_size[3], coords_size[4], coords_size[5]);
    fprintf(stderr, "  coords_min: %d %d %d %d %d %d\n", coords_min[0], coords_min[1], coords_min[2], coords_min[3], coords_min[4], coords_min[5]);
    fprintf(stderr, "  coords_max: %d %d %d %d %d %d\n", coords_max[0], coords_max[1], coords_max[2], coords_max[3], coords_max[4], coords_max[5]);
  }


  // edit here to change the rank map
  //   this is a lexical mapping, may not be optimal

  // (B,A,Y,X) --> x
  int xdim=4;
  int dir_x[4]={DirB, DirA, DirY, DirX};

  // (C,Z) --> y
  int ydim=2;
  int dir_y[2]={DirC, DirZ};

  /*
  // (B,A,Y,X,Z) --> x
  int xdim=5;
  int dir_x[5]={DirB, DirA, DirY, DirX, DirZ};
  // (C) --> y
  int ydim=1;
  int dir_y[1]={DirC};
  */


  // process coordinate in the node
  int tmp=myrank % proc_per_node;
  int px=tmp%pnx_in_node;
  int py=tmp/pnx_in_node;

  int err=0;

  int rcoords[6]; // relative tofu coordinates
  int size[6];    // tofu size
  for(int i=0; i<6; i++){
    if(my_coords[i] < coords_org[i]){
      rcoords[i] = (coords_max[i] + my_coords[i] - coords_min[i]) - coords_org[i];
    } else {
      rcoords[i] = my_coords[i] - coords_org[i];
    }
    size[i]=coords_size[i];
  }

  // size in each direction
  int size_x=pnx_in_node;
  int size_y=pny_in_node;
  for(int i=0; i<xdim; i++){
    size_x*=size[dir_x[i]];
  }
  for(int i=0; i<ydim; i++){
    size_y*=size[dir_y[i]];
  }
  rank_size[0] = size_x;
  rank_size[1] = size_y;

  // logical coordiante
  int rank_x=rcoords[dir_x[0]];
  int rank_y=rcoords[dir_y[0]];
  for(int i=1; i<xdim; i++){
    rank_x+=rcoords[dir_x[i]]*size[dir_x[i-1]];
  }
  for(int i=1; i<ydim; i++){
    rank_y+=rcoords[dir_y[i]]*size[dir_y[i-1]];
  }
  rank_coord[0]=px+pnx_in_node*rank_x;
  rank_coord[1]=py+pny_in_node*rank_y;


  for(int i=0; i<6; i++){
    positive_neighbor_coords[0][i]=rcoords[i];
    positive_neighbor_coords[1][i]=rcoords[i];
    negative_neighbor_coords[0][i]=rcoords[i];
    negative_neighbor_coords[1][i]=rcoords[i];
  }

  //
  // x-direction
  //
  if(px==0){
    // px+1 =  0+1    = 1; positive_neighbor_coords does not change
    // px-1 = (0-1)%2 = 1
    pos_rank_in_node[0] =1 + py*pnx_in_node;
    neg_rank_in_node[0] =1 + py*pnx_in_node;

    // logical node coorditate -=1
    uint8_t *ncoords=negative_neighbor_coords[0];
    const int *Dir=dir_x;
    const int dim=xdim;
    int carry_over= -1;
    for(int i=0; i<dim; i++){
      ncoords[Dir[i]]+= (carry_over + size[Dir[i]]);
      carry_over=(ncoords[Dir[i]]/size[Dir[i]])-1;
      ncoords[Dir[i]] = ncoords[Dir[i]] % size[Dir[i]];
    }
  } else {
    // px+1 = (1+1) % 2 = 0
    // px-1 = (1-1) % 2 = 0; negative_neighbor_coords does not change
    pos_rank_in_node[0] =0 + py*pnx_in_node;
    neg_rank_in_node[0] =0 + py*pnx_in_node;

    // logical node coordinate +=1
    uint8_t *ncoords=positive_neighbor_coords[0];
    const int *Dir=dir_x;
    const int dim=xdim;
    int carry_over=1;
    for(int i=0; i<dim; i++){
      ncoords[Dir[i]]+=carry_over;
      carry_over=ncoords[Dir[i]]/size[Dir[i]];
      ncoords[Dir[i]] = ncoords[Dir[i]] % size[Dir[i]];
    }
  }

  //
  // y-direction
  //
  if(py==0){
    // py+1 =  0+1    = 1; positive_neighbor_coords does not change
    // py-1 = (0-1)%2 = 1
    pos_rank_in_node[1] =px + 1*pnx_in_node;
    neg_rank_in_node[1] =px + 1*pnx_in_node;

    // lgocal node coordinate -=1
    uint8_t *ncoords=negative_neighbor_coords[1];
    const int *Dir=dir_y;
    const int dim=ydim;
    int carry_over= -1;
    for(int i=0; i<dim; i++){
      ncoords[Dir[i]]+= (carry_over + size[Dir[i]]);
      carry_over=(ncoords[Dir[i]]/size[Dir[i]])-1;
      ncoords[Dir[i]] = ncoords[Dir[i]] % size[Dir[i]];
    }
  } else {
    // py+1 = (1+1) % 2 = 0
    // py-1 = (1-1) % 2 = 0; negative_neighbor_coords does not change
    pos_rank_in_node[1] =px + 0*pnx_in_node;
    neg_rank_in_node[1] =px + 0*pnx_in_node;

    // logical node coordinate +=1
    uint8_t *ncoords=positive_neighbor_coords[1];
    const int *Dir=dir_y;
    const int dim=ydim;
    int carry_over=1;
    for(int i=0; i<dim; i++){
      ncoords[Dir[i]]+=carry_over;
      carry_over=ncoords[Dir[i]]/size[Dir[i]];
      ncoords[Dir[i]] = ncoords[Dir[i]] % size[Dir[i]];
    }
  }


  for(int dir=0; dir <2; dir++){
    for(int i=0; i<6; i++){
      positive_neighbor_coords[dir][i]+=coords_org[i];
      if(positive_neighbor_coords[dir][i]>coords_max[i]){
        positive_neighbor_coords[dir][i]-=coords_max[i];
        positive_neighbor_coords[dir][i]+=coords_min[i];
      }
      negative_neighbor_coords[dir][i]+=coords_org[i];
      if(negative_neighbor_coords[dir][i]>coords_max[i]){
        negative_neighbor_coords[dir][i]-=coords_max[i];
        negative_neighbor_coords[dir][i]+=coords_min[i];
      }
    }
  }

  //  for(int dir=0; dir <2; dir++){
  //    uint8_t *p=positive_neighbor_coords[dir];
  //    uint8_t *n=negative_neighbor_coords[dir];
  //    fprintf(stderr, "  rank=%d, dir=%d: positive_ncoords = %d %d %d %d %d, nagative_ncoords = %d %d %d %d %d %d\n",
  //            myrank, dir, p[0],p[1],p[2],p[3],p[4],p[5], n[0],n[1],n[2],n[3],n[4],n[5]);
  //  }

  int map_id=1;
  if(err<0) {
    map_id=err;
  }
  return map_id;
}


// TX=TY=2
int get_tofu_coord_2d_eval2(const int myrank, const uint8_t *my_coords,
                            const uint8_t *coords_org, const uint8_t *coords_size, const uint8_t *coords_min, const uint8_t *coords_max,
                            int *rank_coord, int *rank_size,
                            uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                            uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node){
  if(myrank==0){
    fprintf(stderr, "rankmap for eval2 (v2)\n");
    fprintf(stderr, "  coords_org: %d %d %d %d %d %d\n", coords_org[0], coords_org[1], coords_org[2], coords_org[3], coords_org[4], coords_org[5]);
    fprintf(stderr, "  coords_size: %d %d %d %d %d %d\n", coords_size[0], coords_size[1], coords_size[2], coords_size[3], coords_size[4], coords_size[5]);
    fprintf(stderr, "  coords_min: %d %d %d %d %d %d\n", coords_min[0], coords_min[1], coords_min[2], coords_min[3], coords_min[4], coords_min[5]);
    fprintf(stderr, "  coords_max: %d %d %d %d %d %d\n", coords_max[0], coords_max[1], coords_max[2], coords_max[3], coords_max[4], coords_max[5]);
  }

  int rcoords[6]; // relative tofu coordinates
  int size[6];    // tofu size
  for(int i=0; i<6; i++){
    if(my_coords[i] < coords_org[i]){
      rcoords[i] = (coords_max[i] + my_coords[i] - coords_min[i]) - coords_org[i];
    } else {
      rcoords[i] = my_coords[i] - coords_org[i];
    }
    size[i]=coords_size[i];
  }
  if(myrank==0){
    fprintf(stderr,"tofu size: [A,B,C,X,Y,Z]=[%d,%d,%d,%d,%d,%d]\n",
            size[DirA],size[DirB],size[DirC],size[DirX],size[DirY],size[DirZ]);
  }

  int TA[24]={0,1,1,0,0,1, 1,0,0,1,1,0, 0,1,1,0,0,1, 1,0,0,1,1,0};
  int TB[24]={0,0,1,1,2,2, 2,2,0,0,1,1, 1,1,0,0,2,2, 2,2,1,1,0,0};
  int TX[24]={0,0,0,0,0,0, 1,1,1,1,1,1, 1,1,1,1,1,1, 0,0,0,0,0,0};
  int TY[24]={0,0,0,0,0,0, 0,0,0,0,0,0, 1,1,1,1,1,1, 1,1,1,1,1,1};
  int size_x = 24;
  int size_y = 2*size[DirZ];


  int TC[48];
  int TZ[48];
  for(int i=0;i<48; i++){
    TC[i]=-1;
    TZ[i]=-1;
  }
    int n=0;
    TC[n]=0;
    TZ[n]=0;
    n++;
    for(int z=0; z<size[DirZ]; z++){
      TC[n]=1;
      TZ[n]=z;
      n++;
    }
    for(int z=size[DirZ]-1;  z>0; z--){
      TC[n]=0;
      TZ[n]=z;
      n++;
    }

  // process coordinate in the node
  int tmp=myrank % proc_per_node;
  int px=tmp%pnx_in_node;
  int py=tmp/pnx_in_node;

  int err=0;


  // size in each direction
  if(size[DirA] !=2 || size[DirB] !=3 || size[DirC] !=2 || size[DirX] !=2 || size[DirY] !=2){
    if(myrank==0){
      fprintf(stderr, "bad mpi size: must be [A,B,C,X,Y,Z]=[2,3,2,2,2,*]\n");
    }
    return -9999;
  }

  if(myrank==0){
    //    printf("rankmap for Fugaku evaluation environment 2\n");
    //    printf("tofu size (original): [A,B,C,X,Y,Z]=[%d,%d,%d,%d,%d,%d]\n",
    //            size[DirA_],size[DirB_],size[DirC_],size[DirX_],size[DirY_],size[DirZ_]);
    //    printf("  rotate: X,Y,Z --> %c,%c,%c\n", tofu_char[DirX], tofu_char[DirY], tofu_char[DirZ]);
    printf("tofu size: [A,B,C,X,Y,Z]=[%d,%d,%d,%d,%d,%d]\n",
            size[DirA],size[DirB],size[DirC],size[DirX],size[DirY],size[DirZ]);
    printf("----- node coordinates -----\n");
    printf("map for x:  TA  = ");
    for(int i=0; i<size_x; i++){printf("%3d",TA[i]); }
    printf("\n");
    printf("            TB  = ");
    for(int i=0; i<size_x; i++){printf("%3d",TB[i]); }
    printf("\n");
    printf("            TX  = ");
    for(int i=0; i<size_x; i++){printf("%3d",TX[i]); }
    printf("\n");
    printf("            TY  = ");
    for(int i=0; i<size_x; i++){printf("%3d",TY[i]); }
    printf("\n");

    printf("map for y:  TC  = ");
    for(int i=0; i<size_y; i++){printf("%3d",TC[i]); }
    printf("\n");
    printf("            TZ  = ");
    for(int i=0; i<size_y; i++){printf("%3d",TZ[i]); }
    printf("\n");

    printf("inner node size: 2 2\n");
    printf("----------------------------\n");
    fflush(0);
  }


  rank_size[0] = pnx_in_node*size_x;
  rank_size[1] = pny_in_node*size_y;

  int node_x=-1;
  for(int i=0; i<size_x; i++){
    if(rcoords[DirX] == TX[i]
       && rcoords[DirY] == TY[i]
       && rcoords[DirA] == TA[i]
       && rcoords[DirB] == TB[i]){
      node_x=i;
      break;
    }
  }
  if(node_x<0){
    fprintf(stderr, "cannot happen: no node_x is found\n");
    return -999;
  }

  int node_y=-1;
  for(int i=0; i<size_y; i++){
    if(rcoords[DirZ] == TZ[i]
       && rcoords[DirC] == TC[i] ){
      node_y=i;
      break;
    }
  }
  if(node_y<0){
    fprintf(stderr, "cannot happen: no node_y is found\n");
    return -999;
  }


  // logical coordiante
  rank_coord[0]=px+pnx_in_node*node_x;
  rank_coord[1]=py+pny_in_node*node_y;


  // neighbors
  for(int i=0; i<6; i++){
    positive_neighbor_coords[0][i]=rcoords[i];
    positive_neighbor_coords[1][i]=rcoords[i];
    negative_neighbor_coords[0][i]=rcoords[i];
    negative_neighbor_coords[1][i]=rcoords[i];
  }

  int dir;
  //
  // x-direction : has inner node ranks
  //
  dir=0;
  int xf=(rank_coord[dir]+1) % rank_size[dir];
  int xb=(rank_coord[dir]-1+rank_size[dir]) % rank_size[dir];
  int pxf=xf % pnx_in_node;
  int pxb=xb % pnx_in_node;
  xf/=pnx_in_node;
  xb/=pnx_in_node;
  pos_rank_in_node[dir] = pxf + py*pnx_in_node;
  neg_rank_in_node[dir] = pxb + py*pnx_in_node;

  positive_neighbor_coords[dir][DirA]=TA[xf];
  positive_neighbor_coords[dir][DirB]=TB[xf];
  positive_neighbor_coords[dir][DirX]=TX[xf];
  positive_neighbor_coords[dir][DirY]=TY[xf];
  negative_neighbor_coords[dir][DirA]=TA[xb];
  negative_neighbor_coords[dir][DirB]=TB[xb];
  negative_neighbor_coords[dir][DirX]=TX[xb];
  negative_neighbor_coords[dir][DirY]=TY[xb];

  //
  // y-direction : has inner node ranks
  //
  dir=1;
  int yf=(rank_coord[dir]+1) % rank_size[dir];
  int yb=(rank_coord[dir]-1+rank_size[dir]) % rank_size[dir];
  int pyf=yf % pny_in_node;
  int pyb=yb % pny_in_node;
  yf/=pny_in_node;
  yb/=pny_in_node;
  pos_rank_in_node[dir] = px + pyf*pnx_in_node;
  neg_rank_in_node[dir] = px + pyb*pnx_in_node;

  positive_neighbor_coords[dir][DirC]=TC[yf];
  positive_neighbor_coords[dir][DirZ]=TZ[yf];
  negative_neighbor_coords[dir][DirC]=TC[yb];
  negative_neighbor_coords[dir][DirZ]=TZ[yb];

  for(int dir=0; dir <2; dir++){
    for(int i=0; i<6; i++){
      positive_neighbor_coords[dir][i]+=coords_org[i];
      if(positive_neighbor_coords[dir][i]>coords_max[i]){
        positive_neighbor_coords[dir][i]-=coords_max[i];
        positive_neighbor_coords[dir][i]+=coords_min[i];
      }
      negative_neighbor_coords[dir][i]+=coords_org[i];
      if(negative_neighbor_coords[dir][i]>coords_max[i]){
        negative_neighbor_coords[dir][i]-=coords_max[i];
        negative_neighbor_coords[dir][i]+=coords_min[i];
      }
    }
  }
  //  for(int dir=0; dir <2; dir++){
  //    uint8_t *p=positive_neighbor_coords[dir];
  //    uint8_t *n=negative_neighbor_coords[dir];
  //    const uint8_t *m=my_coords;
  //    fprintf(stderr, "  rank=%d [%d %d %d %d %d %d], dir=%d: positive_ncoords = %d %d %d %d %d %d, nagative_ncoords = %d %d %d %d %d %d\n",
  //            myrank, m[0], m[1], m[2], m[3], m[4], m[5],
  //            dir, p[0],p[1],p[2],p[3],p[4],p[5], n[0],n[1],n[2],n[3],n[4],n[5]);
  //  }

  int map_id=2;
  if(err<0) {
    map_id=err;
  }
  return map_id;
}


int get_tofu_coord_2d(const int myrank, const uint8_t *my_coords,
                      const uint8_t *coords_org, const uint8_t *coords_size, const uint8_t *coords_min, const uint8_t *coords_max,
                      int *rank_coord, int *rank_size,
                      uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                      uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node){
  // process coordinate in the node

  int tofu_vol=1;
  for(int i=0; i<6; i++){
    tofu_vol*=(coords_max[i]-coords_min[i]+1);
  }
  if(tofu_vol >= 48){ // evaluation 2
    return get_tofu_coord_2d_eval2(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                   rank_coord, rank_size,
                                   positive_neighbor_coords, pos_rank_in_node,
                                   negative_neighbor_coords, neg_rank_in_node);

  } else {
    return get_tofu_coord_2d_eval1(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                                   rank_coord, rank_size,
                                   positive_neighbor_coords, pos_rank_in_node,
                                   negative_neighbor_coords, neg_rank_in_node);
  }
}


int get_tofu_coord_and_tni(uint8_t *my_coords, int *rank_coord, int *rank_size,
                           uint8_t (*positive_neighbor_coords)[6], int *pos_rank_in_node,
                           uint8_t (*negative_neighbor_coords)[6], int *neg_rank_in_node,
                           int *tni_list){

  int myrank, np_available;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  const int dim=2;
  // get tofu coordinate
  int rc=utofu_query_my_coords(my_coords);
  if(UTOFU_SUCCESS != rc){
    fprintf(stderr, "rank %d: Failed at utofu_querry_my_coords()! rc=%d\n", myrank, rc);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  // check if the assigned network is a hyper-rectangular
  //   it also obtain the assinged tofu coordiate
  //   N.B:  if the XXX axis use the torus nature,
  //         coords_min[XXX] + coords_size[XXX] != coords_max[XXX]

  uint8_t coords_org[6], coords_size[6];
  uint8_t coords_min[6], coords_max[6];
  int pre_mapid=check_tofu_volume(my_coords, coords_org, coords_size,
                                  coords_min, coords_max, &np_available);
  if(pre_mapid<0){
    if(myrank==0){
      fprintf(stderr, "rank %d: Failed at check_tofu_volume()! err=%d\n", myrank, pre_mapid);
    }
    //    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    return pre_mapid;
  }

  int flag;
  flag=get_tofu_coord_2d(myrank, my_coords, coords_org, coords_size, coords_min, coords_max,
                         rank_coord, rank_size,
                         positive_neighbor_coords, pos_rank_in_node,
                         negative_neighbor_coords, neg_rank_in_node);
  int rc2;
  rc2 = get_tni_list(tni_list, myrank, rank_size, dim, flag);
  if(myrank==0){
    printf("tni_list = %d %d %d %d: dim=%d flag=%d\n", tni_list[0], tni_list[1], tni_list[2], tni_list[3], dim, flag);
  }
  return flag;
}
#endif
