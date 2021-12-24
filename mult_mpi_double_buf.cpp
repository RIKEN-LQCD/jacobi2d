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
  Hopping term with double buffer (persistent MPI version)
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2018 Dec.14, 2019 Feb. 6
*/

#include "mult_mpi_double_buf.h"
#include "comm_helper.h"

namespace MultMpiDoubleBuf {

namespace {
  // for treating fields
  inline int idx(int x, int y){return x+Nx*y;}
}

const char *Laplace2d::name="mpi double buffer";

// initialization
  void Laplace2d::init(const int buffer_extra, const int *neighbors){
  simple_timer timer;
  timer.reset();
  timer.start();

  // information on MPI ranks
  int rankid=communicator::rankid;

  // for receiving
  rnode[DirXp] = neighbors[0];
  rnode[DirXm] = neighbors[1];
  rnode[DirYp] = neighbors[2];
  rnode[DirYm] = neighbors[3];

  // for sending
  snode[DirXp] = rnode[DirXm];
  snode[DirXm] = rnode[DirXp];
  snode[DirYp] = rnode[DirYm];
  snode[DirYm] = rnode[DirYp];

  // buffer size
  m_buffer_extra=buffer_extra;
  int size_x=Ny*buffer_extra;
  int size_y=Nx*buffer_extra;
  size[DirXp]=size_x;
  size[DirXm]=size_x;
  size[DirYp]=size_y;
  size[DirYm]=size_y;
  if(rankid==0){
    printf("size_x size_y: %zd  %zd\n", size_x, size_y);
    fflush(0);
  }

  // allocate memory
  for(int dir=0; dir<4; dir++){
    rbuf_array[0][dir] = new double[size[dir]];
    rbuf_array[1][dir] = new double[size[dir]];
    sbuf[dir] = new double[size[dir]];
  }


  // alias of the recieve buffers
  for(int parity=0; parity<2; parity++){
    rbuf_xp_array[parity]=rbuf_array[parity][DirXp];
    rbuf_xm_array[parity]=rbuf_array[parity][DirXm];
    rbuf_yp_array[parity]=rbuf_array[parity][DirYp];
    rbuf_ym_array[parity]=rbuf_array[parity][DirYm];
  }

  // alias of the send buffers
  sbuf_xp=sbuf[DirXp];
  sbuf_xm=sbuf[DirXm];
  sbuf_yp=sbuf[DirYp];
  sbuf_ym=sbuf[DirYm];

  // initiate persistent communication
  for(int parity=0; parity<2; parity++){
    double **rbuf=rbuf_array[parity];  // send buffer is common for boht parity
    for(int dir=0; dir<4; dir++){
      int tag=dir+4*parity;
      MPI_Recv_init(rbuf[dir], sizeof(double)*size[dir], MPI_BYTE, rnode[dir], tag, MPI_COMM_WORLD, &rrequest[parity][dir]);
      MPI_Send_init(sbuf[dir], sizeof(double)*size[dir], MPI_BYTE, snode[dir], tag, MPI_COMM_WORLD, &srequest[parity][dir]);
    }
  }

  // set initial parity
  buf_parity=0;

  // initiate receiving
  MPI_Startall( 4, rrequest[buf_parity] );

  // reset the timer
  init_timer();

  // initializatoin, done
  timer.stop();
  if(rankid==0){
    double etime=timer.get_elapsed_msec();
    printf("init: %f [msec]\n", etime);
    fflush(0);
  }
}


// finalizaiton
void Laplace2d::finalize(){
  MPI_Barrier(MPI_COMM_WORLD);
  for(int dir=0; dir<4; dir++){
    MPI_Cancel(&rrequest[buf_parity][dir]);
  }

  // free the persistent communications
  for(int parity=0; parity<2; parity++){
    for(int dir=0; dir<4; dir++){
      MPI_Request_free(&rrequest[parity][dir]);
      MPI_Request_free(&srequest[parity][dir]);
      delete [] rbuf_array[parity][dir];
    }
  }
  for(int dir=0; dir<4; dir++){
    delete [] sbuf[dir];
  }
}

// hopping term
//  out(x,y) -= ( in(x+1,y) + in(x-1,y) + in(x,y+1) + in(x,y-1) )
void Laplace2d::hop(double *out, const double *in){
  timer_hop.start();

  // buffers for this parity
  double **rbuf=rbuf_array[buf_parity];
  double *rbuf_xp=rbuf[DirXp];
  double *rbuf_xm=rbuf[DirXm];
  double *rbuf_yp=rbuf[DirYp];
  double *rbuf_ym=rbuf[DirYm];

  // 1. initiate receicing
  //timer_irecv.start();
  //MPI_Startall( 4, rrequest[buf_parity] );
  //timer_irecv.stop();

  // 2. preapre send buffer
#pragma omp master
  {
    timer_pack.start();

  // pack: +x, -x
  for(int y=0; y<Ny; y++){
    sbuf_xp[y]= in[idx(0,y)];
    sbuf_xm[y]= in[idx(Nx-1,y)];
  }

  // pack: +y, -y
  for(int x=0; x<Nx; x++){
    sbuf_yp[x]= in[idx(x,0)];
    sbuf_ym[x]= in[idx(x,Ny-1)];
  }
  timer_pack.stop();
  }
#pragma omp barrier

  // 3. start sending
#pragma omp master
  {
  timer_communicate.start();
  MPI_Startall( 4, srequest[buf_parity] );
  timer_overlap.start();

  // 4. bulk computation
  timer_bulk.start();
  }

#pragma omp for collapse(2)
  for(int x=0; x<Nx; x++){
    for(int y=0; y<Ny; y++){
      int s=idx(x,y);

      // +x
      if(x!=(Nx-1)){
        out[s] -= in[idx(x+1,y)];
      }
      // -x
      if(x!=0){
        out[s] -= in[idx(x-1,y)];
      }

      // +y
      if(y!=(Ny-1)){
        out[s] -= in[idx(x,y+1)];
      }

      // -y
      if(y!=0){
        out[s] -= in[idx(x,y-1)];
      }
    } // y
  } // x
#pragma omp master
  {
  timer_bulk.stop();

  // 5a. intiate receiving
  timer_irecv.start();
  MPI_Startall( 4, rrequest[1-buf_parity] );
  timer_irecv.stop();


  // 5. wait for the receive buffer is ready
  timer_overlap.stop();
  MPI_Waitall(4,rrequest[buf_parity], MPI_STATUSES_IGNORE);
  timer_communicate.stop();


  // 6. boundary computation
  timer_boundary.start();

  for(int y=0; y<Ny; y++){
    // +x
    out[idx(Nx-1,y)] -= rbuf_xp[y];
    // -x
    out[idx(0,y)] -= rbuf_xm[y];
  }

  for(int x=0; x<Nx; x++){
    // +y
    out[idx(x,Ny-1)] -= rbuf_yp[x];
    // -y
    out[idx(x,0)] -= rbuf_ym[x];
  }
  timer_boundary.stop();

  // 7. guarantees the end of sending
  timer_swait.start();
  MPI_Waitall(4,srequest[buf_parity], MPI_STATUSES_IGNORE);
  timer_swait.stop();

  buf_parity=1-buf_parity;
  timer_hop.stop();
  } // omp master
#pragma omp barrier
}

// apply the matrix: (Laplacian + mass2 )
//   mass2 is defined in the confing.h
void Laplace2d::mult(double *out, const double *in){
  double d=4.0+mass2;
  for(int x=0; x<Nx; x++){
    for(int y=0; y<Ny; y++){
      out[idx(x,y)] = d*in[idx(x,y)];
    }
  }
  hop(out,in);
}


// reset the timer
void Laplace2d::init_timer(){
    timer_hop.reset();
    timer_communicate.reset();
    timer_overlap.reset();
    timer_boundary.reset();
    timer_pack.reset();
    timer_bulk.reset();
    timer_irecv.reset();
    timer_swait.reset();
}

// print the elapsed time
void Laplace2d::print_etime(){
    double etime_total=timer_hop.get_elapsed_msec();
    double etime_communicate=timer_communicate.get_elapsed_msec();
    double etime_overlap=timer_overlap.get_elapsed_msec();
    double etime_boundary=timer_boundary.get_elapsed_msec();
    double etime_pack=timer_pack.get_elapsed_msec();
    double etime_bulk=timer_bulk.get_elapsed_msec();
    double etime_irecv=timer_irecv.get_elapsed_msec();
    double etime_swait=timer_swait.get_elapsed_msec();

    printf("hop:      total  %f [msec]\n", etime_total);
    printf("           bulk  %f [msec]\n", etime_bulk);
    printf("       boundary  %f [msec]\n", etime_boundary);
    printf("           pack  %f [msec]\n", etime_pack);
    printf("    communicate  %f [msec]\n", etime_communicate);
    printf("      comm. overlap  %f [msec]\n", etime_overlap);
    printf("  comm. non-overlap  %f [msec]\n", etime_communicate-etime_overlap);
    printf(" initiate recv   %f [msec]\n", etime_irecv);
    printf("     send wait   %f [msec]\n", etime_swait);
    timer_hop.print_resolution();
    fflush(0);

}

  
} //MultMpiDoubleBuf
