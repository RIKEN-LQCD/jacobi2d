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
  Hopping term with single buffer (persistent MPI communication)
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2018 Dec.14, 2019 Feb. 6
*/


#include "mult_mpi_single_buf.h"
#include "comm_helper.h"

namespace MultMpiSingleBuf {


namespace {
  // for treating fields
  inline int idx(int x, int y){return x+Nx*y;}
}

const char *Laplace2d::name="mpi single buffer";

// initialization
  void Laplace2d::init(const int buffer_extra, const int *neighbors){
  simple_timer timer;
  timer.reset();
  timer.start();

    // information of MPI ranks
  int rankid=communicator::rankid;
  int pxp, pxm, pyp, pym;
  if(neighbors == 0){
    if(rankid==0){
      printf("%s: using the default rank map\n", name);
    }
    int npx=rankid % Px;
    int npy=rankid / Px;

    int np_xp=(npx+1) % Px;
    int np_xm=(npx-1+Px) % Px;
    int np_yp=(npy+1) % Py;
    int np_ym=(npy-1+Py) % Py;

    pxp = np_xp + Px*npy;
    pxm = np_xm + Px*npy;
    pyp = npx + Px*np_yp;
    pym = npx + Px*np_ym;
  } else {
    if(rankid==0){
      printf("%s: using the given rank map\n", name);
    }
    pxp = neighbors[0];
    pxm = neighbors[1];
    pyp = neighbors[2];
    pym = neighbors[3];
  }

  // for receiving
  rnode[DirXp] = pxp;
  rnode[DirXm] = pxm;
  rnode[DirYp] = pyp;
  rnode[DirYm] = pym;

  // for sending
  snode[DirXp] = rnode[DirXm];
  snode[DirXm] = rnode[DirXp];
  snode[DirYp] = rnode[DirYm];
  snode[DirYm] = rnode[DirYp];

  // buffer size
  m_buffer_extra=buffer_extra;
  size_t size_x=Ny*buffer_extra;
  size_t size_y=Nx*buffer_extra;
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
    rbuf[dir] = new double[size[dir]];
    sbuf[dir] = new double[size[dir]];
  }

  // alias of the buffers
  rbuf_xp=rbuf[DirXp];
  rbuf_xm=rbuf[DirXm];
  rbuf_yp=rbuf[DirYp];
  rbuf_ym=rbuf[DirYm];
  sbuf_xp=sbuf[DirXp];
  sbuf_xm=sbuf[DirXm];
  sbuf_yp=sbuf[DirYp];
  sbuf_ym=sbuf[DirYm];

  // initialize persistent communication
  for(int dir=0; dir<4; dir++){
    int tag=dir;
    MPI_Recv_init(rbuf[dir], sizeof(double)*size[dir], MPI_BYTE, rnode[dir], tag, MPI_COMM_WORLD, &rrequest[tag]);
    MPI_Send_init(sbuf[dir], sizeof(double)*size[dir], MPI_BYTE, snode[dir], tag, MPI_COMM_WORLD, &srequest[tag]);
  }

  // reset the timer
  init_timer();

  // initialization, done
  timer.stop();
  if(rankid==0){
    double etime=timer.get_elapsed_msec();
    printf("init: %f [msec]\n", etime);
  }

}

// finaization
void Laplace2d::finalize(){
  MPI_Barrier(MPI_COMM_WORLD);
  // free the persisetent communicaiton
  for(int dir=0; dir<4; dir++){
    MPI_Request_free(&rrequest[dir]);
    MPI_Request_free(&srequest[dir]);
    delete [] rbuf[dir];
    delete [] sbuf[dir];
  }
}

// hopping term
//  out(x,y) -= ( in(x+1,y) + in(x-1,y) + in(x,y+1) + in(x,y-1) )
void Laplace2d::hop(double *out, const double *in){
#pragma omp master
  {
    timer_hop.start();

  // 1. initiate receicing
  timer_irecv.start();
  MPI_Startall( 4, rrequest );
  timer_irecv.stop();

  // 2. preapre send buffer
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
  MPI_Startall( 4, srequest );
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

  // 5. wait for the reiceve buffer is ready
  timer_overlap.stop();
  MPI_Waitall(4,rrequest, MPI_STATUSES_IGNORE);
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
  MPI_Waitall(4,srequest, MPI_STATUSES_IGNORE);
  timer_swait.stop();
  // done!
  timer_hop.stop();
  }
#pragma omp barrier
}

// apply the matrix: (Laplacian + mass2 )
//   mass2 is defined in the confing.h
void Laplace2d::mult(double *out, const double *in){
  double d=4.0+mass2;
#pragma omp for collapse(2)
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

} // MultMpiSingleBuf
