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
  solving 2dim Laplace eq. with Jacobi method: serial version
  Issaku Kanamori kanamori@hiroshima-u.ac.jp
  2018 Dec.14, 2019 Feb. 6
*/

#include "config.h"
#include "stdlib.h"
#include "stdio.h"
#include <math.h>

int Px=_Px;
int Py=_Py;

//int NNx=Nx*Px;
//int NNy=Ny*Py;
static const int NNx0=Nx*_Px;
static const int NNy0=Ny*_Py;

// site indexing
inline int idx(int x, int y){return x+NNx0*y;}

// squared norm
double norm2(const double *in){
  double sum2=0.0;
  for(int i=0; i<NNx0*NNy0; i++){
    sum2+=in[i]*in[i];
  }
  return sum2;
}

// set random field
void set_field(double *v){
  for(int i=0; i<NNx0*NNy0; i++){
    v[i]=( (double) rand())/RAND_MAX;
  }
  double n2=norm2(v);
  printf("set_field: n2 before normalizing=%24.15e\n", n2);
  double factor=sqrt(1.0/n2);
  for(int i=0; i<NNx0*NNy0; i++){
    v[i]*=factor;
  }  
}


// hopping term:
//  out(x,y) -= ( in(x+1,y) + in(x-1,y) + in(x,y+1) + in(x,y-1) )
void hop(double *out, const double *in){

  for(int x=0; x<NNx0; x++){
    for(int y=0; y<NNy0; y++){
      int s=idx(x,y);
      // diag
      //out[s] = 0.0;

      // +x
      int xp=x+1;
      if(x==NNx0-1){
	xp-=NNx0;
      }
      out[s] -= in[idx(xp,y)];

      // -x
      int xm=x-1;
      if(x==0){
	xm+=NNx0;
      }
      out[s] -= in[idx(xm,y)];

      // +y
      int yp=y+1;
      if(y==NNy0-1){
	yp-=NNy0;
      }
      out[s] -= in[idx(x,yp)];

      // -y
      int ym=y-1;
      if(y==0){
	ym+=NNy0;
      }
      out[s] -= in[idx(x,ym)];

    } // y
  } // x
}

// apply the matrix: (Laplacian + mass2 )
//   mass2 is defined in the confing.h
void mult(double *out, const double *in){
  double d=4.0+mass2;
  for(int x=0; x<NNx0; x++){
    for(int y=0; y<NNy0; y++){
      out[idx(x,y)] = d*in[idx(x,y)];
    }
  }
  hop(out,in);
}


// solve with jacobi method with fixed number of iteration
//   IterForSolver is defiend in config.h
void solve(double* solution, const double* source){
  int n_iter=IterForSolver;
  double v[NNx0*NNy0]; // work vector
  double factor=-1.0/(4.0+mass2);

  double n2_src=norm2(source);
  for(int xy=0; xy<NNx0*NNy0; xy++){
    solution[xy]=0.0;
  }
  for(int i=0; i<n_iter; i++){
    for(int xy=0; xy<NNx0*NNy0; xy++){
      v[xy]=-source[xy];
    }
    hop(v,solution);
    for(int xy=0; xy<NNx0*NNy0; xy++){
      solution[xy]=factor*v[xy];
    }

    //  mult(v,solution);
    //  for(int xy=0; xy<NNx0*NNy0; xy++){
    //    v[xy]-=source[xy];
    //  }  
    //  double n2_r=norm2(v);
    //  printf("iteration=%d, relative diff2: %g\n",i, n2_r/n2_src);
    if(i%norm2_freq == 0){
      double n2=norm2(solution);
      printf("iteration=%d, norm2(solution): %g\n",i, n2);
    }

  }

  mult(v,solution);
  for(int xy=0; xy<NNx0*NNy0; xy++){
    v[xy]-=source[xy];
  }  
  double n2_r=norm2(v);
  double n2  =norm2(solution);

  printf("iteration=%d, relative diff2: %g\n",n_iter, n2_r/n2_src);
  printf("norm2 of the solution=%24.15e\n", n2);

}


int main(int argc, char **argv){
  printf("Px, Py= %d, %d\n", Px, Py);

  double source[NNx0*NNy0];
  double solution[NNx0*NNy0];
  printf("NNx, NNy= %d, %d\n",NNx0,NNy0);

  set_field(source);

  solve(solution, source);

  double n2=norm2(solution);
  printf("done: norm2 of the solution=%24.15e\n", n2);
}

