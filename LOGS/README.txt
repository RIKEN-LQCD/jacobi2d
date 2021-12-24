# serial
```
#PJM -L node=1
#PJM --mpi proc=1
export OMP_NUM_THREADS=48
mpiexec jacobi2d
```


# utofu double buffering, with rankmap
```
#PJM -L node=4x6x2:torus:strict
#PJM --mpi proc=192
export OMP_NUM_THREADS=12
mpiexec utofu_double_buf_jacobi2d 24 8
```

# mpi double buffering, with rankmap
```
#PJM -L node=4x6x2:torus:strict
#PJM --mpi proc=192
export OMP_NUM_THREADS=12
mpiexec mpi_double_buf_jacobi2d 24 8
```

# mpi single buffering, with rankmap
```
#PJM -L node=4x6x2:torus:strict
#PJM --mpi proc=192
export OMP_NUM_THREADS=12
mpiexec mpi_single_buf_jacobi2d 24 8
```
