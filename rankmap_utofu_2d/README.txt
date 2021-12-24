# Example of rankmap and TNI assignments

rankmap_lib_utofu.c : interface
get_tofu_coord_2d.c : 2-dim rankmap
get_tni_2d.c : TNI assignment

----------------------------
## make
----------------------------
```
make   # uses default TNI pattern
```
or
```
make TNI_FLAG=-D_USE_TNI_PATTERN3  # for pattern 2

```
or
```
make TNI_FLAG=-D_USE_TNI_PATTERN3  # for pattern 3

```
or
```
make TNI_FLAG=-D_USE_TNI_PATTERN4  # for pattern 4

```
See the TNI assignments section below for the details of each pattern.

----------------------------
## Rankmap
----------------------------

Requirement:
- the node allocation must be a hyper-rectangular shape
- each node has 4 MPI processes

If none of the shape bellow applies to the given node shape,
it sets the map_id < 0 and does not give a map.

- number of nodes < 48
   (B*A*Y*X) x (C*Z) node lattice with mesh geometry

- 6-dimensional node allocation matches
  (X,Y,Z,A,B,C)=(2,2,Z,2,3,2)
   it makes
   24 x 2Z node lattice with torus geometry
      = (axis made of X,Y,A and B) x (axis made of Z and C)

The inner node processes divided into 2x2 processes and the relation
between process coordinate and node coordinate is
- proc_x = 2*node_x + inner_x
- proc_y = 2*node_y + inner_y

----------------------------
## TNI assignment
----------------------------

There are 4 assignments defined in get_tni_2d.c and the default is pattern1.
You can specify with
  -D_USE_TNI_PATTERN2, -D_USE_TNI_PATTERN3 or -D_USE_TNI_PATTERN4
to switch to the other pattern.
The given TNI is used to receive the from the specified direction:
- dir=0:  from +x  (and default: send to -x)
- dir=1:  from -x  (and default: send to +x)
- dir=2:  from +y  (and default: send to -y)
- dir=3:  from -y  (and default: send to +y)

by using rdma_comlib_swap_vcq_for_sending() in rdma_utofu/rdma_utofu_comlib.c
one can swap the TNI for sending and make the same TNI to send to +x direction
and receive from +x direction.  See more details README file in utofu_rdma.

The patterns here are the same for all the nodes.  One can change the assignment
node by node in principle.
The pattern 1,2,4 calls some TNIs only 2 times while the other TNIs 4 times.
If the data sizes are the same for all directions, these pattern cause
a load imbalance between TNIs.  The pattern 3 calls each TNI at most 3 times
and thus has the smallest load imbalance if the data sizes are the same.

- PATTERN 1:  2 TNIs for inter-nodes, 4 TNIs for intra-node communication
  (0),...,(3) are MPI process in each node
  1,...,5 are TNI for receiving for each process
  If the data size are the same for all direction, the effective theoretical
  band width is 27.8GB/s, which is 16/24 of the theoretical peak.

       5        5
    4-(2)-1--0-(3)-4
       2        2          y
       |        |          ^
       3        3          |
    4-(0)-1--0-(1)-4       ---->x
       5        5


- PATTERN 2:   4 TNIs for inter-nodes, 2 TNIs for intra-node communication
  (0),...,(3) are MPI process in each node
  1,...,5 are TNI for receiving for each process
  If the data size are the same for all direction, the effective theoretical
  band width is 27.8GB/s, which is 16/24 of the theoretical peak.

       4        5
    3-(2)-0--0-(3)-3
       1        1          y
       |        |          ^
       1        1          |
    2-(0)-0--0-(1)-2       ---->x
       4        5


- PATTERN 3  round robin
  (0),...,(3) are MPI process in each node
  1,...,5 are TNI for receiving for each process
  If the data size are the same for all direction, the effective theoretical
  band width is 36.3GB/s, which is 16/18 of the theoretical peak.

       4        2
    3-(2)-2--1-(3)-0
       5        3          y
       |        |          ^
       2        0          |
    0-(0)-1--5-(1)-4       ---->x
       3        1

- PATTERN 4    4 TNIs for inter-nodes, 2 TNIs for intra-node communication
  (0),...,(3) are MPI process in each node
  1,...,5 are TNI for receiving for each process
  If the data size are the same for all direction, the effective theoretical
  band width is 27.2GB/s, which is 16/24 of the theoretical peak.

       2        2
    1-(2)-4--4-(3)-0
       5        5          y
       |        |          ^
       5        5          |
    1-(0)-4--4-(1)-0       ---->x
       3        3

----
