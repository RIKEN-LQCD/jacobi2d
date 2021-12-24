# double buffering communication with uTofu

-------------------------------
## data struct
-------------------------------

- rdma_comlib_data
  Contains communication buffers and everything.
  The most of the function take rdma_comlib_data* as the argument.
- rdma_comlib_communicator
  Contains Virtual Control Queue (VCQ) related information, rank ids to communicate with,
  and Tofu Network Interface (TNI) ids.
- rbuff_t
  Carries information of the receiving buffer.

-------------------------------
## initialize/finalize
-------------------------------

-  MPI must be initialized before calling rdma_comlib_init().
-  functions:
   - void rdma_comlib_init(void);
   - void rdma_comlib_finalize(void);

-------------------------------
## communication buffers
-------------------------------

- Needs to specify the following to generate
   - to which rank the data to be sent (dst_rank)
   - from which rank the data comes(rcv_rank)
   - TNI to use: both sending and receiving use the same TNI (and the same VCQ) by default.
   - size: data size in byte
   - The VCQ is allocated from the pool.

- functions
   - void rdma_comlib_new(rdma_comlib_data *id, const int *tni_id, const int *dst_rank, const int *rcv_rank, const size_t *size);
   - void rdma_comlib_delete(rdma_comlib_data *id);
   - void rdma_comlib_swap_vcq_for_sending(rdma_comlib_data *id1, rdma_comlib_data *id2);

- example:
```
    rdma_comlib_new(id1, tni1,  rank1,  rank2, size);
    rdma_comlib_new(id2, tni2,  rank2,  rank1, size);

//    rank2 --> tni1[vcq1] (ThisRank:id1)                       recv
//                         (ThisRank:id1) tni1[vcq1] --> rank1  send
//    rank2 <-- tni2[vcq2] (ThisRank:id2)                       send
//                         (ThisRank:id2) tni2[vcq2] <-- rank1  recv
```
   For each id, 1 VCQ is allocated, which uses the given TNI.  To communicate with
   rank1, e.g., both tni1 and tni2 are used.

rdma_comlib_swap_vcq_for_sending(id1, id2) swaps the VCQ (and thus TNI) for sending
between id1 and id2.
```
    rdma_comlib_new(id1, tni1,  rank1,  rank2, size);
    rdma_comlib_new(id2, tni2,  rank2,  rank1, size);
    rdma_comlib_swap_vcq_for_sending(id1, id2);
//    rank2 --> tni1[vcq1] (ThisRank:id1)                        recv
//                         (ThisRank:id1) tni2[vcq2] --> rank1   send
//    rank2 <-- tni1[vcq1] (ThisRank:id2)                        send
//                         (ThisRank:id2) tni2[vcq2] <-- rank1   recv
```
   To communicate with rank1, only tni2 is used.

-------------------------------
## communication
-------------------------------

- Starts/stops non-blocking communication.
- Internally switch the receive buffer for double buffering
- functions:
    - void rdma_comlib_irecv(rdma_comlib_data *id);
    - void rdma_comlib_isend(rdma_comlib_data *id);
    - void rdma_comlib_send_wait(rdma_comlib_data *id);
    - void rdma_comlib_recv_wait(rdma_comlib_data *id);


-------------------------------
## accesser
-------------------------------

- Pointer to send/receive buffers
- Do not modify the receive buffer, otherwise the cache injection mechanism does not work.
- functions
    - void *get_sbuff_ptr_rdma_comlib_data(rdma_comlib_data *id);
    - void *get_rbuff_ptr_rdma_comlib_data(rdma_comlib_data *id);

-----
