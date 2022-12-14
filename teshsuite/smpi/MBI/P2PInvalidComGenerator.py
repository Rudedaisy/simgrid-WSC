#! /usr/bin/python3
import os
import sys
import generator_utils as gen

template = """// @{generatedby}@
/* ///////////////////////// The MPI Bugs Initiative ////////////////////////

  Origin: @{origin}@

  Description: @{shortdesc}@
    @{longdesc}@

   Version of MPI: Conforms to MPI 1.1, does not require MPI 2 implementation

BEGIN_MPI_FEATURES
  P2P!basic: @{p2pfeature}@
  P2P!nonblocking: @{ip2pfeature}@
  P2P!persistent: @{persfeature}@
  COLL!basic: Lacking
  COLL!nonblocking: Lacking
  COLL!persistent: Lacking
  COLL!tools: Yes
  RMA: Lacking
END_MPI_FEATURES

BEGIN_MBI_TESTS
  $ mpirun -np 2 ${EXE}
  | @{outcome}@
  | @{errormsg}@
END_MBI_TESTS
//////////////////////       End of MBI headers        /////////////////// */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char **argv) {
  int nprocs = -1;
  int rank = -1;
  int src=0, dest=1;
  int stag = 0, rtag = 0;
  int buff_size = 1;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("Hello from rank %d \\n", rank);

  MPI_Datatype type = MPI_INT;
  MPI_Comm newcom = MPI_COMM_WORLD;

  @{init1}@
  @{init2}@
  if (rank == 0) {
    @{change_com1}@
    @{operation1}@ /* MBIERROR1 */
    @{start1}@
    @{fini1}@
  }else if (rank == 1) {
    @{change_com2}@
    @{operation2}@ /* MBIERROR2 */
    @{start2}@
    @{fini2}@
  }
  @{free1}@
  @{free2}@

  if(newcom != MPI_COMM_NULL && newcom != MPI_COMM_WORLD)
    MPI_Comm_free(&newcom);

  MPI_Finalize();
  printf("Rank %d finished normally\\n", rank);
  return 0;
}
"""


for p1 in gen.send + gen.isend + gen.psend:
    for p2 in gen.recv + gen.irecv + gen.precv:
        patterns = {}
        patterns = {'p1': p1, 'p2': p2}
        patterns['origin'] = "MBI"
        patterns['generatedby'] = f'DO NOT EDIT: this file was generated by {os.path.basename(sys.argv[0])}. DO NOT EDIT.'
        patterns['p2pfeature'] = 'Yes' if p1 in gen.send or p2 in gen.recv  else 'Lacking'
        patterns['ip2pfeature'] = 'Yes' if p1 in gen.isend or p2 in gen.irecv  else 'Lacking'
        patterns['persfeature'] = 'Yes' if p1 in gen.psend or p2 in gen.precv  else 'Lacking'
        patterns['p1'] = p1
        patterns['p2'] = p2
        patterns['init1'] = gen.init[p1]("1")
        patterns['init2'] = gen.init[p2]("2")
        patterns['start1'] = gen.start[p1]("1")
        patterns['start2'] = gen.start[p2]("2")
        patterns['fini1'] = gen.fini[p1]("1")
        patterns['fini2'] = gen.fini[p2]("2")
        patterns['operation1'] = gen.operation[p1]("1") #send
        patterns['operation2'] = gen.operation[p2]("2") #recv
        patterns['free1'] = gen.free[p1]("1")
        patterns['free2'] = gen.free[p2]("2")
        patterns['change_com1'] = ""
        patterns['change_com2'] = ""

        replace = patterns.copy()
        replace['origin'] = "inspired from MPI-Corrbench"
        replace['shortdesc'] = 'Point to point @{p2}@ has an invalid communicator'
        replace['longdesc'] = 'MPI_COMM_NULL used in point to point @{p2}@'
        replace['outcome'] = 'ERROR: InvalidCommunicator'
        replace['errormsg'] = 'Invalid Communicator. @{p2}@ at @{filename}@:@{line:MBIERROR2}@ uses a null communicator.'
        replace['change_com2'] = 'newcom = MPI_COMM_NULL;'
        gen.make_file(template, f'InvalidParam_ComNull_{p2}_{p1}nok.c', replace)

        replace = patterns.copy()
        replace['shortdesc'] = 'Point to point @{p2}@ has an invalid communicator'
        replace['longdesc'] = 'MPI_COMM_NULL used in point to point @{p2}@'
        replace['outcome'] = 'ERROR: InvalidCommunicator'
        replace['errormsg'] = 'Invalid Communicator. @{p1}@ at @{filename}@:@{line:MBIERROR1}@ uses a null communicator.'
        replace['change_com1'] = 'newcom = MPI_COMM_NULL;'
        replace['change_com2'] = ""
        gen.make_file(template, f'InvalidParam_ComNull_{p1}_{p2}nok.c', replace)
