#! /usr/bin/python3

# Copyright 2021-2022. The MBI project. All rights reserved.
# This program is free software; you can redistribute it and/or modify it under the terms of the license (GNU GPL).

import os
import sys
import generator_utils as gen

template = """// @{generatedby}@
/* ///////////////////////// The MPI Bugs Initiative ////////////////////////

  Origin: MBI

  Description: @{shortdesc}@
    @{longdesc}@

   Version of MPI: Conforms to MPI 1.1, does not require MPI 2 implementation

BEGIN_MPI_FEATURES
  P2P!basic: Lacking
  P2P!nonblocking: Lacking
  P2P!persistent: Lacking
  COLL!basic: @{collfeature}@
  COLL!nonblocking: @{icollfeature}@
  COLL!persistent: Lacking
  COLL!tools: Lacking
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

#define buff_size 128

int main(int argc, char **argv) {
  int nprocs = -1;
  int rank = -1;
  int root = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("Hello from rank %d \\n", rank);

  if (nprocs < 2)
    printf("MBI ERROR: This test needs at least 2 processes to produce a bug.\\n");

  MPI_Comm newcom = MPI_COMM_WORLD;
  MPI_Datatype type = MPI_INT;
  MPI_Op op = MPI_SUM;

  int dbs = sizeof(int)*nprocs; /* Size of the dynamic buffers for alltoall and friends */
  @{init1}@
  @{init2}@

  if (@{change_cond}@) {
    @{operation1a}@ /* MBIERROR1 */
    @{fini1a}@
    @{operation2a}@
    @{fini2a}@
  } else {
    @{operation1b}@ /* MBIERROR2 */
    @{fini1b}@
    @{operation2b}@
    @{fini2b}@
  }

  @{free1}@
  @{free2}@

  MPI_Finalize();
  printf("Rank %d finished normally\\n", rank);
  return 0;
}
"""

for c1 in gen.coll + gen.icoll + gen.ibarrier:
    for c2 in gen.coll + gen.icoll + gen.ibarrier:
        patterns = {}
        patterns = {'c1': c1, 'c2': c2}
        patterns['generatedby'] = f'DO NOT EDIT: this file was generated by {os.path.basename(sys.argv[0])}. DO NOT EDIT.'
        patterns['collfeature'] = 'Yes' if c1 in gen.coll or c2 in gen.coll else 'Lacking'
        patterns['icollfeature'] = 'Yes' if c1 in gen.icoll + gen.ibarrier or c2 in gen.icoll + gen.ibarrier else 'Lacking'
        patterns['c1'] = c1
        patterns['c2'] = c2
        patterns['init1'] = gen.init[c1]("1")
        patterns['init2'] = gen.init[c2]("2")
        patterns['fini1a'] = gen.fini[c1]("1")
        patterns['fini2a'] = gen.fini[c2]("2")
        patterns['fini1b'] = gen.fini[c1]("1")
        patterns['fini2b'] = gen.fini[c2]("2")
        patterns['free1'] = gen.free[c1]("1")
        patterns['free2'] = gen.free[c2]("2")
        patterns['operation1a'] = gen.operation[c1]("1")
        patterns['operation1b'] = gen.operation[c1]("1")
        patterns['operation2a'] = gen.operation[c2]("2")
        patterns['operation2b'] = gen.operation[c2]("2")
        patterns['change_cond'] = 'rank % 2'
        shortdesc = ' collective ordering'

        if c1 == c2:
            # Generate the correct code using the same collective twice
            replace = patterns.copy()
            replace['shortdesc'] = 'Correct' + shortdesc
            replace['longdesc'] = f'All ranks call {c1} twice'
            replace['outcome'] = 'OK'
            replace['errormsg'] = ''
            replace['change_cond'] = 'rank < nprocs'
            replace['operation1b'] = ''
            replace['operation2b'] = ''
            replace['fini1b'] = ''
            replace['fini2b'] = ''
            gen.make_file(template, f'CallOrdering_{c1}_{c2}_ok.c', replace)
            # Generate the correct code using the collective once
            replace = patterns.copy()
            replace['shortdesc'] = 'Correct' + shortdesc
            replace['longdesc'] = f'All ranks call {c1} once'
            replace['outcome'] = 'OK'
            replace['errormsg'] = ''
            replace['init2'] = ''
            replace['change_cond'] = 'rank < nprocs'
            replace['operation2a'] = ''
            replace['operation1b'] = ''
            replace['operation2b'] = ''
            replace['fini2a'] = ''
            replace['fini1b'] = ''
            replace['fini2b'] = ''
            replace['free2'] = ''
            gen.make_file(template, f'CallOrdering_{c1}_ok.c', replace)
        else:
            # Generate the correct ordering with two different collectives
            replace = patterns.copy()
            replace['shortdesc'] = 'Correct' + shortdesc
            replace['longdesc'] = f'All ranks call {c1} and then {c2}'
            replace['outcome'] = 'OK'
            replace['errormsg'] = ''
            replace['change_cond'] = 'rank < nprocs'
            replace['operation1b'] = ''
            replace['operation2b'] = ''
            replace['fini1b'] = ''
            replace['fini2b'] = ''
            gen.make_file(template, f'CallOrdering_{c1}_{c2}_ok.c', replace)
            # Generate the incorrect ordering with two different collectives
            replace = patterns.copy()
            replace['shortdesc'] = 'Incorrect' + shortdesc
            replace['longdesc'] = f'Odd ranks call {c1} and then {c2} while even ranks call these collectives in the other order'
            replace['outcome'] = 'ERROR: CallMatching'
            replace['errormsg'] = 'Collective mistmatch. @{c1}@ at @{filename}@:@{line:MBIERROR1}@ is matched with @{c2}@ line @{filename}@:@{line:MBIERROR2}@.'
            replace['operation1b'] = gen.operation[c2]("2")  # Inversion
            replace['operation2b'] = gen.operation[c1]("1")
            replace['fini1a'] = gen.fini[c1]("1") # Inversion
            replace['fini2a'] = gen.fini[c2]("2")
            replace['fini1b'] = gen.fini[c2]("2") # Inversion
            replace['fini2b'] = gen.fini[c1]("1")
            replace['free1'] = gen.free[c2]("2")
            replace['free2'] = gen.free[c1]("1")

            gen.make_file(template, f'CallOrdering_{c1}_{c2}_nok.c', replace)

    # Generate the incorrect ordering with one collective
    replace = patterns.copy()
    replace['shortdesc'] = 'Incorrect' + shortdesc
    replace['longdesc'] = f'Odd ranks call {c1} while even ranks do not call any collective'
    replace['outcome'] = 'ERROR: CallMatching'
    replace['errormsg'] = 'Collective mistmatch. @{c1}@ at @{filename}@:@{line:MBIERROR1}@ is not matched.'
    replace['operation1b'] = ''  # Remove functions
    replace['operation2b'] = ''
    replace['operation2a'] = ''
    replace['init2'] = ''
    replace['fini1b'] = ''
    replace['fini2a'] = ''
    replace['fini2b'] = ''
    replace['free1'] = gen.free[c1]("1")
    replace['free2'] = ''
    gen.make_file(template, f'CallOrdering_{c1}_none_nok.c', replace)
    # Generate a correct ordering with a conditional not depending on ranks
    replace = patterns.copy()
    replace['shortdesc'] = 'Correct' + shortdesc
    replace['longdesc'] = f'All ranks call {c1}'
    replace['outcome'] = 'OK'
    replace['errormsg'] = ''
    replace['change_cond'] = 'rank < nprocs'
    replace['operation1b'] = '' # Remove functions
    replace['operation2b'] = ''
    replace['operation2a'] = ''
    replace['init2'] = ''
    replace['fini1b'] = ''
    replace['fini2a'] = ''
    replace['fini2b'] = ''
    replace['free1'] = gen.free[c1]("1")
    replace['free2'] = ''
    gen.make_file(template, f'CallOrdering_{c1}_none_ok.c', replace)
