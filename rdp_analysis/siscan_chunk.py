#!/usr/bin/env python3
"""
Run OpenRDP's SiScan over one stride of the triplet set.

The `openrdp` CLI scans all C(37,3) = 7,770 triplets in a single serial process.
For SiScan that took >54 h on the smallest of our three alignments, and the 0.1
and 0.5 runs were killed by the kernel before finishing. This driver reproduces
exactly what `openrdp -m siscan` does to each triplet, but handles only the
triplets where `index % nchunks == chunk`, so the work can be spread over
several processes and each one keeps a small, bounded footprint.

Chunking does not change the numbers. Siscan.execute() calls random.seed() and
np.random.seed() with the same seed at the top of every triplet, so a triplet's
result depends only on the triplet itself, never on how many were scanned
before it.

Each chunk writes its raw (pre-merge) events to a pickle. merge_siscan_chunks.py
concatenates those and applies merge_breakpoints() once over the whole set, which
is what the single-process run would have done.
"""
import argparse
import os
import pickle
import sys
import time

from openrdp import Scanner
from openrdp.common import TripletGenerator, setup_upgma, upgma
from openrdp.siscan import Siscan


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('infile', help='alignment FASTA')
    p.add_argument('-c', '--cfg', required=True, help='OpenRDP config .ini')
    p.add_argument('-o', '--outfile', required=True, help='output pickle of raw events')
    p.add_argument('--chunk', type=int, required=True, help='this chunk index, 0-based')
    p.add_argument('--nchunks', type=int, required=True, help='total number of chunks')
    p.add_argument('-s', '--seed', type=int, default=3,
                   help='SiScan permutation seed; must match across chunks (default 3, '
                        'the same default the openrdp CLI uses)')
    p.add_argument('--progress-every', type=int, default=10,
                   help='report progress every N triplets handled by this chunk')
    p.add_argument('--start-index', type=int, default=0,
                   help='skip triplets with a global index below this (for resuming '
                        'or for isolating one stretch of the triplet set)')
    p.add_argument('--max-triplets', type=int, default=None,
                   help='stop after scanning this many triplets (diagnostics only; '
                        'produces a deliberately incomplete chunk)')
    args = p.parse_args()

    if not 0 <= args.chunk < args.nchunks:
        sys.exit(f'chunk {args.chunk} out of range for nchunks {args.nchunks}')

    # Scanner does the FASTA reading, duplicate-sequence removal and validation.
    # Reusing it keeps the alignment matrix byte-identical to a normal openrdp run.
    scanner = Scanner(cfg=args.cfg, methods=('siscan',), verbose=False, seed=args.seed)
    scanner._import_data(args.infile)

    settings = dict(scanner.config.items('Siscan'))
    siscan = Siscan(scanner.alignment, settings=settings, verbose=False)

    # SiScan needs the UPGMA tree to tell the major from the minor parent.
    tree, upgma_mat = setup_upgma(scanner.alignment, scanner.seq_names)
    tree, upgma_mat = upgma(tree, upgma_mat)

    triplets = TripletGenerator(scanner.alignment, scanner.seq_names)

    t0 = time.time()
    n_done = 0
    for idx, triplet in enumerate(triplets):
        if idx % args.nchunks != args.chunk or idx < args.start_index:
            continue
        if args.max_triplets is not None and n_done >= args.max_triplets:
            print(f'chunk {args.chunk}/{args.nchunks}: stopping at --max-triplets '
                  f'{args.max_triplets} (global index {idx}) -- INCOMPLETE CHUNK',
                  flush=True)
            break
        siscan.execute(triplet, tree=tree, random_seed=args.seed)
        n_done += 1
        if n_done % args.progress_every == 0:
            rate = n_done / max(time.time() - t0, 1e-9) * 3600
            print(f'chunk {args.chunk}/{args.nchunks}: {n_done} triplets '
                  f'(global index {idx}), {len(siscan.raw_results)} raw events, '
                  f'{rate:.1f} triplets/h', flush=True)

    # Deliberately unmerged -- the merge is global, see module docstring.
    with open(args.outfile, 'wb') as handle:
        pickle.dump({'chunk': args.chunk, 'nchunks': args.nchunks,
                     'n_triplets': n_done, 'seed': args.seed,
                     'infile': os.path.abspath(args.infile),
                     'raw_results': siscan.raw_results}, handle)

    print(f'chunk {args.chunk}/{args.nchunks} DONE: {n_done} triplets, '
          f'{len(siscan.raw_results)} raw events, '
          f'{(time.time() - t0) / 3600:.2f} h', flush=True)


if __name__ == '__main__':
    main()
