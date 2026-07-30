#!/usr/bin/env python3
"""
Merge the per-chunk SiScan pickles from siscan_chunk.py into one CSV.

Concatenates the raw events from every chunk, then applies OpenRDP's
merge_breakpoints() once over the whole set, so the output matches what a
single-process `openrdp -m siscan` run would have written. Refuses to run if
any chunk is missing, since a silently short merge would look like a clean
"no events here" result.
"""
import argparse
import pickle
import sys

from openrdp.common import merge_breakpoints


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('pickles', nargs='+', help='per-chunk pickles')
    p.add_argument('-o', '--outfile', required=True, help='output CSV')
    p.add_argument('--nchunks', type=int, required=True,
                   help='expected number of chunks; all must be present')
    args = p.parse_args()

    raw, seen, total_triplets = [], set(), 0
    for path in args.pickles:
        with open(path, 'rb') as handle:
            d = pickle.load(handle)
        if d['nchunks'] != args.nchunks:
            sys.exit(f'{path}: nchunks {d["nchunks"]} != expected {args.nchunks}')
        if d['chunk'] in seen:
            sys.exit(f'{path}: duplicate chunk {d["chunk"]}')
        seen.add(d['chunk'])
        raw.extend(d['raw_results'])
        total_triplets += d['n_triplets']

    missing = sorted(set(range(args.nchunks)) - seen)
    if missing:
        sys.exit(f'ERROR: missing chunks {missing} -- refusing to write a partial merge')

    events = merge_breakpoints(raw)

    with open(args.outfile, 'w') as out:
        out.write('Method,Start,End,Recombinant,Parent1,Parent2,Pvalue\n')
        for e in events:
            out.write(f'Siscan,{e[2]},{e[3]},{e[0]},{e[1][0]},{e[1][1]},{e[4]}\n')

    print(f'{args.outfile}: {total_triplets} triplets, {len(raw)} raw events, '
          f'{len(events)} after merging')


if __name__ == '__main__':
    main()
