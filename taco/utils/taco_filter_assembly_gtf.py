import argparse
import logging
import sys

from taco.lib.gtf import GTF

def main():
    logging.basicConfig(level=logging.DEBUG)
    parser = argparse.ArgumentParser()
    parser.add_argument('--frac', type=float, default=0.0)
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    gtf_file = args.gtf_file
    frac_cutoff = args.frac

    all_t_ids = set()
    t_ids = set()
    for f in GTF.parse(open(gtf_file)):
        if f.feature == 'transcript':
            t_id = f.attrs['transcript_id']
            frac = float(f.attrs['rel_frac'])
            keep = (frac >= frac_cutoff)
            all_t_ids.add(t_id)
            if keep:
                t_ids.add(t_id)
                print str(f)
        elif f.feature == 'exon':
            t_id = f.attrs['transcript_id']
            assert t_id in all_t_ids
            if t_id in t_ids:
                print str(f)

if __name__ == '__main__':
    sys.exit(main())
