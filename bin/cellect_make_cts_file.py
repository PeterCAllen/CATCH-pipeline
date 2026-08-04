#!/usr/bin/env python
"""Port of CELLECT scripts/make_cts_file_snake.py to an argv interface."""
import argparse
import os


def main(args):
    with open(args.annotations) as fh:
        annotations = [line.strip() for line in fh if line.strip()]

    if not annotations:
        raise ValueError(f"No annotations listed in {args.annotations}")

    with open(args.out, "w") as fh_out:
        for annotation in annotations:
            pre__annot = f"{args.run_prefix}__{annotation}"
            # The trailing dot is required: LDSC appends "<chr>.l2.ldscore.gz".
            pre__annot_path = os.path.join(args.ldscore_dir, pre__annot) + "."
            fh_out.write(f"{pre__annot}\t{pre__annot_path}\n")

    print(f"Wrote {args.out}: {len(annotations)} annotations")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotations', required=True)
    parser.add_argument('--run_prefix', required=True)
    parser.add_argument('--ldscore_dir', required=True)
    parser.add_argument('--out', required=True)
    main(parser.parse_args())
