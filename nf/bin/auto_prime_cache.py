#!/usr/bin/env python3

import sys
import subprocess


def main(args):
    """ Try to find the format of indexed file automatically.
        Works only with single input, not manifest """
    new_args = ["prime_cache"]
    is_input = False
    drop_arg = False
    for arg in args:
        if is_input:
            new_args.append("-i")
            new_args.append(arg)
            with open(arg, "rb") as f:
                start = f.read(9)
            if start == b'Seq-entry':
                format = "asn-seq-entry"
            elif start[0] == b'>':
                format = "fasta"
            else:
                format = "asnb-seq-entry"
            new_args.append("-ifmt")
            new_args.append(format)
            is_input = False
        elif drop_arg:
            drop_arg = False
        elif arg == "-i" or arg == "-input":
            is_input = True
        elif arg == "-ifmt":
            drop_arg = True
        else:
            new_args.append(arg)
    return subprocess.run(new_args).returncode


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))