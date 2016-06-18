import argparse

from sv_pipeline import *

parser = argparse.ArgumentParser()
parser.add_argument("type", help="Choose type of computation")
parser.add_argument("dir", help="Choose directory")
parser.add_argument("out_pdf", help="Name of pdf to save files to")
parser.add_argument("--L", help="Minimimum Matching Length")

if __name__ == '__main__':
    args = parser.parse_args()
    DIR = args.dir
    out_pdf = args.out_pdf
    min_matching_length = args.L if args.L else 600

    print('DIR = %s' % DIR)

    if DIR == '':
        sys.exit('ERROR: DIRECTORY with SV calls required.')
    if args.type == "four":

        four_graphs(DIR, min_matching_length)
    elif args.type == "sixteen":
        sixteen_graphs(DIR)
    elif args.type == "simple":
        simple(DIR)
    elif args.type == "classify":
        classify_reads(DIR)
    else:
        print("Unknown command")
    os.system("pdftk figs/*.pdf cat output {}".format(out_pdf))


