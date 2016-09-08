import argparse

from sv_pipeline import *

parser = argparse.ArgumentParser()
parser.add_argument("type", help="Choose type of computation")
parser.add_argument("dir", help="Choose directory")
parser.add_argument("--out_pdf", help="Name of pdf to save files to")
parser.add_argument("--output_prefix", help="prefix of results.txt and figs file")
parser.add_argument("--L", help="Minimimum Matching Length")
parser.add_argument("--sw_window_size", help="Window size for smith waterman")
parser.add_argument("--sw_threshold", help="Threshold for smith waterman")

if __name__ == '__main__':
    args = parser.parse_args()
    DIR = args.dir
    out_pdf = args.out_pdf
    min_matching_length = int(args.L) if args.L else 600
    output_prefix = args.output_prefix if args.output_prefix else "output"
    sw_window_size = int(args.sw_window_size) if args.sw_window_size else 0
    sw_threshold = float(args.sw_threshold) if args.sw_threshold else 0.0

    print('DIR = %s' % DIR)

    if DIR == '':
        sys.exit('ERROR: DIRECTORY with SV calls required.')
    if args.type == "four":
        os.system("mkdir -p " + output_prefix + "_figs")
        four_graphs(DIR, min_matching_length, output_prefix, sw_window_size, sw_threshold)
    elif args.type == "sixteen":
        sixteen_graphs(DIR)
    elif args.type == "simple":
        simple(DIR)
    elif args.type == "classify":
        classify_reads(DIR)
    else:
        print("Unknown command")

    # combine pdfs if the system has pdftk
    if os.system("which pdftk") == 0:
        os.system('rm {}_figs/all.pdf'.format(output_prefix))
        os.system("pdftk {}_figs/*.pdf cat output {}".format(output_prefix, out_pdf))


