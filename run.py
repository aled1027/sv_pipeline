import argparse

from sv_pipeline import *

parser = argparse.ArgumentParser()
parser.add_argument("type", help="Choose type of computation")
parser.add_argument("dir", help="Choose directory")

if __name__ == '__main__':
    args = parser.parse_args()
    DIR = args.dir
    print('DIR = %s' % DIR)

    if DIR == '':
        sys.exit('ERROR: DIRECTORY with SV calls required.')
    if args.type == "four":
        four_graphs(DIR)
    elif args.type == "sixteen":
        sixteen_graphs(DIR)
    elif args.type == "simple":
        simple(DIR)
    else:
        print("Unknown command")


