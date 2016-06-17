
import os

temp_dir = "/tmp/sv_pipeline"
os.system("mkdir -p {}".format(temp_dir))
os.system("mkdir -p figs")

from .the_pipeline import *






