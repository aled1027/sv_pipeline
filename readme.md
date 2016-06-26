# Structual Variant Pipeline

## Usage

The basic pipeline may be ran with the follow command. 

```
python run.py four data/ambig_calls/ figs/all.pdf --L 600
```

The first argument specifies the type of pipeline that is being run. 
In this case `four` indicates that the pdf with the four plots per page should be generated.
Other options include `simple`, which does (TODO), and `sixteen` which generates a pdf page with 16 subplots, where each subplot is the overlap graph for a different `L`.

The second argument, `data/ambig_calls`, specifies the folder where the input files are found.
The `four` pipeline requires 3 files per region: `x-refcoords.bed`, `x_merged.m4`, `x.fa`, where `x` is the name of the region.
The pipeline uses all `*.m4` files in the specified directory as regions to examine.

The third argument specifies the output pdf. 
If pdftk is not installed (and in which case the pipeline will run correctly), then `figs/all.pdf` will not exist, but `figs/` will be populated with a pdf for each region.

The fourth argument, `--L`, is an optional argument specifying the minimum matching length to say that a read overlaps another read. The `L` parameter is fed into `minimap`; see `minimap` for more information.

## Notes
- If running pipeline via `ssh`, login with `ssh -X`

## Dependencies
- See requirements.txt for python dependencies
- Pdftk
    - Not necessary, but it combines many pdfs into a single pdf.
    - On debian, install with `sudo apt-get install pdfkt`

## Profiling
```
python -m cProfile -s time the_pipeline.py > profile.text 2>&1
```

