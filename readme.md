# Structual Variant Pipeline

## Notes
- If running pipeline via `ssh`, login with `ssh -X`


## Profiling
```
python -m cProfile -s time the_pipeline.py > profile.text 2>&1
```

## TODO
- double check that four pipeline works with bioinf0
- optimize `the_pipeline.py`
    - Reduce the number of times that we need to read from disk

