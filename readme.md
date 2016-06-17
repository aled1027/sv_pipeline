# Structual Variant Pipeline

## Notes
- If running pipeline via `ssh`, login with `ssh -X`

## Profiling
```
python -m cProfile -s time the_pipeline.py > profile.text 2>&1
```

## TODO
- separate out io calls, maybe into a separate module

