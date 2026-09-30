# Calculation code

Start with the [repository guide](../README.md#reproducing-the-calculations).
The standard-library-only [reproduction driver](../tools/reproduce.py) supplies
the documented configuration paths and parameter grids.

```bash
python3 tools/reproduce.py check
python3 tools/reproduce.py primary
python3 tools/reproduce.py focused
python3 tools/reproduce.py sensitivity
```

Run these commands from the repository root. New output goes to `reproduced/`,
never to the archived [Results](../Results) directory. The main guide contains
the [script-by-script reference](../README.md#script-reference).

NumPy optionally accelerates the exact factor contractions. Without NumPy the
original sparse implementation is used; neither backend uses sampling or an
approximate likelihood. Install the acceleration with `python3 -m pip install numpy`.

The [provenance check](../tools/result_provenance.py) detects changed calculation
code, JSON inputs or numerical outputs. Only record a new calculation after
reproduction and scientific review, not merely to silence a failed check.
