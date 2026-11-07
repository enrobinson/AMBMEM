# Installation

## Requirements

- AmberTools 26 or newer
- Python 3.10 or newer
- ParmEd from the AmberTools environment
- NumPy
- PyYAML
- pytest for development testing

CPPTRAJ and tLEaP are supplied by AmberTools. AMBMEM does not require OpenMM or sander.

## Install from GitHub

```bash
conda create -n ambmem -c conda-forge \
  ambertools=26 python=3.12 numpy pyyaml pytest
conda activate ambmem

git clone https://github.com/enrobinson/AMBMEM.git
cd AMBMEM
python -m pip install .
```

For development, use `python -m pip install -e .`.

Verify with:

```bash
ambmem --version
ambmem doctor
python -m pytest
```
