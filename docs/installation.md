# Installation

```bash
conda create -n ambmem -c conda-forge ambertools=26 python=3.12 numpy pyyaml pytest
conda activate ambmem
git clone https://github.com/enrobinson/AMBMEM.git
cd AMBMEM
python -m pip install .
ambmem doctor
```

AmberTools 26 or newer supplies `tleap`, ParmEd, Lipid21, OPC, ions, and GAFF2. AMBMEM supplies validated coordinate references for packing.
