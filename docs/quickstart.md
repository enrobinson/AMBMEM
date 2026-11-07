# Quick start

```bash
conda activate ambmem
cp /path/to/AMBMEM/examples/membrane_only_control.yaml build.yaml
ambmem validate build.yaml
ambmem build build.yaml
```

Require the safety marker before minimization or MD:

```bash
test -f output_directory/SIMULATION_READY \
  && echo "PASS: simulation-ready" \
  || echo "FAIL: do not simulate"
```

The requested lateral size is approximate. AMBMEM reports the selected commensurate periodic box and realized integer composition in `report.json`.
