
# Henrys law coefficients

Henrys law coefficients are all stored in `henry.yaml`. The coefficient is computed in the following way. `H` has SI units, [mol/(kg Pa)]. Below T is in Kelvin.

```python
H = A*np.exp(B*(1.0/298.15 - 1.0/T))
```

All coefficients are from the Sanders 2015 database: https://doi.org/10.5194/acp-15-4399-2015 . They attach all the data in the supplementary materials. The important data are in `HbpSI.f90` and `henry.bib`, in this folder. This database has multiple entries for each molecule. I use the most recently published entry.

