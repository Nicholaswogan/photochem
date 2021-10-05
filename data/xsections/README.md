
# Cross Section Data

This directory contains cross section data for photolysis reactions. Data is in the folders, and citations for the data are in `metadata.yaml`.

The python script `check_consistency.py` checks for consistency between the metadata and the data in the folders.

# Notes

## 9/7/21

Kevin gave me his list of photolysis reactions, cross section data, and quantum yield data. All are located in the folder `Zahnle_Kevin_data`.

Here, I implemented Kevin Zahnle's cross sections for most species. Here are some of Kevin's notes on this cross sections.

- *C4H2* xs is 2 times C2H2 xs.
- *HSO* set to HO2. Kevin says "HSO is more stable than HO2, but, you know, photolysis never really matters for these short-lived species"
- *S8L*. Estimate for cross section of linear S8 chains in visible. set equal to S8R in uv
- *CH3OH*. CH3OH + hv  => CH3O +  H this branch is reported 86±10%. The other processes are undetected
- *CH3CN*. CH3CN + HV  -> CH3 +  CN the H + CH2CN is in fact somewhat more important

I did not use Kevin's cross sections for the following species: 'CO2', 'H2CO', 'O2', 'O3', 'SO2', and 'NO'. Instead I use the Phidrates database.

I omit the following photolysis reactions

```yaml
- equation: SO2 + hv => S + O + O
  type: photolysis
- equation: HNO + hv => NO + H
  rate-constant:
    A: 1.0e-3
    b: 0.0
    Ea: 0.0
- equation: HNO2 + hv => NO + OH
  rate-constant:
    A: 1.0e-3
    b: 0.0
    Ea: 0.0
```

I omit `SO2 + hv => S + O + O` because Phidrates has no data for this reaction. I omit `HNO + hv => NO + H` and `HNO2 + hv => NO + OH` because Kevin hard-codes rates instead of using cross sections. I am unable to find cross sections. 

## 10/3/21

I noticed that when modeling Mars, CO2 photolysis was much higher than previous models. It looks like the Phidrates cross sections are too large for CO2 + hv => CO + O branch. So I switched to Kevin Zahnle's CO2 cross sections. NOW i use Kevin's cross sections for all but

'H2CO', 'O2', 'O3', 'SO2', 'NO'

Which uses Phidrates.

## 10/5/21

I added `metadata.yaml`, which contains citations for all the cross sections. I also added `check_consistency.py` which check for consistency between data and metadata.


