

# particle plan

In the reaction file, there will be a new entry for particles. There will be two ways, for now to make particles

1. reaction
2. saturation of gas phase

I need to revisit the math for particle size evolution.

```yaml
particles:
  - name: HCAER
    composition:
      C: 6
      H: 4
    optical-properties: khare1984
    optical-type: mie
    formation: reaction
    nucleation-size: 0.001 # microns
    
  - name: SO4AER
    composition:
      S: 1
      O: 4
    optical-properties: khare1984
    optical-type: mie
    formation: saturation
    gas-phase: SO4
    nucleation-size: 0.001 # microns
    saturation-parameters: ??? 
```

Then in the settings file, we might add the ability to optionally overwrite the optical properties? Maybe leave this for later.

```yaml
particles:
  - name: HCAER 
    optical-properties: khare1984
    type: mie 
```

