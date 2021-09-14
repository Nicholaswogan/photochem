

# particle plan

We can deal with two problems sort-of separately

1. Problem of making aerosols
2. the problem of evolving their size distributions

Lets only worry about (1), and leave (2) for a rainy day. So we will assume fixed particle sizes for now.

We can not actually separate these problems completely. Modeling full particle size distributions, and proper physics, ties these two things together. Because species can condense + evaporate onto and off of aerosols.

Ultimately this can all get really complicated.

Suppose we have a single aerosol, which multiple species can condense onto. We must track the full distribution of particle sizes, and also the composition as a function of size. So this adds (number of size bins)*(number of atoms in aerosol) new equations per altitude. This should be (20 bins)(4 atoms) = 80 new things to track, which is insane.

Let's look at Seinfeld (2006), to gain intuition about this problem.

```yaml
particles:
  # approach used in photochempy
  - name: HCAER
    formation: reaction 
    composition:
      C: 6
      H: 4
    density: 1.0 # g/cm3
    optical-properties: khare1984
    optical-type: mie
    
    # approach for S in photochempy
  - name: SO4AER
    formation: saturation
    composition:
      S: 1
      O: 4
    gas-phase: SO4
    saturation-parameters: ??? 
    density: 1.0 # g/cm3
    optical-properties: khare1984
    optical-type: mie
    
    # Lavvas approach
  - name: haze
    formation: reactions
    fixed-composition: false
    track-composition: false
    density: 1.0 # g/cm3
    monomer-radius: 7.35e-4 # microns
    reactions:
    - equation: H2CN + HCN => polymer
      rate-constant:
        A: 1.1e-12
        b: 0.0
        Ea: 900
    


  # zahnle approach
  - name: HCAER2
    formation: multiple-saturation
    gas-phase:
      - name: C4H2
        condensation-rate: 1e-3
        floor: 1e-8
      - name: C4H4
        condensation-rate: 1e-3
        floor: 1e-8
      - name: HCCCN
        condensation-rate: 1e-3
        floor: 1e-8
      - name: HCN
        condensation-rate: 1e-3
        floor: 1e-8
    density: 1.0 # g/cm3
    optical-properties: khare1984
    optical-type: mie
    
    
```

I will take the Lavvas approach. Except initially, I will not consider microphysics. I will assume fixed particle sizes vs altitude.


Lets start simple with sulfur condensation.

```yaml
- name: H2SO4aer
  formation: saturation
  density: 1.0 # g/cm3
  optical-properties: khare1984
  optical-type: mie
  composition:
    H: 2
    S: 1
    O: 4
  gas-phase: H2SO4
  saturation-parameters:
    A: 16.722
    B: -10156.0
    C: 0.0
  
```








