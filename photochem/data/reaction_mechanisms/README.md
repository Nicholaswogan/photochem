`zahnle_earth.yaml` is Kevin Zahnle's network of reactions. He calls the network `earth_125.rx`. The meaning of "125" is unclear. I have documented the small changes to the network [at this Github repo.](https://github.com/Nicholaswogan/ImpactAtmosphere/blob/main/ImpactAtmosphere/data/notes.md). But from now, onward, I will document the changes here.

# 10/27/21

- Several species had valid polynomial Gibbs energy fits down to 100 K. But this was not compatible with Titan, because Titan is really cold. Therefore, I decreased the lower bounds of the fits to 50 K.
- I added many particles. Particles are made in various way from species. 

I added

```yaml
- equation: C2H4 + N2D <=> CH3CN + H
  rate-constant: {A: 2.3e-10, b: 0.0, Ea: 503}
```

I was modeling Titan, and was getting far to little CH3CN. Lavvas et al. (2008) uses this reaction rate.

I will keep the reaction for now, but note that this reaction might not be correct. From Loison et al (2015),

> The CH3CN molecule is not produced in our model by the N(2D) + C2H4 reaction but through the association reaction H+CH2CN (Balucani et al., 2012; Lee et al., 2011). We calculated the rate constant from our semi-empirical model for association reactions (Hébrard et al., 2013), the CH2CN radical being produced by the H + C2H4CN reaction (this work, C2H4CN itself being the product of the H + C2H3CN association) and the N + C2H3 reaction (Payne et al., 1996).

I would like to follow Loison et al.'s advice, but our model is missing two things

1. Radiative association reactions: A + B => AB + hv. These reactions are important for some speices in Titan's atmosphere (Vuitton et al. 2012)
2. We do not have C2H4CN in the network.

# 12/16/21

I added

```yaml
- equation: Cl + O3 <=> ClO + O2
  rate-constant: {A: 2.8e-11, b: 0.0, Ea: 250}
```

From Atkinson et al. (2007). Full citation:

> R. Atkinson, D. L. Baulch, R. A. Cox, J. N. Crowley, R. F. Hampson, et al.. Evaluated kinetic and photochemical data for atmospheric chemistry: Volume III ? reactions of inorganic halogens. Atmospheric Chemistry and Physics Discussions, European Geosciences Union, 2006, 6 (2), pp.2281-2702. ⟨hal-00301101⟩

# 12/23/21

I updated the following reaction. Rate from JPL-15

```yaml
- equation: S + O2 <=> SO + O
  rate-constant: {A: 1.6e-12, b: 0, Ea: -100.0}
```

Also added NH3aer and He.


