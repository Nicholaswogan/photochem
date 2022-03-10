# Rayleigh scattering cross sections

Main secondary source is [Ranjan and Sasselov (2017)](https://doi.org/10.1089/ast.2016.1519). Primary sources are listed in `rayleigh.yaml`

## Formalisms

 `vardavas`

From [Varavadas and Carver (1984)](https://doi.org/10.1016/0032-0633(84)90074-6), Equation 26

$$\sigma_R = 4.577 \times 10^{-21} \frac{1}{(\lambda_{nm}10^{-3})^4} \left(\frac{6+3\Delta}{6-7\Delta}\right)  \left( A \left(1+\frac{B}{(\lambda_{nm}10^{-3})^2} \right) \right)^2$$

Where, $\sigma_R$ is in $\mathrm{cm}^2$, and $\Delta$ is called the depolarization factor. The terms $A$ and $B$ are empirical constants for the refractive index of a given molecule.

$$n-1 =  A \left(1+\frac{B}{(\lambda_{nm}10^{-3})^2} \right) $$

where $n$ is the refractive index.