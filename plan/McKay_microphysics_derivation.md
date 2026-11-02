# McKay Microphysics

## Introduction

Here, I derive a microphysical model that is motivated by the one used in [McKay et al. (1989)](https://doi.org/10.1016/0019-1035(89)90160-7) for Titan. The main difference with the McKay model is that this derivation assumes that new haze production generates monomers. McKay, on the other hand, assumed that all new haze production added mass to existing particles. I made this change because in the McKay approach, you have to make assumptions about the upper boundary conditions that seem hard to justify for an exoplanet model. The model presented here has the following assumptions:

- The atmosphere is one-dimensional.
- The model is steady state.
- Haze mass is produced at a rate $q$, and production introduces particles of a prescribed monomer radius $r_0$.
- The only haze transport mechanism is sedimentation. We neglect any lofting (i.e. through eddy diffusion).
- Particle radius evolves from coagulation and haze production. There is no mechanism for mass removal from aerosols (i.e. evaporation).
- At a given altitude, particles are represented by a single radius (i.e., a monodisperse approximation).
- Aerosols are compact spherical particles with constant material density.

I'm going to generally use the McKay et al. (1989) notation. All variables and their definitions are in a table at the end of the document.

## Main Derivation

We begin our derivation with a statement of mass conservation.

$$\frac{\partial \rho_h}{\partial t} + \frac{\partial \Phi_m}{\partial z} = q$$

Here, $\rho_h$ is the haze mass per volume in the atmosphere. It is altered by the flux of mass entering/leaving that volume ($\Phi_m$) and the production or loss of mass, $q$. We can also make a statement of particle conservation:

$$\frac{\partial n}{\partial t} + \frac{\partial \Phi_N}{\partial z} = S_N$$

Here, the particles/cm$^3$ change in time because of the flux $\Phi_N$ and the source term $S_N$. If we assume that the source term is only production and coagulation, then

$$S_N = -Kn^2 + J$$

where $J$ is the production rate of new particles and $Kn^2$ is the loss rate of particle number through coagulation. McKay assumes $J=0$, so production adds mass to the existing particle population without independently creating new particles. Here, we retain $J$ until specifying the monomer source. Note that in a later subsection we describe why coagulation is $Kn^2$, and derive an expression for $K$ based on Brownian motion.

Now, assume steady state, so $\partial\rho_h/\partial t=\partial n/\partial t=0$. Also, we define positive downward mass and particle-number fluxes, $C=-\Phi_m$ and $F=-\Phi_N$. The continuity equations become

$$\frac{dC}{dz} = -q$$
$$\frac{dF}{dz} = Kn^2 - J$$

Integrating mass conservation from altitude $z$ to the upper boundary gives

$$C(z) = C(z_{\rm top}) + \int_z^{z_{\rm top}} q(z') dz'$$

For a given specified production rate $q$ as a function of altitude, and the incoming mass flux at the top of the model domain $C(z_{\rm top})\geq0$ (often set to zero), we can easily compute $C$ everywhere in the atmosphere through the above integral.

We now make a key approximation: **all particles at a given altitude are represented by the same radius $r(z)$**. We also treat them as compact spheres with constant material density $\rho$. Thus, the representative particle has mass

$$m=\frac{4\pi}{3}\rho r^3,$$

and so the haze mass density is

$$\rho_h=nm$$

Further, we assume that particles are transported only from falling (no lofting from eddy diffusion). Given that all particles have the same radius, then they also have the same settling velocity $v$. Therefore,

$$F=-nv, \qquad C=-nmv$$

These imply the useful relation

$$C=mF$$

Now, we further assume that all produced haze mass enters as new particles of a prescribed monomer radius $r_0$

$$m_0=\frac{4\pi}{3}\rho r_0^3.$$

For nonnegative mass production $q$, the corresponding particle-number production is

$$J=\frac{q}{m_0}$$

The mass of these new particles is already included in $q$; $J$ is not an additional mass source. The number conservation equation is then

$$\frac{dF}{dz} = Kn^2 - \frac{q}{m_0}$$

Substituting in $n=-F/v$,

$$\frac{dF}{dz}=K\frac{F^2}{v(r)^2}-\frac{q}{m_0}.$$

We can simply integrate the above equation downward from the top of the atmosphere. Here, $C$ is known from the production integral. Wherever haze is present, we can compute the radius with:

$$r=\left(\frac{3C}{4\pi\rho F}\right)^{1/3}$$

The above relation comes from combining $C=mF$ and $m=\frac{4\pi}{3}\rho r^3$. We evaluate the above formula for radius during integration with the current $C$ and $F$, then use that radius to compute settling and coagulation. Also, particle number density can be derived with $n=-\frac{F}{v(r)}$.

At the upper boundary, we assume that any incoming mass flux at the top of the model domain consists of monomers, giving the number-flux boundary condition

$$F(z_{\rm top})=\frac{C(z_{\rm top})}{m_0}.$$

A positive incoming monomer flux gives $r(z_{\rm top})=r_0$. With zero incoming flux, $C/F\to m_0$ at the onset of production, so $r\to r_0$. In any empty region above production, $C=F=n=\rho_h=0$; the monomer radius $r_0$ can be used as a radius placeholder without introducing particles.

## Fall velocity

To estimate the fall velocity of aerosols, we use Stokes' law with a Cunningham correction:

$$v = -\frac{2(\rho-\rho_a)gr^2}{9\eta}C_c$$

Here, $\eta$ is the gas dynamic viscosity. We will use a temperature-dependent viscosity parameterization for Earth air. The correction term $C_c$ is then,

$$C_c = 1+\mathrm{Kn}\left[1.257+0.4\exp\left(-\frac{1.1}{\mathrm{Kn}}\right)\right], \qquad \mathrm{Kn}=\frac{\lambda_{\rm gas}}{r}$$

Here, $\lambda_{\rm gas} = 1/(n_a A_{\rm eff})$, with $A_{\rm eff}=6.0\times10^{-15}\ \mathrm{cm^2}$ is the effective molecular collision area that we adopt. This differs from McKay's original settling prescription, which blends a free-molecular settling coefficient toward a continuum expression rather than using the Cunningham correction.

## Coagulation

Brownian motion brings particles into contact, and we assume they stick when they touch. First define an uncharged pair-encounter kernel $K_{\rm enc}$, with units cm$^3$/s. This is different from McKay's particle-number loss coefficient $K$, which we explain below.

In the **free-molecular limit**, typically high up in the atmosphere, particles move approximately ballistically between encounters. For identical spheres, the collision cross section is $\pi(2r)^2$, and the mean relative thermal speed is $\sqrt{16k_BT/(\pi m)}$. Their product gives

$$K_{\rm fm} = \pi(2r)^2\sqrt{\frac{16k_BT}{\pi m}} = 8\sqrt{\frac{3k_BT}{\rho}}\sqrt r,$$

where we substituted $m=(4\pi/3)\rho r^3$.

In the **continuum limit**, typically lower in the atmosphere, encounters are controlled by particle diffusion. For two identical spheres with diffusion coefficient $D$, the encounter kernel is $4\pi(r+r)(D+D)$. Using the Stokes-Einstein relation $D=k_BT/(6\pi\eta r)$ gives

$$K_{\rm cont} = 16\pi rD = \frac{8k_BT}{3\eta}.$$

We blend these limits using the Knudsen number $\mathrm{Kn}=\lambda_{\rm gas}/r$:

$$K_{\rm enc} = \frac{\mathrm{Kn}\,K_{\rm fm}+K_{\rm cont}}{1+\mathrm{Kn}}.$$

This approximate interpolation approaches the free-molecular kernel at large $\mathrm{Kn}$ and the continuum kernel at small $\mathrm{Kn}$. McKay's Fortran code uses this blend for $\mathrm{Kn}<100$ and the free-molecular kernel alone for $\mathrm{Kn}\geq100$.

Each particle encounters others at a rate $nK_{\rm enc}$, with units s$^{-1}$: this is the encounter rate per particle. Multiplying by the number density gives $K_{\rm enc}n^2$, but counts every pair twice: once from each particle's perspective. The number of distinct encounters per unit volume per unit time is therefore $K_{\rm enc}n^2/2$, with units encounters cm$^{-3}$ s$^{-1}$. Each merger turns two particles into one, so it removes one particle from the count. Thus, the corresponding particle-number loss rate has units particles cm$^{-3}$ s$^{-1}$, and the uncharged number-loss coefficient is $K=K_{\rm enc}/2$. Particle and encounter counts are treated as dimensionless in the dimensional analysis; the kernel itself has units cm$^3$ s$^{-1}$.

Finally, similarly charged particles repel one another, reducing the rate at which they reach contact. Following McKay, we crudely parameterize this inhibition with

$$f_{\rm charge}=\exp(-r/r_c),$$

where the prescribed charging radius $r_c$ is the particle radius at which the rate is reduced by a factor of $e$. McKay nominally adopts $r_c = 0.09$ $\mu$m for Titan's atmosphere. Combining the blend, pair-counting factor, and charge suppression gives

$$K = \frac{e^{-r/r_c}}{2}\frac{\mathrm{Kn}\,K_{\rm fm}+K_{\rm cont}}{1+\mathrm{Kn}}$$

The particle-number sink in the continuity equation is then $-Kn^2$.

## Notation

| Variable | Description | Units |
|---|---|---|
| $t$ | Time. | s |
| $z$ | Altitude, increasing upward. | cm |
| $z_{\rm top}$ | Upper boundary of the haze model. | cm |
| $r(z)$ | Representative particle radius used in the local single-radius approximation. | cm |
| $n(z)$ | Particle number density. | particles cm$^{-3}$ |
| $\rho$ | Material density of the compact spherical particles, taken to be constant. | g cm$^{-3}$ |
| $m(r)$ | Mass of one particle, $m=(4\pi/3)\rho r^3$. | g particle$^{-1}$ |
| $\rho_h(z)$ | Haze mass per unit atmospheric volume, $\rho_h=nm$; distinct from particle material density $\rho$. | g cm$^{-3}$ |
| $q(z)$ | Local haze mass production rate per unit atmospheric volume. | g cm$^{-3}$ s$^{-1}$ |
| $S_N(z)$ | Net particle-number source rate per unit atmospheric volume, including production and loss. | particles cm$^{-3}$ s$^{-1}$ |
| $J(z)$ | Production rate of new particles per unit atmospheric volume. | particles cm$^{-3}$ s$^{-1}$ |
| $v(z)$ | Signed particle settling velocity; negative because particles fall toward decreasing $z$. | cm s$^{-1}$ |
| $\Phi_m(z)$ | Signed upward haze mass flux, $\Phi_m=nmv$; negative for sedimentation. | g cm$^{-2}$ s$^{-1}$ |
| $\Phi_N(z)$ | Signed upward particle-number flux, $\Phi_N=nv$; negative for sedimentation. | particles cm$^{-2}$ s$^{-1}$ |
| $C(z)$ | Positive downward haze mass flux, $C=-\Phi_m$; incoming mass flux plus integrated mass production above $z$. | g cm$^{-2}$ s$^{-1}$ |
| $F(z)$ | Positive downward particle-number flux, $F=-\Phi_N$. | particles cm$^{-2}$ s$^{-1}$ |
| $r_0$ | Prescribed monomer radius, used for newly produced and incoming particles. | cm |
| $m_0$ | Mass of one monomer, $(4\pi/3)\rho r_0^3$; relates number and mass sources through $J=q/m_0$. | g particle$^{-1}$ |
| $K$ | Coagulation number-loss coefficient: the particle-number sink is $Kn^2$. Includes charge suppression when applied. Particle counts are treated as dimensionless in the coefficient's units. | cm$^3$ s$^{-1}$ |
| $K_{\rm enc}$ | Uncharged Brownian pair-encounter kernel; with perfect sticking and no charge suppression, $K=K_{\rm enc}/2$. | cm$^3$ s$^{-1}$ |
| $K_{\rm fm}$ | Uncharged pair-encounter kernel in the free-molecular limit. | cm$^3$ s$^{-1}$ |
| $K_{\rm cont}$ | Uncharged pair-encounter kernel in the continuum limit. | cm$^3$ s$^{-1}$ |
| $D$ | Particle diffusion coefficient in the continuum limit. | cm$^2$ s$^{-1}$ |
| $T(z)$ | Atmospheric temperature. | K |
| $n_a(z)$ | Atmospheric gas number density. | molecules cm$^{-3}$ |
| $\rho_a(z)$ | Atmospheric gas mass density, used for buoyancy in the settling formula. | g cm$^{-3}$ |
| $g(z)$ | Magnitude of gravitational acceleration. | cm s$^{-2}$ |
| $k_B$ | Boltzmann constant; written $k$ in McKay. | erg K$^{-1}$ |
| $\eta(z)$ | Gas dynamic viscosity, needed for the continuum limits of settling and coagulation. | g cm$^{-1}$ s$^{-1}$ |
| $\lambda_{\rm gas}(z)$ | Gas mean free path, $\lambda_{\rm gas}=1/(n_a A_{\rm eff})$. | cm |
| $A_{\rm eff}$ | Effective molecular collision area used to estimate the mean free path; adopted value $6.0\times10^{-15}$ cm$^2$. | cm$^2$ |
| $\mathrm{Kn}$ | Knudsen number, defined here as $\lambda_{\rm gas}/r$, used to interpolate between free-molecular and continuum limits. | dimensionless |
| $C_c$ | Cunningham slip correction factor, which reduces drag relative to the no-slip Stokes expression. | dimensionless |
| $r_c$ | Prescribed charging radius: the coagulation coefficient is multiplied by $\exp(-r/r_c)$. It parameterizes electrostatic inhibition of encounters. | cm |
| $f_{\rm charge}$ | Charge-suppression factor, $\exp(-r/r_c)$. | dimensionless |
