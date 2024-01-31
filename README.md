# Holographic Helical Mathieu Gauss Vector Modes
This was done as part of the final project for an undergraduate optics course. 

<p align="center">
  <img src="/images/hologram.PNG" />
</p>

Using _Experimental generation of Helical
Mathieu-Gauss vector modes_ as main source, a computational simulation for the production of holograms was implemented, allowing for a new experimental
generation of such modes.

<p align="center">
  <img src="/images/experimental_vs_simulation.PNG" />
</p>

For the full report, see: [Mathematical simulation and experimental generation of Helical
Mathieu-Gauss vector modes](https://drive.google.com/file/d/1Vz3pPY1zT6H6XSvn4EKt_pxGjdzZ9dft/view).

## Brief Background
The solutions of the Helmholtz equation in a cylindrical orthogonal coordinate system yield distinct families of non-diffracting beams. For the case of elliptic cylindrical coordinates, one obtains Mathieu beams, whose corresponding differential equations are:
```math
\begin{align}
\left[\frac{d^2}{d\eta^2} + (a - 2q\cos 2\eta)\right]\Theta(\eta) &= 0\\
\left[\frac{d^2}{d\xi^2} - (a - 2q\cosh 2\xi)\right]R(\xi) &= 0
\end{align}
```
where $q = \left({fk_{t}}/{2}\right)^2$, and $a$ is the major axis of the ellipse. Their solutions are the angular and radial Mathieu functions respectively.

Non-diffracting beams can be written in the form:
```math
\begin{align}
U(\mathbf{r}) = e^{-iKz}\text{GB}(\mathbf{r})W\left(\frac{x}{\mu},\frac{y}{\mu};k_{t}\right)
\end{align}
```
where $\text{GB}$ is the fundamental Gaussian beam and $W$ describes the transverse distribution. For our case, $W$ corresponds to the Mathieu functions and so we can construct helical Mathieu-Gauss vector beams as:
```math
\begin{equation}
    \textbf{HMGV}_{m_1, m_2} = \cos \theta \text{ HMG}^+_{m_1} \hat{e}_R + \sin \theta \text{ } e^{i \alpha}\text{ HMG}^-_{m_2} \hat{e}_L  
\end{equation}
```
which are non-separable superpositions of the spatial and polarization degrees of freedom, weighted by $\theta$ and $\alpha$ respectively. Superscripts $+$ and $-$ denote polarization handedness.


## Computational Simulation
The characteristic values and expansion coefficients of the desired Mathieu
functions were calculated using [E. Cojocaru's MATLAB toolbox](https://www.mathworks.com/matlabcentral/fileexchange/22081-mathieu-functions-toolbox-v-1-0).

Both angular and radial functions of the given parity are multiplied by a Gaussian beam and normalized.

To create the binary hologram, we define $\theta _{xy}$ and $\theta _z$, which measure the deviation from the $z$-axis in terms of an azimuthal and polar angle respectively. This deviation is expressed with an exponential function having dependency on these angles that is multiplied with the helical Mathieu-Gauss beam of interest.
The phase is encoded as follows:
```math
\begin{align}
\varphi &= \tan ^{-1}\left(\frac{\Im V}{\Re V}\right)\\
\phi &= \sin ^{-1}\left(\frac{|\varphi|}{\max |\varphi|}\right)\\
T &= \frac{1}{2} + \frac{1}{2}\text{sgn}(\cos\varphi + \cos \phi)
\end{align}
```

$T$ gives us the binary holograms, the negatives of which one has to print. 

For experimental generation of the vector modes, both scalar Helical Mathieu-Gauss modes ($+$/ $-$) have to interfere.
