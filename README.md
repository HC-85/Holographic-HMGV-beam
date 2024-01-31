# Holographic Helical Mathieu Gauss Vector Modes
This was done as part of the final project for an undergraduate optics course. 

Using _Experimental generation of Helical
Mathieu-Gauss vector modes_ as main source, a computational simulation for the production of holograms was implemented, allowing for a new experimental
generation of such modes.

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
## Computational Simulation
The characteristic values and expansion coefficients of the desired Mathieu
functions were calculated using [E. Cojocaru's MATLAB toolbox](https://www.mathworks.com/matlabcentral/fileexchange/22081-mathieu-functions-toolbox-v-1-0).

Both angular and radial functions of the given parity are multiplied a Gaussian beam and normalized.

