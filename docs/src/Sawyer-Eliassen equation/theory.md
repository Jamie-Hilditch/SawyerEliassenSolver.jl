# The Sawyer-Eliassen equation 

The Sawyer-Eliassen equation [Sawyer_Sutton_1956,Eliassen_1962](@citep) describes the ageostrophic overturning circulations at fronts. The derivation here largely follows [Mooers_1975](@citet) although we allow for arbitary momentum and buoyancy forcing. Consider a 2D set-up with a buoyancy field ``B(x,z)`` in thermal wind balance ``\partial B/\partial x = f\partial V/\partial z`` with a geostrophic velocity ``V(x,z)``. Assuming that perturbations from this state are also 2D the Boussinesq equations can be written.

```math
\begin{align}
    \frac{\partial u}{\partial t} - fv + \frac{1}{\rho_0}\frac{\partial p}{\partial x} & = \mathcal{F}^{(x)} \\     
    \frac{\partial v}{\partial t} + u\frac{\partial V}{\partial x} + w\frac{\partial V}{\partial z} + fu & = \mathcal{F}^{(y)} \\
    \frac{\partial w}{\partial t} + \frac{1}{\rho_0}\frac{\partial p}{\partial z} - b & = \mathcal{F}^{(z)} \\
    \frac{\partial b}{\partial t} + u\frac{\partial B}{\partial x} + w\frac{\partial B}{\partial z} & = \mathcal{B} \\
    \frac{\partial u}{\partial x} + \frac{\partial w}{\partial z} & = 0
\end{align}
```
where the non-linear terms have been absorbed into the arbitary RHS forcing. Introducing a streamfunction ``\psi`` such that ``u = -\partial\psi/\partial z`` and ``w = \partial\psi/\partial x`` we form the evolution equation for the component of vorticity into the page.

```math
\begin{equation}
    \left(\frac{\partial^2 ~}{\partial x^2} + \frac{\partial^2 ~}{\partial z^2}\right)\frac{\partial\psi}{\partial t} = \frac{\partial b}{\partial x} - f\frac{\partial v}{\partial z} + \frac{\partial \mathcal{F}^{(z)}}{\partial x} - \frac{\partial \mathcal{F}^{(x)}}{\partial z}
\end{equation}
```

and the ``v`` and ``b`` perturbations are determined by 
```math
\begin{align}
    \frac{\partial v}{\partial t} & = -\mathcal{J}(\psi,V + fx) + \mathcal{F}^{(y)} \\ 
    \frac{\partial b}{\partial t} & = -\mathcal{J}(\psi,B) + \mathcal{B}
\end{align}
```
where ``\mathcal{J}(\psi, \cdot) \equiv (\partial{\psi}/\partial x) \partial/ \partial z - (\partial{\psi}/\partial z) \partial/ \partial x`` is advection by the perturbations.

The dynamics are controlled by a single second-order in time equation for the stream-function
```math
\begin{equation}
    \left(\frac{\partial^2 ~}{\partial x^2} + \frac{\partial^2 ~}{\partial z^2}\right)\frac{\partial^2\psi}{\partial t^2} = -\frac{\partial ~}{\partial x}\mathcal{J}(\psi,B) + f\frac{\partial ~}{\partial z}\mathcal{J}(\psi,V + fx) + \mathfrak{F}(x,z,t)
\end{equation}
```
where the forcing is 
```math
\begin{equation}
    \mathfrak{F}(x,z,t) = \frac{\partial^2 \mathcal{F}^{(z)}}{\partial x\partial t} - \frac{\partial^2 \mathcal{F}^{(x)}}{\partial z\partial t} - f\frac{\partial\mathcal{F}^{(y)}}{\partial z} + \frac{\partial\mathcal{B}}{\partial x}
\end{equation}
```
In expanding the RHS some cancellations can be made thanks to the assumption of thermal wind balance and the Sawyer-Eliassen equation is reached
```math
\begin{equation}
    \left(\frac{\partial^2 ~}{\partial x^2} + \frac{\partial^2 ~}{\partial z^2}\right)\frac{\partial^2\psi}{\partial t^2} = -\frac{\partial B}{\partial z}\frac{\partial^2\psi}{\partial x^2} + 2\frac{\partial B}{\partial x}\frac{\partial^2\psi}{\partial x\partial z} - f\left(f + \frac{\partial V}{\partial x}\right)\frac{\partial^2\psi}{\partial z^2} + \mathfrak{F}(x,z,t)
\end{equation}
```
