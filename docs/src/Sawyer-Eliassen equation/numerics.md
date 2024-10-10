# Numerical Solution

The Sawyer-Eliassen equation can be written in the form
```math
\begin{equation}
    \zeta_{tt} = -\mathcal{L}\zeta + F \\
\end{equation}
```
with the accompanying ``v`` and ``b`` equations
```math
\begin{align}
    v_t & = -u(f + V_x) - wV_z + \mathcal{F}^{(y)} \\ 
    b_t & = -uB_x - wBz + \mathcal{B}
\end{align}
```
where 
```math
\begin{equation*}
        u = -\psi_z, \quad w = \psi_x, \quad \psi = \nabla^{-2}\zeta 
\end{equation*}
```

## Pseudospectral discretisation ##
The Sawyer-Eliassen operator ``\mathcal{L}`` varies in space but not time whereas the forcings are functions of space and time. We use a pseudospectral discretisation solving for ``\zeta`` in spectral space. Enforcing the boundary conditions in ``z`` we express ``\zeta`` in Fourier-sine modes. ``u`` and ``w`` are then computed in Fourier-cosine and Fourier-sine space respectively. Products between the background flow and waves are computed in physical space and a user prescribed number of the highest wavenumber modes are zeroed-out in order to dealias the solution. ``v`` and ``b`` only ever exist in physical space.

## Diagonally Implicit Runge-Kutta Nyström ##

We utilise a 2-stage 3rd-order accurate diagonally implicit Runge-Kutta Nyström (DIRKN) timestepping scheme due to [Sharp_Fine_Burrage_1990](@citet). The state variables are advanced according to 
```math
\begin{align}
    \zeta^{n+1} & = \zeta^n + h \zeta_t + h^2 \sum_{j=1}^2 b_j \zeta_{tt}^{n+c_j} \\ 
    \zeta_t^{n+1} & = \zeta_t^n  + h \sum_{j=1}^2 b_j^\prime \zeta_{tt}^{n+c_j} \\ 
    v^{n+1} & = v^n + h \sum_{j=1}^2 b_j^\prime v_t^{n+c_j} \\ 
    b^{n+1} & = b^n + h \sum_{j=1}^2 b_j^\prime b_t^{n+c_j}
\end{align}
```
with the intermediate stages defined by the implicit equations 
```math
\begin{equation}
    \zeta^{n + c_j} = \zeta_n + c_j h \zeta_t + h^2 \sum_{k=1}^j a_{jk} \zeta_{tt}^{n + c_k}
\end{equation}
```
To solve, we define the implicit Sawyer-Eliassen operator ``\mathcal{L}^I \equiv I + h^2 a_{ii} \mathcal{L}``. As is typical of DIRKN schemes we take the diagonal elements ``a_{11}`` and ``a_{22}`` to be the same. Then we proceed by solving 
```math 
\begin{equation}
    \mathcal{L}^I\zeta^{n + c_1} = \zeta^n + c_1 h \zeta_t^n + h^2 a_{11} F^{n + c_1}
\end{equation}
```
and then 
```math 
\begin{equation}
    \mathcal{L}^I\zeta^{n + c_2} = \zeta^n + c_2 h \zeta_t^n + h^2 a_{21} \left( - \mathcal{L}\zeta^{n + c_1} + F^{n + c_1}\right) + h^2 a_{22} F^{n + c_2}
\end{equation}
```

### Coefficients ###

[Sharp_Fine_Burrage_1990](@citet) derive conditions under which the numerical scheme is 3rd order accurate. The coefficients form a one-parameter family which we write in terms of ``c = c_1``.
```math
\begin{equation}
\begin{gathered}
    c_2 = \frac{3c - 2}{3(2c-1)}, \quad
    a_{11} = a_{22} = \frac{1}{2}c^2, \quad
    a_{21} = \frac{-2(9c^4 - 9c^3 + 3c - 1)}{9(2c-1)^2}, \\
    b_1 = \frac{1 - c}{4(3c^2 - 3c + 1)}, \quad 
    b_2 = \frac{(3c-1)(2c-1)}{4(3c^2 - 3c + 1)}, \\ 
    b_1' = \frac{1}{4(3c^2 - 3c + 1)}, \quad 
    b_2' = \frac{3(2c-1)^2}{4(3c^2 - 3c + 1)}    
\end{gathered}
\end{equation}
```
The numerical scheme is only R-stable (unconditionally stable when ``\mathcal{L}`` is positive definite) for certain ranges of ``c``. By default we use ``c = 17/14``.
