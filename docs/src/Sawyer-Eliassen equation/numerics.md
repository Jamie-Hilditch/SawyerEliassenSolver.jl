# Numerical Solution

The Sawyer-Eliassen equation and the accompanying ``v`` and ``b`` equations can be written in the form

```math
\begin{align}
    \psi_{tt} & = -\mathcal{L}_\psi\psi + \nabla^{-2}\mathfrak{F} \\
    v_t & = \mathcal{L}_v\psi + \mathcal{F}^{(y)} \\ 
    b_t & = \mathcal{L}_b\psi + \mathcal{B}
\end{align}
```
where ``\mathcal{L}_\psi``, ``\mathcal{L}_v`` and ``\mathcal{L}_b`` are spatially varying but constant in time linear differential operators. The forcing functions ``\mathfrak{F}``, ``\mathcal{F}^{(y)}`` and ``\mathcal{B}`` can be functions of space and time.