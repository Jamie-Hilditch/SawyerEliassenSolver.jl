# Variably Stratified Frontal Zone #

## Background Flow ##

```math
\begin{aligned}
    V_x &= 0 \\
    B_x &= fV_z = \Gamma f^2 \\ 
    B_z &= \left[\Pi_0^2 + \left(\Pi_\infty^2 - \Pi_0^2\right)\frac{1}{1 + \mathrm{e}^{z + D}}\right]f^2
\end{aligned}
```

## Initial Conditions ##

The streamfunction is specified as a wavepacket of the form
```math
\psi = A\cos(kx + mz -\omega t)
```
where ``A`` is a slowly varying amplitude that is neglected when constructing the initial conditions. The initial amplitude is
```math
A(x,z) = \mathrm{e}^{-(\lambda_x^2(x - x_0)^2 + \lambda_z^2(z - z_0))/2}
```
The wavepacket is initialised in the top part of the domain where ``B_z = \Pi_0^2f^2``. ``k`` and ``m`` are chosen such that the wave has frequency ``\omega = \frac{1}{\sqrt{2}}f``. That is 
```math
\frac{m}{k} = \frac{\Gamma \pm \sqrt{\Gamma^2 - \left(\Pi_0^2 - f^{-2}\omega^2\right)\left(1 - f^{-2}\omega^2\right)}}{\left(1 - f^{-2}\omega^2\right)}
```

Given the streamfunction we can compute ``u``, ``w``, ``v``, and ``b`` by
```math
\begin{aligned}
    u = -\psi_z &= mA\sin(kx + mz - \omega t) \\ 
    w = \psi_x &= -kA\sin(kx + mz - \omega t) \\ 
    v_t = -u(f + V_x) - wV_z \implies v &= \left[\frac{f\Gamma k}{\omega} - \frac{fm}{\omega}\right]A\cos(kx + mz - \omega t) \\ 
    b_t = -uB_x - wB_z \implies b &= \left[\frac{f^2\Pi_0^2 k}{\omega} - \frac{f^2\Gamma m}{\omega}\right]A\cos(kx + mz - \omega t)
\end{aligned}
```

``\psi_t`` is then computed from the thermal wind imbalance
```math
\psi_t = \left(\frac{\partial^2~}{\partial x^2} + \frac{\partial^2~}{\partial z^2}\right)^{-1}\left(b_x - fv_z\right)
```
