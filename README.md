---
author:
- Marin Lauber
title: 'Pseudo-Spectral Method for 2D Homogeneous Isotropic Decaying Turbulence'
---

*This paper describes the theory behind obtaining and solving the
vorticity-stream function formulation of the incompressible
Navier-Stokes equations. A basic `Fortran` (and `Python`) implementation
of the method and algorithm described herein are available at:
<https://github.com/marinlauber/PS_HIT>.*

Vorticity equation {#vorticity-equation .unnumbered}
------------------

The vorticity (or *Helmholtz's*) equation describes the transport of
vorticity and is obtained by applying the curl operator to the
*Navier-Stokes* equations
$$\nabla \times\left(  \frac{\partial{\bm u}}{\partial t} + ({\bm u}\cdot \nabla){\bm u} = -\frac{1}{\rho}\nabla p + \nu \nabla^2 {\bm u}\right),$$
where ${\bm u}$ is the velocity field and $p$ is the pressure field. The
non-linear convective term can be rearranged, making use of the property
$${\bm u}\times\left(\nabla\times{\bm u}\right) + \left({\bm u}\cdot\nabla\right){\bm u} = \nabla\left(\frac{{\bm u}^2}{2}\right) \rightarrow  \left({\bm u}\cdot\nabla\right){\bm u} = \nabla\left(\frac{{\bm u}^2}{2}\right) + \underbrace{\left(\nabla\times{\bm u}\right)}_{={\bm\omega}}\times{\bm u},$$
where we have made use of the definition of the vorticity
(${\bm\omega}=\nabla\times{\bm u}$). This leaves us with
$$\nabla \times\left(  \frac{\partial{\bm u}}{\partial t} + \nabla\left(\frac{{\bm u}^2}{2}\right) + {\bm \omega} \times{\bm u} = -\frac{1}{\rho}\nabla p + \nu \nabla^2 {\bm u}\right).$$
Distributing the curl operator and making use of the definition of the
vorticity, we have
$$\frac{\partial{\bm \omega}}{\partial t} + \nabla \times\nabla\left(\frac{{\bm u}^2}{2}\right) + \nabla \times\left({\bm \omega} \times{\bm u}\right) = -\nabla \times\left( \frac{1}{\rho}\nabla p\right) + \nu \nabla^2 {\bm \omega}.$$
The third term can be expanded using the *curl of the curl* vector
identity
$$\frac{\partial{\bm \omega}}{\partial t} + \nabla \times\nabla\left(\frac{{\bm u}^2}{2}\right) + {\bm u}\cdot\nabla{\bm\omega} - {\bm u}\nabla\cdot{\bm\omega} - {\bm\omega}\cdot\nabla{\bm u} + {\bm\omega}\nabla\cdot{\bm u} = -\nabla \times\left( \frac{1}{\rho}\nabla p\right) + \nu \nabla^2 {\bm \omega}.$$
Because we are dealing with an incompressible flow, the continuity
equation provides us with $\nabla\cdot{\bm u} =0$ and we also know that
$\nabla\cdot{\bm\omega}\equiv 0$, by definition[^1]. We are thus left
with only
$$\frac{\partial{\bm \omega}}{\partial t}  + {\bm u}\cdot\nabla{\bm\omega} - {\bm\omega}\cdot\nabla{\bm u}  = -\nabla \times\left( \frac{1}{\rho}\nabla p\right) + \nu \nabla^2 {\bm \omega},$$
in which the second term is also zero, this can be proved using the
definition of the curl of a tensor field
$$(\nabla\times{\bm S})\cdot{\bm a} = \nabla\times({\bm S}^\top \cdot{\bm a}) \qquad \forall {\bm a.}$$
For our case, this gives
$$\left(\nabla \times\nabla\left(\frac{{\bm u}^2}{2}\right)\right)\cdot{\bm a} = \nabla\times\left(\left(\frac{{\bm u}^2}{2}\right)^\top\cdot{\bm a}\right) = \nabla\times\left(\nabla\left(\frac{{\bm u}^2}{2}\cdot{\bm a}\right)\right) = \nabla\times\nabla\psi\equiv0.$$
As the curl of the gradient of a scalar field is always zero. Finally,
the pressure term can be expanded using the properties of the *curl*
operator
$$\nabla \times\left( \frac{1}{\rho}\nabla p\right) = \nabla\left( \frac{1}{\rho}\right) \times\nabla p + \left( \frac{1}{\rho}\right)\nabla \times\nabla p = \frac{1}{\rho^2}\nabla\rho\times\nabla p.$$
For a flow with constant entropy (incompressible, homentropic), the
pressure is solely a function of density, this means that the two vector
$\nabla\rho$ and $\nabla p$ are aligned and the pressure (more
precisely, the baroclinic term which account for the change in the
vorticity field due to the interaction of density and pressure surfaces
[@bellot]) term vanishes. The vorticity equation is therefore given by
$$\frac{\partial{\bm \omega}}{\partial t} + {\bm u}\cdot\nabla{\bm\omega} = {\bm\omega}\cdot\nabla{\bm u} + \nu \nabla^2 {\bm \omega}.$$
In *2D*, these equations can be simplified further by noting that the
component of the velocity field are
$${\bm u} = (u, v, 0) \qquad\text{and}\qquad \frac{\partial}{\partial z} =0.$$
The vorticity is, by definition, normal to the $x$-$y$ plane
(${\bm\omega}=(0, 0, \omega_z)$). This means that the term
${\bm\omega}\cdot\nabla{\bm u}$ is zero, as is easily demonstrated
$${\bm\omega}\cdot\nabla{\bm u} = \left(\underbrace{\omega_x}_{0}\frac{\partial}{\partial x} + \underbrace{\omega_y}_{0}\frac{\partial}{\partial y} + \omega_z\underbrace{\frac{\partial}{\partial z}}_{0}\right){\bm u}.$$
This term is the so-called vortex-stretching term. It is the mean by
with the energy is cascaded from the large scales to the smaller scales
following the $k^{-5/3}$ rule. In *2D* this means that the angular
velocity of a particle is conserved, and the length of a vortex tube
cannot change due to continuity [@mitnotes].\
*In *3D*, the vortex-stretching term has the form
$\omega_{j}\frac{\partial u_i}{\partial x_j}$, thus looking at, say the
second component, we have
$$\omega_{j}\frac{\partial u_i}{\partial x_j} = \omega_{1}\frac{\partial u_2}{\partial x_1} + \omega_{2}\frac{\partial u_2}{\partial x_2} + \omega_{3}\frac{\partial u_2}{\partial x_3},$$
the first and the third components can be seen as *vortex-turning* and
the second component is the vortex-stretching.*\
We also note that because the vorticity field has only one non-zero
component, the only equation for the transport of vorticity is given by
$$\label{equ:vorticity}
    \frac{\partial\omega}{\partial t} +u\frac{\partial\omega}{\partial x} + v\frac{\partial\omega}{\partial y} =  \frac{1}{Re} \left(\frac{\partial^2\omega}{\partial x^2}+ \frac{\partial^2\omega}{\partial y^2}\right).$$
Where $Re$ is the Reynold's number of the flow.

Stream function {#stream-function .unnumbered}
---------------

To remove the implicit dependency on the velocity field from the
vorticity equation, we introduce a vector-valued stream function
${\bm\psi}$. We will now show that for *2D* flows this is in fact a
scalar function. Using the definition of the stream function, we have
$${\bm u} = \nabla\times{\bm\psi}, \qquad u_i = \varepsilon_{ijk}\frac{\partial \psi_k}{\partial x_j},$$
where the only two non-zero components of the velocity field are
$$u_1 = \varepsilon_{123}\frac{\partial \psi_3}{\partial x_2}=\frac{\partial \psi_3}{\partial x_2}, \qquad u_2 = \varepsilon_{213}\frac{\partial \psi_3}{\partial x_1}=-\frac{\partial \psi_3}{\partial x_1}.$$
With $\varepsilon_{ijk}$, the Levi-Civita permutation symbol. Thus only
one component of the stream-function is non-zero, and it can be seen as
a scalar field $\psi$, with the following relationship to the velocity
field
$$u = \frac{\partial \psi}{\partial y}, \,\,\,\, v = -\frac{\partial \psi}{\partial x}.$$
We can substitute these relation into the vorticity equation
([\[equ:vorticity\]](#equ:vorticity){reference-type="ref"
reference="equ:vorticity"}), giving $$\label{equ:4}
    \frac{\partial\omega}{\partial t} +\frac{\partial \psi}{\partial y}\frac{\partial\omega}{\partial x} - \frac{\partial \psi}{\partial x}\frac{\partial\omega}{\partial y} =  \frac{1}{Re} \left(\frac{\partial^2\omega}{\partial x^2}+ \frac{\partial^2\omega}{\partial y^2}\right).$$
A very useful property of the stream function is that is existence
grantees the flow to be solenoidal (satisfies continuity). Indeed,
substitution of the definition of the stream function into the
continuity equation gives
$$\nabla\cdot{\bm u} = \nabla\cdot(\nabla\times{\bm\psi})\equiv 0,$$
which is zero by definition.\
Equation [\[equ:vorticity\]](#equ:vorticity){reference-type="ref"
reference="equ:vorticity"} still has 2 unknowns, $\psi$ and $\omega$,
fortunately we can build a Poisson equation for the vorticity by
substituting the velocity components, in terms of the stream function,
in the definition of the only non-zero component of the vorticity
$$\omega_{i}=\varepsilon_{ijk}\frac{\partial u_k}{\partial x_j} \qquad \rightarrow \qquad \omega =\omega_{z}=\frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}$$
giving $$\label{equ:stream}
\begin{split}
    \omega =-\frac{\partial }{\partial x}\frac{\partial \psi}{\partial x} - \frac{\partial }{\partial y}\frac{\partial \psi}{\partial y} &= -\frac{\partial^2 \psi}{\partial x^2} - \frac{\partial^2 \psi}{\partial x^2},\\
    -\omega&=\frac{\partial^2 \psi}{\partial x^2} + \frac{\partial^2 \psi}{\partial x^2}.
\end{split}$$
Equation [\[equ:vorticity\]](#equ:vorticity){reference-type="ref"
reference="equ:vorticity"} together with
equation [\[equ:stream\]](#equ:stream){reference-type="ref"
reference="equ:stream"} can be solved to simulate the transport of a
vorticity field.

Spectral collocation method {#spectral-collocation-method .unnumbered}
---------------------------

For periodic problems, such as the one of homogeneous isotropic decaying
turbulence in a periodic domain, collocation methods using trigonometric
polynomial expansions are unmatched in terms of accuracy and precision.
However, they tend to be an order of magnitude more expensive
(computationally) when compared to high-order finite difference schemes.
fortunately, when using trigonometric basis functions, the use of Fast
Fourier Transform (*FFT*'s) routines allows for a reduction in the
computational time. We describe here the implementation of such a method
to solve equation [\[equ:4\]](#equ:4){reference-type="ref"
reference="equ:4"} and
[\[equ:stream\]](#equ:stream){reference-type="ref"
reference="equ:stream"}[^2]\
We define the forward Fourier transform of a discrete function $u$ as
(transforming from Fourier to physical space), following [@article]
$$\label{equ:fft}
    u_{i,j} = \sum_{m=-\frac{N_x}{2}}^{\frac{N_x}{2}-1}\sum_{n=-\frac{M_y}{2}}^{\frac{M_y}{2}-1}\tilde{u}_{m,n}e^{{\bm i}\left( \frac{2\pi m}{L_x} x_i + \frac{2\pi n}{L_y}y_j \right)},$$
with its opposite transform, the backward Fourier transform
(transforming from physical to Fourier space)
$$\tilde{u}_{m,n} = \frac{1}{N_xN_y}\sum_{i=0}^{N_x-1}\sum_{j=0}^{N_y-1} u_{i,j}e^{{\bm -i}\left( \frac{2\pi m}{L_x} x_i + \frac{2\pi n}{L_y}y_j \right)}.$$
We note here that the transforms are normalised, a forward, followed by
a backward transform recover the initial discrete field. In this case we
have chosen to normalise the backward transform. The way numbers are
defined as
$$k_x = \frac{2\pi m}{L_x}, \qquad k_y = \frac{2\pi n}{L_y}.$$ The
natural ordering of Fourier coefficients is as shown in
equation [\[equ:fft\]](#equ:fft){reference-type="ref"
reference="equ:fft"}, however, *FFT* routines tend to organise them in a
different way, for a reason of algorithm efficiency and are ordered as
follows
$$k_x = \left[0, 1, \cdots, \frac{N_x}{2}-1, -\frac{N_x}{2}, \cdots, -2, -1 \right] \,\,\, \text{for even} \,\, N_x.$$
The $N_x$ discrete collocation (grid) points are uniformly space in our
$[L_x, L_y]$ domain as
$$x_i = \frac{iL_x}{N_x} ,\qquad y_i = \frac{iL_y}{N_y},$$ with
$i=0, 1, \cdots, N_x/N_y$. Because the problem is periodic,
$x_0 = x_{N_x}$, we do not include the last grid point in the
simulations, as collocation methods using periodic basis functions
automatically impose the periodicity on the problem.\
The high precision of the spectral collocation method comes from its
ability to threat spatial derivatives. By transforming a discrete
physical field to Fourier-space, what would otherwise be done using a
differential scheme (truncated to a certain accuracy), derivatives are
evaluated by multiplying the Fourier coefficients by the corresponding
complex wave number, for example, the $n^{th}$ derivative of a discrete
function $u$ is given by
$$\frac{\partial^{(n)} u_{i, j}}{\partial x^{(n)}} = \sum_{m=-\frac{N_x}{2}}^{\frac{N_x}{2}-1}\sum_{n=-\frac{M_y}{2}}^{\frac{M_y}{2}-1}\tilde{u}_{m,n}({\bm i}k_x)^{(n)}e^{{\bm i}\left( \frac{2\pi m}{L_x} x_i + \frac{2\pi n}{L_y}y_j \right)},$$
this is the power of spectral collocation methods. The only errors
associated with this operation are the interpolation and aliasing errors
[@spectral]. Applying these operators to the vorticity transport
equation result in a transport equation for the discrete Fourier
coefficients
$$\frac{\partial\tilde{\omega}_{m,n}}{\partial t} + ({\bm i}k_y\tilde{\psi}_{m,n}\circ{\bm i}k_x\tilde{\omega}_{m,n}) - ({\bm i}k_x\tilde{\psi}_{m,n}\circ{\bm i}k_y\tilde{\omega}_{m,n}) = \frac{1}{Re}\left[\left(-k_x^2 - k_y^2\right)\tilde{\omega}_{m,n}\right]$$
where the $\circ$ operator represents a convolution sum. Treating of
this term in Fourier space involves solving triad interactions, this is
very expansive, $O(N^2)$. In the pseudo-spectral approach, this term is
treated in physical space with the help of the convolution theorem and
*FFT's*, which reduces the cost to $O(45/4N\log_2 (3/2N))$. To prevent
aliasing of higher frequencies generated by the non-linear term,
de-aliasing must be done. Here we use $3/2$ rule to pad the wave number
during the transforms. Instead of using the inverse *FFT* with $N$
points, we use $M = 3/2N$ points. The new Fourier coefficients are
zero-padded, that isotropic $$\begin{split}
    \tilde{\omega}_{m, n} &= \tilde{\omega}_{m, n} \quad \text{if} \quad |m|\le N_x \quad |n|\le N_y,\\
    \tilde{\omega}_{m, n} &= 0 \quad \text{otherwise}.\\
\end{split}$$ Once we have obtained the two terms in the convolution sum
in physical space via the backward transform and using $M$ points, we
simply proceed to a multiplication. The result is then transformed back
in Fourier space using the forward transform, still with $M$ points. We
then discard all wave numbers such that $|m|>N_x$ and $|n|>N_y$. In
summary, the convolution sums are treated as
$$({\bm i}k_y\tilde{\psi}_{m,n}\circ{\bm i}k_x\tilde{\omega}_{m,n}) = \left[\mathcal{F}_{M}\left(\mathcal{F}^{-1}_{M}({\bm i}k_y\tilde{\psi}_{m,n})\mathcal{F}^{-1}_{M}({\bm i}k_x\tilde{\omega}_{m,n})\right)\right]_{N}\,,$$
where $\mathcal{F}_{M}$ represent a *Fast Fourier Transform* using $M$
points, but we are ultimately interested only in $N$ of them, and the
remaining are discarded. The Fourier coefficients of the stream function
are obtained by the explicit relationship between vorticity and
stream-function
$$-\tilde{\omega}_{m, n} = (-k_x^2 - k_y^2)\tilde{\psi}_{m,n} \qquad\rightarrow\qquad \tilde{\psi}_{m,n} = \frac{\tilde{\omega}_{m, n}}{(k_x^2 + k_y^2)}.$$
Which is easily obtained at every evaluation of the Navier-Stokes
operator.\
A pseudo-algorithm for a single iteration of the above procedure, or for
one step of a Runge-Kutta method, can be written as follows:

Get spectral stream-function:
$\qquad\qquad \tilde{\psi}_{m,n} = \tilde{\omega}_{m, n}/(k_x^2 + k_y^2)$
Set the un-defined stream function coefficient to zero
$\qquad\qquad \tilde{\psi}_{0, 0} = 0.0$ Build both convolution sums
with zero-padding:
$\qquad\qquad ({\bm i}k_y\tilde{\psi}_{m,n}\circ{\bm i}k_x\tilde{\omega}_{m,n}) = \left[\mathcal{F}_{M}\left(\mathcal{F}^{-1}_{M}({\bm i}k_y\tilde{\psi}_{m,n})\mathcal{F}^{-1}_{M}({\bm i}k_x\tilde{\omega}_{m,n})\right)\right]_{N}$
$\qquad\qquad ({\bm i}k_x\tilde{\psi}_{m,n}\circ{\bm i}k_y\tilde{\omega}_{m,n}) = \left[\mathcal{F}_{M}\left(\mathcal{F}^{-1}_{M}({\bm i}k_x\tilde{\psi}_{m,n})\mathcal{F}^{-1}_{M}({\bm i}k_y\tilde{\omega}_{m,n})\right)\right]_{N}$
Where all transforms use $M$ points, but we ultimately keep $N$ of them.
Construct the non-linear operator by summing both convolution sums:
$\qquad\qquad \mathcal{L}=({\bm i}k_y\tilde{\psi}_{m,n}\circ{\bm i}k_x\tilde{\omega}_{m,n}) + ({\bm i}k_x\tilde{\psi}_{m,n}\circ{\bm i}k_y\tilde{\omega}_{m,n})$
Build the dissipation term
$\qquad\qquad \mathcal{D} = \left[\left(-k_x^2 - k_y^2\right)\tilde{\omega}_{m,n}\right]$
Assemble the complete Navier-Stokes operator
$\qquad\qquad \mathcal{N} = -\mathcal{L} + \mathcal{D}$

This has to be solved at each call of the $\mathcal{N}$ operator.

Time integration {#time-integration .unnumbered}
----------------

Using spectral operators for the spatial derivatives ensures that the
discretisation error of spatial terms is kept to a minimum. In order to
keep the global discretisation errors low, high-order temporal
integrator are required. We implement both a third and a fourth-order
Runge-Kutta method. The third-order Runge-Kutta method id defined
as [@spectral] $$\begin{split}
    & K_1 = \mathcal{N}(\tilde{\omega}_n),\\
    & K_2 = \mathcal{N}(\tilde{\omega}_n + \frac{1}{2}\Delta t K_1),\\
    & K_3 = \mathcal{N}(\tilde{\omega}_n + \frac{3}{4}\Delta t K_2),\\
    & \tilde{\omega}_{n+1} = \tilde{\omega}_n + \frac{1}{9}\Delta t\left(2K_1 + 3K_2 + 4K_3\right),\\
\end{split}$$ where $\mathcal{N}(\tilde{\omega}_n)$ represents all the
spatial operators of the transport equation of the Fourier coefficients.
Similarly, the fourth-order Runge-Kutta is defined as $$\begin{split}
    & K_1 = \mathcal{N}(\tilde{\omega}_n),\\
    & K_2 = \mathcal{N}(\tilde{\omega}_n + \frac{1}{2}\Delta t K_1),\\
    & K_3 = \mathcal{N}(\tilde{\omega}_n + \frac{1}{2}\Delta t K_2),\\
    & K_4 = \mathcal{N}(\tilde{\omega}_n + \Delta t K_3),\\
    & \tilde{\omega}_{n+1} = \tilde{\omega}_n + \frac{1}{6}\Delta t\left(K_1 + 2K_2 + 2K_3 + K_4\right).\\
\end{split}$$ These schemes can be applied to physical and spectral
quantities. To avoid using two additional *FFT*'s per iteration,
time-integration is made in spectral space.

Validation {#validation .unnumbered}
----------

We validate the numerical method described above using the
*Taylor-Green* vortex, a known solution to the *Navier-Stokes* equation.
Expressed in its *2D* vorticity formulation, it reads
$$\omega^e(x, y, t) = 2\kappa\cos(\kappa x)\cos(\kappa y)e^{-2\kappa^2t/Re},$$
where $\kappa=4$ is the number of vortices in each direction and $Re$ is
the Reynolds number of the flow and is set to 1.0 for those simulations
(dissipation is important in this flow). Clearly, the flow that has this
as a solution is purely dissipative, we should not expect the non-linear
(or convolution sum in our case) to be important[^3]. This means that
the vorticity is invariant under convection. We use the $L_2$ norm of
the error between the exact solution and the simulation to quantify the
error
$$L_2 = \sqrt{\sum_{i=1}^{N_x}\sum_{j=1}^{N_y}|\omega^{e}_{i, j} - \omega_{i, j}|^2}.$$
Four different spatial resolutions are used to investigate the
convergence in the $l_2$-norm of the error. The results are presented on
the left of figure [\[fig:l2\]](#fig:l2){reference-type="ref"
reference="fig:l2"}. Convergence close to machine accuracy is obtained
for all spatial resolutions. The increase of error is due to round-off
errors accumulating in the solution.

![image](valid/L2norm.pdf){width="\textwidth"}

 

![image](valid/L2normRK.pdf){width="\textwidth"}

The temporal convergence of the two time integrator used is presented on
the right of figure [\[fig:l2\]](#fig:l2){reference-type="ref"
reference="fig:l2"}. The actual convergence rate of both method is found
to be in excellent agreement with the expected rate.

![image](valid/vorticity.pdf){width="60%"}

Decaying Turbulence {#decaying-turbulence .unnumbered}
-------------------

To simulate the decay of homogeneous isotropic turbulence, or more
specifically, the decay of the turbulent kinetic energy of the flow, we
must first generate an appropriate vorticity field. This vorticity field
must respect a given energy spectrum and satisfy continuity. We base our
initial energy spectrum on, following [@article]
$$E(k) = \frac{a_s}{2}\frac{1}{k_p}\left(\frac{k}{k_p}\right)^{2s+1}e^{-\left(s+\frac{1}{2}\right)\left(\frac{k}{k_p}\right)^2},$$
where $k_p$ is the wave number at which the peak turbulent energy is
located, $s$ is a shape factor for the spectrum and $a_s$ is a
normalisation parameter, defined as
$$a_s = \frac{(2s+1)^{s+1}}{2^ss!}.$$ This one dimensional spectrum is
the used to generate the *2D* vorticity field by a random distribution
of phase to the Fourier coefficients in the $k_x$-$k_y$ plane
$$\tilde{\omega}({\bm k}) = \sqrt{\frac{k}{\pi}E(k)}e^{{\bm i}\zeta({\bm k})},$$
where the random phase function, $\zeta({\bm k}, t)$ is defined as
follows for the first, second, third and fourth quadrants in the
$k_x$-$k_y$ plane, respectively $$\begin{split}
    \zeta(k_x, k_y) & = \xi(k_x, k_y) + \eta(k_x, k_y) , \qquad \xi({\bm k}),  \zeta({\bm k}) \in \mathcal{N}[0, 2\pi],\\
    \zeta(-k_x, k_y) & = -\xi(k_x, k_y) + \eta(k_x, k_y), \\
    \zeta(-k_x, -k_y) & = -\xi(k_x, k_y) - \eta(k_x, k_y), \\
    \zeta(k_x, -k_y) & = \xi(k_x, k_y) - \eta(k_x, k_y), \\
\end{split}$$ this ensures (approximately) that the Hermitan
symmetry[^4] of the Fourier coefficients of a real-valued signal is
respected when constructing the spectral vorticity distribution.

![image](vorticity.pdf){width="\textwidth"}

 

![image](vorticity1.pdf){width="\textwidth"}

Energy Spectrum {#energy-spectrum .unnumbered}
---------------

For homogeneous turbulence, we can define the energy spectrum to be a
function of the stream-function
$$E(k, t) = \sum_{k-\Delta k \le |{\bm k}| \le k+\Delta k} \frac{1}{2} k^2 |\tilde{\psi}({\bm k}, t) |^2,$$
where $\tilde{\psi}$ is readily obtained from the known vorticity field
$\tilde{\omega}$. Each contribution to the final energy spectrum consist
of the sum of all the energy contribution at each wave number falling
within a prescribed bin size, or spectral band, $\Delta k$. The actual
value of the spectral band is determined by comparing theoretical and
computed spectrum, see right of
figure [\[fig:spectrum\]](#fig:spectrum){reference-type="ref"
reference="fig:spectrum"}. The *KBL* theory predicts that, in the limit
of infinite Reynolds number, the inertial range approaches a $k^{-3}$
scaling. This cannot be observed in
figure [\[fig:spectrum\]](#fig:spectrum){reference-type="ref"
reference="fig:spectrum"} as this is the initial spectrum, and thus the
*TKE* has not yet been cascaded down the inertial range.

![image](spectrum.pdf){width="60%"}

Running the Code {#running-the-code .unnumbered}
----------------

The code provided requires the use of the `FFTPack 5.1` fast Fourier
transform library, available at
<https://www2.cisl.ucar.edu/resources/legacy/fft5>. Once the library has
been installed, the `Makefile` procided with the code will create the
executable `main`. The `main.f90` file contains the main loop of the
solver, and calls different subroutines stored in the other `.f90` file.
This file is also used to specify the simulation size, time, initial
conditions, file IO, etc. The following table describe the subroutines
contained in all the file, and their use.

  Fortran file        Subroutines                   function
  ------------------- ----------------------------- ----------------------------------------------
  `main.f90`          None                          Main loop of program
  `fields.f90`        `taylorgreen(...)`            Generate *Taylor Green* vortex
                      `vorticity_field(...)`        Generate decaying vorticity field
                      `compute_spectrum(...)`       Compute *TKE* spectrum
  `spatial.f90`       `stream_function(...)`        Compute $\tilde{\psi}$ from $\tilde{\omega}$
                      `discrete_NS_operator(...)`   Build discrete Navier-Stokes operator
                      `convolve(...)`               Concolve two fields with 3/2 padding
  `runge_kutta.f90`   `rk3_step(...)`               One step of Runge-Kutta 3
                      `rk4_step(...)`               One step of Runge-Kutta 4
  `fft.f90`           `ifft2(...)`                  Wrapper for backward FFT
                      `fft2(...)`                   Wrapper for forward FFT
  `utils.f90`         `hello(...)`                  Print welcome message
  `constants.f90`     None                          Stores constants ($\pi$, $e$, ${\bm i}$)
  `types.f90`         None                          Define variable types (dp)

  : List of `.f90` file and their content[]{label=""}

Final Comment {#final-comment .unnumbered}
-------------

Copyright (c) 2019 Marin Lauber\
Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
\"Software\"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:\
The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.\
THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.\

1

C. Bailly, G. Comte-Bellot. *Turbulence*. Springer Series in:
Experimental Fluid Mechanics. 2015.

MIT. *Marine Hydrodynamics Lecture 9*. MIT Marine Hydrodynamics lectures
notes. 2005.

O. San, A. E. Staples. *High-order methods for decaying two-dimensional
homogeneous isotropic turbulence*. Computer & Fluids 63. p: 105--127.
2012.

C. Canuto, M. Y. Hussaini, A.Quarteroni, T. A. Zang. *Spectral Methods
in Fluid Dynamics*. Springer Series in: Computational Physics. 1988.

[^1]: The divergence of the curl of a vector field is always zero.

[^2]: Spectral collocation methods solve for the strong form of the
    governing equations, as opposed to Galerkin methods.

[^3]: in fact we can solve it by setting $\mathcal{L}=0$ is the above
    algorithm and still recover the right answer.

[^4]: negative wave number Fourier coefficients are the complex
    conjugate of positive ones:
    $\tilde{\omega}(-k) = \overline{\tilde{\omega}(k)}$.
