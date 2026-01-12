# Flow of a generalized Newtonian fluid

In this chapter we consider the flow of incompressible generalized Newtonian fluids where inertia can be neglected. These flows can be modelled by the Stokes equations, Equations {eq}`eq49-chap1`, where
the constant viscosity $\mu$ is replaced by a flow-dependent viscosity. Typical applications where a generalized Newtonian fluid model can be useful are: flows of polymeric melts and solutions with a dominant shear component, such as pipes, channels with only mildly varying cross section.

The finite element discretization will be similar to the Stokes equations in {numref}`Chap3`. However, the resulting discretized equations are nonlinear and require iteration to obtain a solution. Two iteration methods will be discussed: Picard and Newton-Raphson.

## Viscous fluids

Viscous fluids are defined as fluids that have an instantaneous stress response with respect to the velocity gradient. Starting the flow gives an immediate stress response and the stress is removed immediately without delay if the flow is stopped. The most general viscous fluid model is 

$$
  \ten \tau = \ten{\mathcal{F}}(\ten L)
$$ (eq1-chap4)

where $\ten L=(\nao \vek u)^T$ and $\ten{\mathcal{F}}$ is a function with $\ten{\mathcal{F}}(\ten 0)=\ten 0$. The principle of super-imposed rigid body motions, also called principle of material objectivity, requires that the function $\ten{\mathcal{F}}$ must fulfill

$$
  \ten Q\cdot\ten{\mathcal{F}}(\ten L)\cdot\ten Q^T = \ten{\mathcal{F}}(\ten\Omega+\ten Q\cdot\ten L\cdot\ten Q^T)
$$ (eq2-chap4)

for any rotation tensor $\ten Q$ (property: $\ten Q\cdot\ten Q^T=\ten I$, $\det\ten Q=1$) and rotation rate tensor $\ten \Omega$ (property: $\ten \Omega^T=-\ten \Omega$). For the special case $\ten Q=\ten I$ and $\ten\Omega=-\ten W$, we see that

$$
 \ten{\mathcal{F}}(\ten L) =\ten{\mathcal{F}}(\ten D)
$$ (eq3-chap4)

and therefore

$$
  \ten \tau = \ten{\mathcal{F}}(\ten D)
$$ (eq4-chap4)

where $\ten{\mathcal{F}}$ must fulfill:

$$
  \ten Q\cdot\ten{\mathcal{F}}(\ten D)\cdot\ten Q^T = \ten{\mathcal{F}}(\ten Q\cdot\ten D\cdot\ten Q^T)
$$ (eq5-chap4)

From this equation we can conclude, that $\ten{\mathcal{F}}(\ten D)$ must be an isotropic tensor function. The most general form of the viscous constitutive model is therefore given by

$$
  \ten \tau = \mathcal{F}_0\ten I+\mathcal{F}_1\ten D+\mathcal{F}_2\ten D^2
$$ (eq6-chap4)

where $\mathcal{F}_0$, $\mathcal{F}_1$ and $\mathcal{F}_2$ are functions of the three invariants only. This is called a *Reiner-Rivlin fluid*. 
````{exercise}
:label: ex:4.1
Show that for the Reiner-Rivlin fluid with $\mathcal{F}_2\neq 0$, in simple shear flow the *first normal stress difference* $\tau_{xx}-\tau_{yy}=0$ and the *second normal stress difference* $\tau_{yy}-\tau_{zz}\neq0$. Here $x$ is the flow direction, $y$ the gradient direction and $z$ the vorticity direction.

````

The Reiner-Rivlin fluid given by Equation {eq}`eq6-chap4` is only relevant for the history of rheology, since no real fluids have been found that obey this model with $\mathcal{F}_2\neq 0$. Therefore, for all practical purposes we can assume $\mathcal{F}_2=0$. Also, the isotropic term (the first term in the right-hand side of Equation {eq}`eq6-chap4`)
can be included in the pressure term of the Cauchy stress tensor ($-p\ten I$) and the most general viscous model for all practical purposes becomes

$$
  \ten \tau = 2\eta(\bar I_D,\bar{I\!I}_D, \bar{I\!I\!I}_D)\ten D
$$ (eq7-chap4)

where we have substituted $\mathcal{F}_1=2\eta$ and 
we chose the moments of $\ten D$ as the invariants (see Eq.~\eqref{eq:moments} for the definition of the moments of a tensor). 
Since we have assumed incompressible flow $\bar I_D=\nao\cdot\vek u=0$, and we have only a dependence on two invariants left:

$$
  \ten \tau = 2\eta(\bar{I\!I}_D, \bar{I\!I\!I}_D)\ten D
$$ (eq8-chap4)

Assuming the (viscosity) function $\eta$ is determined from shear data only, 
there is no possibility to obtain the dependence on the third invariant, since
$\bar{I\!I\!I}_D=0$ for shear flow. In fact, $\bar{I\!I\!I}_D=0$ for any 2D flow. 

````{exercise}
:label: ex:4.2
Consider a uni-axial elongation flow, where 
$\ten D=\dot\epsilon\vek e_x\vek e_x-\frac12\dot\epsilon\vek e_y\vek e_y-\frac12\dot\epsilon\vek e_z\vek e_z$. Show that $\bar{I\!I}_D=\frac32\dot\epsilon^2$ and $\bar{I\!I\!I}_D=\frac34\dot\epsilon^3$. Argue that 
$\dot\epsilon_\text{e}=2\bar{I\!I\!I}_D/\bar{I\!I}_D$ can be used as an *effective elongational rate* in general 3D incompressible flows. Show that $\dot\epsilon_\text{e}=0$ in 2D incompressible flows. 

````

Assuming that the function $\eta$ is independent from $\bar{I\!I\!I}_D$, we finally arrive at
the *generalized Newtonian* fluid model: 

$$
 \ten\tau =2\eta(\dot\gamma_\text{e})\ten D = \eta(\dot\gamma_\text{e}) (\nao\vek u+(\nao\vek u)^T) 
$$ (eq9-chap4)

where $\eta(\dot\gamma)$ is the viscosity function[^1] and $\dot\gamma_\text{e}$ is the
*equivalent shear rate*, which is a function of $\ten D$ and defined as

$$
  \dot\gamma_\text{e} = \sqrt{2\bar{I\!I}_D}=\sqrt{2\ten D:\ten D}
$$ (eq10-chap4)


````{admonition} Remark 4.1
:class: note

Although it is possible to use data measured in a simple shear flow directly for defining the function  $\eta(\dot\gamma)$, it is quite common to use “models“ to fit the actual shear data. Popular models are 

$$
\begin{alignat*}{3}
\eta(\dot\gamma) &~= \dfrac{m}{\dot\gamma^{1-n}} & \quad&\text{Power-law}\\
\eta(\dot\gamma) &= \eta_\infty + 
\dfrac{\eta_0-\eta_\infty}{\big(1+(\lambda\dot\gamma)^2\big)^{(1-n)/2}} &
        \quad&\text{Carreau}\\
\eta(\dot\gamma) &= \eta_\infty + 
\dfrac{\eta_0-\eta_\infty}{\big(1+(\lambda\dot\gamma)^a\big)^{(1-n)/a}} &
        \quad&\text{Carreau-Yasuda}
\end{alignat*}
$$ (eq11-chap4)

where $m$, $n$, $\eta_0$, $\eta_\infty$, $\lambda$ and $a$ are material constants. In {numref}`fig1-chap4` $\eta(\dot\gamma)$ for various models is plotted. The typical slope of the straight line in $\eta(\dot\gamma)$ over several decades is $n-1$.

```{figure-md} fig1-chap4

<img src="media/chap4/fig1-chap4.png"  width="600px">

Plot of $\eta(\dot\gamma)$  for various generalized Newtonian models ($n=0.3$).
```
````


````{exercise}
:label: ex:4.3

Assume $\eta_\infty=0$.
Show, that for $\lambda\dot\gamma \gg 1$ the Carreau and Carreau-Yasuda model can be approximated by a power law model having a value of $m=\eta_0/\lambda^{1-n}$.

````


````{exercise}
:label: ex:4.4

For $\dot\gamma\rightarrow0$ the viscosity of the power law model is unbounded. We want to limit the viscosity to some given value $\eta_\text{max}$. Show that by using a Carreau 
model with $\eta_0=\eta_\text{max}$, $\eta_\infty=0$ and $\lambda=(\eta_\text{max}/m)^{1/(1-n)}$, we can achieve just that and still approximate the original power-law model for $\lambda\dot\gamma \gg 1$.


````


````{admonition} Remark 4.2
:class: note

 The shear stress in a shear flow for a generalized Newtonian fluid is

$$
 \tau(\dot\gamma)=\eta(\dot\gamma)\dot\gamma
$$ (eq12-chap4)

In these lecture notes, we will only consider fluids where the shear stress curve is monotonically increasing with shear rate, i.e. 

$$
 \deriv{\tau}{\dot\gamma}=\deriv{\eta}{\dot\gamma}\dot\gamma+\eta=\eta(\deriv{\ln\eta}{\ln\dot\gamma}+1) > 0, \quad\text{for all }\dot\gamma > 0
$$ (eq13-chap4)

In {numref}`fig2-chap4` $\tau(\dot\gamma)$ for various models is plotted. The typical slope of the straight line in $\tau(\dot\gamma)$ over several decades is $n$. Therefore, the monotonicity requirement Equation {eq}`eq13-chap4` is being fulfilled only if $n>0$ is chosen for these models.

```{figure-md} fig2-chap4

<img src="media/chap4/fig2-chap4.png"  width="600px">

Plot of $\tau(\dot\gamma)$  for various generalized Newtonian models ($n=0.3$).
```

````

```{admonition} Remark 4.3
:class: note

If the slope of stress curve $\tau(\dot\gamma)$ becomes negative, i.e. $\lderiv{\tau}{\dot\gamma}<0$, the flow problem given by Equations {eq}`eq18-chap4` becomes of a different nature (hyperbolic instead of elliptic) in the regions where this happens {cite}`Regirer1968`. The nature of the solution, if it exists at all for given boundary conditions, also becomes different with possible extreme localization of the shear rate. The study of these flows is beyond the scope of these lecture notes. 

```

```{admonition} Remark 4.4
:class: note

The stress curve for a power law is given by $\tau(\dot\gamma)=m\dot\gamma^n$, i.e. the slope $\lderiv{\tau}{\dot\gamma}=nm \dot\gamma^{n-1}=n\eta(\dot\gamma)\rightarrow\infty$ for $\dot\gamma\rightarrow0$ and $0<n<1$. The unbounded stress slope (and viscosity $\eta$) for small shear rates can lead to numerical convergence problems (see {numref}`Chap4.5.1`).

```

## Problem formulation

The momentum and mass balance are given by Equations {eq}`eq42-chap1` and {eq}`eq46-chap1` read

$$
\begin{alignat*}{2}
 -\nao\cdot\ten\sigma^T&=\vek f&&\text{in }\Omega\\
   \nao\cdot\vek u&=0&\qquad&\text{in }\Omega
\end{alignat*}
$$ (eq14-chap4)

The Cauchy stress tensor $\ten\sigma$ is split into a pressure part and the extra-stress tensor $\ten\tau$, according to Equation {eq}`eq44-chap1`: 

$$
 \ten\sigma =-p\ten I+\ten \tau
$$ (eq15-chap4)

For a generalized Newtonian fluid the expression for $\ten\tau$ is given by the constitutive Equation {eq}`eq9-chap4`, instead of Equation {eq}`eq47-chap1` for a Newtonian one:

$$
 \ten\tau =2\eta(\dot\gamma_\text{e})\ten D = \eta(\dot\gamma_\text{e}) (\nao\vek u+(\nao\vek u)^T)
$$ (eq16-chap4)

where $\eta(\dot\gamma)$ is the viscosity function and $\dot\gamma_\text{e}$ is the equivalent shear rate Equation {eq}`eq10-chap4`:

$$
  \dot\gamma_\text{e} = \sqrt{2\bar{I\!I}_D}=\sqrt{2\ten D:\ten D}
$$ (eq17-chap4)

The equations for a generalized Newtonian fluid now become:

$$
\begin{alignat*}{2}
 -\nao\cdot(2\eta(\dot\gamma_\text{e}) \ten D)+\nao p &=\vek f&&\text{in }\Omega \\
   \nao\cdot\vek u&=0&\qquad&\text{in }\Omega
\end{alignat*}
$$ (eq18-chap4)

The boundary conditions (for now) are assumed to be

$$
\begin{alignat*}{2}
    \vek u &= \vek u_\text{D}&&\qquad\text{ on }\Gamma_D\\
    \ten\sigma\cdot\vek n=-p\vek n+2\eta(\dot\gamma_\text{e}) \ten D\cdot\vek n&=\vek t_{\text{N}} &&\qquad\text{ on }\Gamma_{\text{N}}
\end{alignat*}
$$ (eq19-chap4)

where the boundary $\Gamma=\Gamma_\text{D}\cup\Gamma_{\text{N}}$ and has been split into a *Dirichlet* and a *Neumann* part (like in the Stokes equations, see {numref}`fig1-chap3`). The Neumann boundary condition is equivalent to an imposed traction of $\vek t_{\text{N}}$.

## Weak formulation

The generalized Newtonian fluid model is basically constructed by replacing the constant viscosity coefficient $\mu$ with a flow dependent viscosity function $\eta(\dot\gamma_\text{e})$. Therefore, deriving the weak form is the same as for the Stokes equations in {numref}`Chap3.2`, except for the replacement of $\mu$ with $\eta(\dot\gamma_\text{e})$. After doing that, the following *weak form of the generalized Newtonian fluid equations* is found: find $\vek u\in \vek U$, $p\in P$ such that

$$
\begin{align*}
 \int_\Omega \eta(\dot\gamma_\text{e})(\nao\vek v)^T:(\nao\vek u+(\nao\vek u)^T)\,d\Omega-\int_\Omega\nao\cdot\vek vp\,d\Omega
        &=\int_{\Gamma_\text{N}} \vek v\cdot\vek t_\text{N}\,d\Gamma +\int_\Omega \vek v\cdot\vek f\,d\Omega\\
  \int_\Omega q\nao\cdot\vek u\,d\Omega &=0
\end{align*}
$$ (eq20-chap4)

for all $\vek v\in \vek V$ and $q\in Q$. In this equation $\dot\gamma_\text{e}$ is a function of $\ten D$ according to Equation {eq}`eq10-chap4`.

## The Galerkin method

Also for the discretization with the Galerkin method we can follow the same approach as in {numref}`Chap3.3` and replace the constant viscosity coefficient $\mu$ with a flow dependent viscosity function $\eta(\dot\gamma_\text{e})$.

The discretized equations for the flow of a generalized Newtonian fluid can be obtained by replacing the general spaces with the finite dimensional approximation spaces in the weak form Equations {eq20-chap4}: 
find $\vek u_h\in \vek U_h$, $p_h\in P_h$ such that

$$
\begin{align}
 \int_\Omega \eta(\dot\gamma_\text{e}^h)(\nao\vek v_h)^T:(\nao\vek u_h+(\nao\vek u_h)^T)\,d\Omega-\int_\Omega\nao\cdot\vek v_hp_h\,d\Omega
        &=\int_{\Gamma_\text{N}} \vek v_h\cdot\vek t_\text{N}\,d\Gamma +\int_\Omega \vek v_h\cdot\vek f\,d\Omega\label{eq:mom_genweak2}\\
  \int_\Omega q_h\nao\cdot\vek u_h\,d\Omega &=0\label{eq:contgenweak2}
\end{align}
$$ (eq21-chap4)

for all $\vek v_h\in \vek V_h$ and $q_h\in Q_h$. In this equation $\dot\gamma_\text{e}^h$ is given by

$$
  \dot\gamma_\text{e}^h =\sqrt{\big(\nao\vek u_h+(\nao\vek u_h)^T):(\nao\vek u_h+(\nao\vek u_h)^T\big)/2}
$$ (eq22-chap4)

The approximation spaces for velocity and pressure can be taken the same as for the Stokes equations, at least if the monotonicity requirement Equation {eq}`eq13-chap4` is being fulfilled. Using Equations {eq}`eq15-chap3`-{eq}`eq16-chap3` together with Equations {eq}`eq20-chap3` the following set of discretized equations can be derived:

$$
\begin{align}
    \vek g_{k}(\vek u_h) + \sum_{m=1}^{N_p} \vek b_{mk}p_m
  &=\vek f_k, \qquad k=1,\dots N_{u}\\
  \sum_{m=1}^{N_u}\vek b_{km}\cdot \vek u_m&=0\label{eq:cont_gen_disc0}, \qquad k=1,\dots N_{p}
\end{align}
$$ (eq23-chap4)

where $\vek b_{km}$ and $\vek f_k$ are the same as for Stokes (Equations {eq}`eq22-chap3` (two last equations)) and the vectors $\vek g_{k}(\vek u_h)$ are defined as

$$
 \vek g_{k}(\vek u_h) = \int_\Omega \eta(\dot\gamma_\text{e}^h)
  \nao\phi_k\cdot\big(\nao\vek u_h+(\nao\vek u_h)^T\big)\,d\Omega, \qquad
  k=1,\dots N_{u}
$$ (eq24-chap4)

Note, that because $\dot\gamma_\text{e}^h$ depends on $\vek u_h$, we now have obtained a *nonlinear set of equations*.



````{exercise}
:label: ex:4.5

Derive the expression for $\vek g_{k}(\vek u_h)$ in Equation {eq}`eq24-chap4`.
````

We can now define the column vector $\col g(\col u)$ as follows

$$
          \col{g}(\col{u})=  
  \begin{pmatrix}
    \col{g_{1}} \\
    \col{g_{2}} \\
      \vdots  \\
    \col{g_{N_u}} 
  \end{pmatrix}\label{eq:gcol}
$$ (eq25-chap4)

with $\col{g_{k}}$ the column of size $d$ of the components of $\vek g_{k}(\vek u_h)$ with respect to a Cartesian coordinate system given by the unit vectors $\vek e_i$, $i=1,\dots,d$. The column vector $\col u$ is the vector of unknown velocities defined indentically to the Stokes equations (see Equation {eq}`eq24-chap3` (second matrix)). The system given by Equations {eq}`eq23-chap4` can now be written in a matrix form as follows

$$
\begin{align*}
\col{g}(\col{u}) +\mat{B^T}\col{p} &= \col{f}\\
   \mat{B}\col{u}  &= \col{0} 
\end{align*}
$$ (eq26-chap4)

The column vectors $\col{p}$, $\col{f}$ and matrix $\mat{B}$ are the same as for the Stokes equations (see Equations {eq}`eq24-chap3`).

For later use, we write $\vek g_{k}(\vek u_h)$ as follows

$$
    \vek g_{k}(\vek u_h)=\sum_{m=1}^{N_u} \ten K_{km}(\vek u_h)\cdot \vek u_m, \qquad k=1,\dots N_{u}
$$ (eq27-chap4)

where the tensors $\ten K_{km}(\vek u_h)$ are defined as

$$
 \ten K_{km}(\vek u_h) = \int_\Omega \eta(\dot\gamma_\text{e}^h)
  [(\nao\phi_k\cdot\nao\phi_m)\ten I+\nao\phi_m\nao\phi_k]\,d\Omega, \qquad
  k,m=1,\dots N_{u}
$$ (eq28-chap4)

The column vector $\col g(\col u)$ can then be written as

$$
  \col{g}(\col{u})=\mat{K}(\col{u})\col{u}
$$ (eq29-chap4)

where $\mat K(\col u)$ has been defined as

$$
   \mat{K}(\col{u})=  
  \begin{pmatrix}
    \mat{K}_{11} & \mat{K}_{12} & \cdots & \mat{K}_{1N_u} \\
    \mat{K}_{21} & \mat{K}_{22} & \cdots & \mat{K}_{2N_u} \\
      \vdots & \vdots & \ddots & \vdots \\
    \mat{K}_{N_u1} & \mat{K}_{N_u2} & \cdots & \mat{K}_{N_uN_u}
  \end{pmatrix}
$$ (eq30-chap4)

with $\mat{K_{km}}$ the $d\times d$ matrices of the components of $\ten K_{km}(\vek u_h)$ with respect to a Cartesian coordinate system given by the unit vectors $\vek e_i$, $i=1,\dots,d$. The system Equations {eq}`eq26-chap4` can now be written as 

$$
\begin{align*}
\mat{K}(\col{u})\col{u} +\mat{B^T}\col{p} &= \col{f} \\
   \mat{B}\col{u}  &= \col{0}
\end{align*}
$$ (eq31-chap4)

or in full matrix form

$$
\begin{pmatrix}
  \mat{K}(\col{u}) & \mat{B^T}\\
  \mat{B} & \mat{0}
\end{pmatrix}
\begin{pmatrix}
 \col{u} \\
 \col{p}
\end{pmatrix}=
\begin{pmatrix}
 \col{f} \\
 \col{0}
\end{pmatrix}
$$ (eq32-chap4)

Note, that for a constant viscosity function $\eta(\dot\gamma)=\mu$, we have $\mat K(\col u)=\mat A$ and the discretized Stokes system Equation {eq}`eq27-chap3` is obtained.

The discretized equations for the flow of generalized Newtonian fluid model as obtained in this section are nonlinear in the velocity unknowns $\col{u}$. Therefore the solution cannot be obtained in `one go' and iteration is needed. In the next section we will introduce the Picard iteration scheme for that.


## Picard iteration
(Chap4.5.1)=
### Theory

The *Picard iteration scheme* is defined by a starting vector $\col{u^{(0)}}$ for the velocity  and a sequence of solutions $\{\col{u^{(j)}},\col{p^{(j)}\}}$, $j=1,2,3,\dots$, given by

$$
\begin{align*}
\mat{K}(\col{u^{(j-1)}})\col{u^{(j)}} +\mat{B^T}\col{p^{(j)}} &= \col{f}\\
   \mat{B}\col{u^{(j)}} &= \col{0}
\end{align*}
$$ (eq33-chap4)

Note, that we do not need a starting value for the pressure. 
Basically we solve a Stokes system every iteration with a viscosity obtained using the solution from the previous iteration. Since, it is just a Stokes solve with an adaptive viscosity it is straightforward to implement if you already have a Stokes solver.

```{admonition} Remark 4.5
:class: note

The starting vector $\col u^{(0)}$ should be chosen as close as possible to the actual solution and 
preferably already fulfilling the Dirichlet boundary conditions. For that, we use the Stokes solution

$$
\begin{align}
\mat{A}\col{u^{(0)}} +\mat{B^T}\col{p^{(0)}} &= \col{f} \\
\mat{B}\col{u^{(0)}} &= \col{0}
\end{align}
$$ (eq34-chap4)

In case there are Neumann boundary conditions, the viscosity $\mu$ should be chosen to represent a good estimate of the actual viscosity in the flow. The Stokes pressure solution $\col{p^{(0)}}$ is discarded.


```

(Chap4.5.2)=
### Convergence of the Picard iteration: lid-driven cavity flow

In order to study the convergence of the Picard iteration, we consider a flow on a square domain $\Omega=[-1,1]\times[-1,1]$. We impose the following Dirichlet conditions on the boundary of the domain:

$$
\begin{equation}
   \vek u_\text{D}=
         \begin{cases}
            \big(1-[\frac12-\frac12\cos(\pi x)]^{10}\big)\vek e_x & \text{upper boundary (`the lid', $y=1$)} \\
            \vek 0 & \text{otherwise}
         \end{cases}
\end{equation}
$$ (eq35-chap4)

The pressure level is imposed by imposing $p$ to be zero in the lower left corner. 

````{admonition} Remark 4.6
:class: note


The imposed “lid velocity” is plotted in {numref}`fig3-chap4`. At the ends $x=-1$ and $x=1$ both $u_x=0$ and $\plderiv{u_x}{x}=0$ to avoid singular behaviour for pressure and velocity gradients in the two upper corners.

```{figure-md} fig3-chap4

<img src="media/chap4/fig3-chap4.png"  width="700px">

Plot of the velocity in $x$-direction on the upper boundary ($y=1$).
```

````

We solve this problem using an $N\times N$ mesh of Taylor-Hood $Q_2$-$Q_1$ elements (see {numref}`fig13-chap3`). In {numref}`fig4-chap4` - {numref}`fig6-chap4` we show the Stokes solution of this problem: streamlines and $\dot\gamma_\text{e}$ based on $\col u^{(0)}$ and pressure contours[^2] of $p^{(0)}$ using $\mu=1$ and $N=40$. The imposed motion of the upper boundary generates a vortex with a strength (flow rate) of 0.199. Most of the “shearing” of the flow occurs near the lid with a maximum shear rate around 10.2, as shown by the contours of $\dot\gamma_\text{e}$. The pressure is symmetric up to machine precision as shown by the contour halfway between the max and min value of the solution.


```{figure-md} fig4-chap4

<img src="media/chap4/fig4-chap4.png"  width="600px">

Streamfunction contours of the Stokes solution $\vek u^{(0)}$.
```

```{figure-md} fig5-chap4

<img src="media/chap4/fig5-chap4.png"  width="600px">

Contours of $\dot\gamma_\text{e}$ for the Stokes solution $\vek u^{(0)}$.
```

```{figure-md} fig6-chap4

<img src="media/chap4/fig6-chap4.png"  width="600px">

Pressure contours of the Stokes solution $p^{(0)}$ for $\mu=1$.
```

Next, we solve the flow for a power law fluid ($m=1$, $n=0.2$, see {numref}`fig7-chap4`) using a Picard iteration scheme with the Stokes flow solution as a starting vector. The same mesh ($N=40$) and elements (Taylor-Hood $Q_2$-$Q_1$) are used.
In {numref}`fig8-chap4`, we show the difference between iterations of the solution sequence $\{\col u^{(j)},\col p^{(j)}\}$, $j=1,2,3,\dots$, where $\Delta_u^{(j)}$ and $\Delta_p^{(j)}$ are defined by

$$
 \Delta_u^{(j)} = \max_i ( |u_i^{(j)}- u_i^{(j-1)}|), \quad 
  \Delta_p^{(j)} = \max_i ( |p_i^{(j)}- p_i^{(j-1)}|), \quad  j=1,2,3,\dots
$$ (eq36-chap4)

Some things to note:

1. The velocity differences $\Delta_u^{(j)}$ converge. Moreover, the convergence is linear, i.e.:

   $$\Delta_u^{(j+1)} = C\,\Delta_u^{(j)} $$  (eq37-chap4)

   where $C$ is a constant and the iteration number $j$ is sufficiently large (around 10).

2. After an initial convergence, the value of $\Delta_p^{(j)}$ fluctuates around $10^{-4}$ and convergence of the pressure solution is lost.

3. The fluctuations in $\Delta_p^{(j)}$ seem to be random, indicating a possible loss of precision in floating-point representation of the matrix and in the subsequent solution of the system.


```{figure-md} fig7-chap4

<img src="media/chap4/fig7-chap4.png"  width="600px">

Viscosity $\eta(\dot\gamma)$ of a power law fluid with $m=1$, $n=0.2$.
```


```{figure-md} fig8-chap4

<img src="media/chap4/fig8-chap4.png"  width="600px">

Difference between iterations for a power law fluid with $m=1$, $n=0.2$.
```

In {numref}`fig9-chap4`-{numref}`fig11-chap4` we show the solution of this problem after 100 iterations: streamlines and $\dot\gamma_\text{e}$ based on $\col{u}^{(100)}$ and pressure contours[^3] of $p^{(100)}$. The imposed motion of the upper boundary generates a vortex with a strength (flow rate) of 0.0635, which is three times smaller than for Stokes. The centre of the vortex is shifted upwards and streamlines become closer near the lid. The “shearing”ß of the flow now becomes even more localised near the lid with a maximum shear rate around 23.9, as shown by the contours of $\dot\gamma_\text{e}$. The pressure is symmetric, but not up to machine precision as shown by the contour halfway between the max and min value of the solution. The latter supports the suspicion of loss of precision in the build and solve of the system matrix, as discussed above.

```{figure-md} fig9-chap4

<img src="media/chap4/fig9-chap4.png"  width="600px">

Streamfunction contours for a power law fluid with $m=1$, $n=0.2$
```

```{figure-md} fig10-chap4

<img src="media/chap4/fig10-chap4.png"  width="600px">

Contours of $\dot\gamma_\text{e}$ for a power law fluid with $m=1$, $n=0.2$.
```

```{figure-md} fig11-chap4

<img src="media/chap4/fig11-chap4.png"  width="600px">

Pressure contours for a power law fluid with $m=1$, $n=0.2$.
```

The source of the loss of precision is the wide range of viscosities in the flow. In fact, from {numref}`fig10-chap4` we see that the range of $\dot\gamma_\text{e}$ is from 0 to 23.9 and therefore the range of viscosities (see {numref}`fig7-chap4`) is from 0.079 to $\infty$(!). The very high viscosities in some regions of the flow only act as a “penalty” for the flow and it seems wise to limit the maximum value. Therefore, we solve the same flow with a Carreau model having $\eta_0=1000$ and a $\lambda$ value according to {numref}`ex:4.4`. The model is shown in {numref}`fig12-chap4` and compared to the power law model.



```{figure-md} fig12-chap4

<img src="media/chap4/fig12-chap4.png"  width="700px">

Viscosity $\eta(\dot\gamma)$ of a power law fluid with $m=1$, $n=0.2$ and a Carreau fluid with $\eta_0=1000$, $n=0.2$ and $\lambda=(\eta_0/m)^{1/(1-n)}$.
```

In {numref}`fig13-chap4`, we show the difference between iterations of the solution sequence $\{\col u^{(j)},\col p^{(j)}\}$, $j=1,2,3,\dots$ for the Carreau model with $\eta_0=1000$. It is clear that we now have also linear convergence for the pressure field, except near very small values of $\Delta_p^{(j)}$ of $10^{-9}$, where we hit another precision loss due to the, still, relatively large differences in viscosities. The plots for streamfunction and $\dot\gamma_\text{e}$ are the same as for the power law model (not shown) and the new pressure field is given in {numref}`fig14-chap4`. The only difference with {numref}`fig11-chap4` is that the contour halfway between the max and min value is now indeed again near machine precision. In {numref}`fig15-chap4` we show the viscosity $\eta(\dot\gamma_\text{e})$ in the flow. Considering that (see {numref}`fig12-chap4`) the Carreau models starts to deviate from the power law model near a viscosity of 500, a considerable part of the centre of the vortex and the lower left and right corner have greatly reduced viscosities as compared to the power law model.

Finally, in {numref}`fig16-chap4` we show the condition number (the ratio of the maximum and minimum eigenvalue) of the system matrix for the first twenty iterations of the power law and the Carreau model as compared to the Stokes system for $\mu=1$. Clearly, the condition of the matrix is much worse for the power-law model than for the Carreau model, supporting the claim that the convergence problems for the pressure are related to precision loss.


````{exercise}
:label: ex:4.6
How do you expect, the following actions:   

1. Refining the mesh.
2. Scaling the system of equations.
3. Using 128 bit floating point (quadruple) instead of 64 bit (double).

will affect the convergence of the pressure for a power law fluid model?

````

```{figure-md} fig13-chap4

<img src="media/chap4/fig13-chap4.png"  width="600px">

Difference between iterations for a Carreau fluid with $\eta_0=1000$, $n=0.2$ and $\lambda=(\eta_0/m)^{1/(1-n)}$.
```

```{figure-md} fig14-chap4

<img src="media/chap4/fig14-chap4.png"  width="600px">

Pressure contours for a Carreau fluid with $\eta_0=1000$, $n=0.2$ and $\lambda=(\eta_0/m)^{1/(1-n)}$.
```

```{figure-md} fig15-chap4

<img src="media/chap4/fig15-chap4.png"  width="600px">

Viscosity contours for a Carreau fluid with $\eta_0=1000$, $n=0.2$ and $\lambda=(\eta_0/m)^{1/(1-n)}$.
```

```{figure-md} fig16-chap4

<img src="media/chap4/fig16-chap4.png"  width="600px">

Condition number of the system matrix for a Stokes fluid ($\mu=1$), a power law fluid with $m=1$, $n=0.2$ and for a Carreau fluid with $\eta_0=1000$, $n=0.2$ and $\lambda=(\eta_0/m)^{1/(1-n)}$. The latter two are shown as a function of the interation number.
```

In {numref}`fig17-chap4` we show the value of $\Delta_u^{(j)}$ for a Carreau model with different value of $n$. Lowering $n$ (more shear-thinning) clearly increases the number of required iterations. In the next section we will discuss the Newton-Raphson iteration to possibly accelerate the iteration process.

```{figure-md} fig17-chap4

<img src="media/chap4/fig17-chap4.png"  width="600px">

Difference between iterations $\Delta_u^{(j)}$ for a Carreau fluid with $\eta_0=1000$, $\lambda=1000^{1.25}=5623.41$, and for various values of $n$.
```
(Chap4.6)=
## Newton-Raphson iteration


(Chap4.6.1)=
### Theory

Newton-Raphson iteration is basically an iteration technique for finding the root(s) of a function $f(x)$, i.e. finding the value(s) of $x$ where $f(x)=0$. Starting from an approximation of a root $x^{(j)}$, a new approximation $x^{(j+1)}$ is obtained by linearizing $f(x)$ around $x^{(j)}$ and finding the root of this linear equation. This is equivalent to constructing the tangent line at $(x^{(j)},f(x^{(j)}))$ and intersecting this line with the $x$-axis (see {numref}`fig18-chap4`).

```{figure-md} fig18-chap4

<img src="media/chap4/fig18-chap4.png"  width="600px">

Newton-Raphson iteration process for $f(x)=0$.
```

Linearization of $f(x)$ around $x^{(j)}$ gives:

$$
  f(x) \approx f(x^{(j)}) + f'(x^{(j)})(x-x^{(j)})
$$ (eq38-chap4)

Setting the linear approximation to zero, gives the new approximation of the root:

$$
   f(x^{(j)}) + f'(x^{(j)})(x^{(j+1)}-x^{(j)})=0 \qquad\text{or}\qquad x^{(j+1)}=x^{(j)}-\frac{f(x^{(j)})}{f'(x^{(j)})}
$$ (eq39-chap4)

The Newton-Raphson iteration process is now defined by a starting value $x^{(0)}$, the sequence of solutions $x^{(j)}$, $j=1,2,3,\dots$ and a stopping criterium, e.g. $|x^{(j)}-x^{(j-1)}|\leq\epsilon$, with $\epsilon$ a small positive value. Convergence is obtained if the starting solution $x^{(0)}$ is close enough to the root. If the multiplicity of the root is one, the convergence is quadratic provided the iterations $x^{(j)}$ are sufficiently close  to the root and the function $f(x)$ is sufficiently smooth. Quadratic convergence means that the error can be written as

$$
  |x^{(j+1)}-\alpha| < C_1 |x^{(j)}-\alpha|^2
$$ (eq40-chap4)

where $\alpha$ is the root and $C_1$ is a constant. In practise, it also means that the difference between iterations converges quadratically:

$$
  |x^{(j+1)}-x^{(j)}| < C_2 |x^{(j)}-x^{(j-1)}|^2
$$ (eq41-chap4)

where $C_2$ is a constant.

```{admonition} Remark 4.7
:class: note

As we shall see in the next subsection, the window of convergence and/or quadratic convergence can be quite small. Finding a starting vector that leads to convergence is usually a process of trial and error.

```

In practice, we do not have a function of a single variable, but multiple functions of multiple variables, i.e.\ we are searching for the solution of the nonlinear system of equations $\col{f}(\col{u})=\col{0}$. Again we can define the Newton-Raphson iteration process by linearization around an approximation $\col{u^{(j)}}$ of the nonlinear system of equations:

$$
  \col{f}(\col{u}) \approx \col{f}(\col{u^{(j)}}) + \mat{J}(\col{u^{(j)}})(\col{u}-\col{u^{(j)}})
$$ (eq42-chap4)

where the matrix $\mat J(\col u^{(j)})$ is the Jacobian at $\col u^{(j)}$:

$$
   \mat{J}(\col{u^{(j)}}) = \pderiv{\smash{\col{f}(\col{u})}}{\col{u}}(\col{u}^{(j)})
$$ (eq43-chap4)

Setting the linear approximation to zero, gives the new approximation of the root:

$$
  \col{f}(\col{u^{(j)}}) + \mat{J}(\col{u^{(j)}})(\col{u^{(j+1)}}-\col{u^{(j)}}) = \col{0} \qquad\text{or}\qquad 
  \col{u}^{(j+1)}=\col{u^{(j)}}-\mat{J}(\col{u^{(j)}})^{-1}\col{f}(\col{u^{(j)}}).
$$ (eq44-chap4)

The Newton-Raphson iteration process is now defined by a starting value $\col{u^{(0)}}$, the sequence of solutions $\col{u^{(j)}}$, $j=1,2,3,\dots$ and a stopping criterium, e.g. $\Delta_u^{(j)}=\max_i|u_i^{(j)}-u_i^{(j-1)}|\leq\epsilon$, with $\epsilon$ a small positive value. Convergence is obtained if the starting solution $\col{u^{(0)}}$ is close enough to the actual solution.

### Linearization of the weak form

In order to find the linearized set of equations we need an expression for the Jacobian matrix $\mat{J}(\col{u})$. One way of doing this is taking direct derivatives with respect to $\col{u}$ of the discretized nonlinear set of equations, here given by Equations {eq}`eq23-chap4`-{eq}`eq24-chap4`. However, it is often easier to linearize the weak form (in this case Equations {eq}`eq20-chap4`), before the discretization step. That is what we will do here.

To make the linearization easier to read, we use the short-hand notation:

$$
\begin{alignat*}{2}
   \vek u^{(j)}&\rightarrow \vek u &p^{(j)}&\rightarrow p \\
   \vek u^{(j+1)}&\rightarrow \vek u +\Delta\vek u \qquad& p^{(j+1)}&\rightarrow p + \Delta p 
\end{alignat*}
$$ (eq45-chap4)

and linearize all terms with respect to $\Delta \vek u$ and $\Delta p$. In order to expand and linearize the first term in the weak form Equation {eq}`eq20-chap4` (first equation), we need to expand $\eta(\dot\gamma_\text{e})$ first. Since, $\dot\gamma_\text{e}$ depends directly on $\ten D$, we need the derivative of $\eta(\dot\gamma_\text{e})$ with respect to $\ten D$:

$$
 \pderiv{\dot\gamma_\text{e}}{\ten D}=\pderiv{\sqrt{2\bar{I\!I}(\ten D)}}{\ten D}=
     \frac12\big(2\bar{I\!I}(\ten D)\big)^{-\frac12}2\pderiv{\bar{I\!I}(\ten D)}{\ten D}=
     \dfrac{2}{\dot\gamma_\text{e}}\ten D
$$ (eq46-chap4)

where we have used the expression in Equation {eq}`eq2-A` for the derivative of the second moment of $\ten D$. Now we can expand the viscosity function $\eta(\dot\gamma_\text{e})$ as follows[^4]:

$$
  \begin{split}
 \eta\big(\dot\gamma_\text{e}(\ten D+\Delta \ten D)\big)&=\eta\big(\dot\gamma_\text{e}(\ten D)\big)+
     \eta'(\dot\gamma_\text{e})\pderiv{\dot\gamma_\text{e}}{\ten D}:\Delta \ten D+ \text{h.o.t.}\\
     &=\eta\big(\dot\gamma_\text{e}(\ten D)\big)+
     2\frac{\eta'(\dot\gamma_\text{e})}{\dot\gamma_\text{e}}\ten D:\Delta \ten D+ \text{h.o.t.}\\
     &=\eta\big(\dot\gamma_\text{e}(\ten D)\big)+
     2\frac{\eta'(\dot\gamma_\text{e})}{\dot\gamma_\text{e}}\ten D:\nao(\Delta \vek u)+ \text{h.o.t.}
  \end{split}
$$ (eq47-chap4)

where $\eta'=\lderiv{\eta}{\dot\gamma_\text{e}}$ and

$$
  2\ten D=2\ten D(\vek u)=\nao\vek u+(\nao\vek u)^T, \qquad 2\Delta \ten D=2\ten D(\Delta \vek u)=\nao\Delta \vek u+(\nao\Delta \vek u)^T
$$ (eq48-chap4)

Applying the weak form Equation {eq}`eq20-chap4` to the solution $\vek u +\Delta\vek u$ and $p + \Delta p$, substituting Equation {eq}`eq47-chap4` and neglecting higher-order terms, we get the following linearized weak form: find $\Delta \vek u\in \vek U$, $\Delta p\in P$ such that

$$
\begin{gather*}
\begin{split}
\int_\Omega 4\frac{\eta'(\dot\gamma_\text{e})}{\dot\gamma_\text{e}}(\nao\vek v)^T:\ten D(\vek u)
\ten D(\vek u):\nao\Delta \vek u\,d\Omega\hspace{5cm}&\\
+\int_\Omega 2\eta(\dot\gamma_\text{e})(\nao\vek v)^T:\ten D(\Delta \vek u)\,d\Omega
-\int_\Omega\nao\cdot\vek v\Delta p\,d\Omega=\\
-\int_\Omega 2\eta(\dot\gamma_\text{e})(\nao\vek v)^T:\ten D(\vek u)\,d\Omega
+\int_\Omega\nao\cdot\vek vp\,d\Omega
        +\int_{\Gamma_\text{N}} \vek v\cdot\vek t_\text{N}\,d\Gamma +\int_\Omega \vek v\cdot\vek f\,d\Omega
\end{split}
\\
  \int_\Omega q\nao\cdot\Delta\vek u\,d\Omega = -\int_\Omega q\nao\cdot\vek u\,d\Omega
\end{gather*}
$$ (eq49-chap4)

for all $\vek v\in \vek V$ and $q\in Q$.

```{admonition} Remark 4.8
:class: note

The system is linear in pressure. Therefore, we can combine the two pressure terms in Equation {eq}`eq49-chap4` into one term: 

$$
 \dots -\int_\Omega\nao\cdot\vek v(p+\Delta p)\,d\Omega = \dots
$$ (eq50-chap4)

In this way, we solve for $p^{(j+1)}=p+\Delta p$, instead of the pressure difference $\Delta p$.
```
```{admonition} Remark 4.9
:class: note

If the starting velocity $\vek u^{(0)}$ of the iteration process is a solution of the Stokes problem the right-hand side of Equation {eq}`eq49-chap4` (last line) will be zero (why?). In that case, the right-hand side might be set to zero to make the system simpler.
However, in practise it will not be exactly zero and the error will propagate in the iteration process. Therefore, we will keep the right-hand side of Equation {eq}`eq49-chap4` (last line).

```

````{exercise}
:label: ex:4.7

 Show that for the models given in Remark~\ref{rem:genmodels}, $\eta'/\dot\gamma$ is given by:
 $$
  \begin{alignat*}{3}
   \frac{\eta'(\dot\gamma)}{\dot\gamma} &~= -\dfrac{m(1-n)}{\dot\gamma^{2-n}} & \quad&\text{Power-law}\\
   \frac{\eta'(\dot\gamma)}{\dot\gamma} &=  
-\dfrac{(\eta_0-\eta_\infty)(1-n)\lambda^2}{\big(1+(\lambda\dot\gamma)^2\big)^{(3-n)/2}} &
          \quad&\text{Carreau}\\
   \frac{\eta'(\dot\gamma)}{\dot\gamma} &= 
-\dfrac{(\eta_0-\eta_\infty)(1-n)\lambda^a \dot\gamma^{a-2}}{\big(1+(\lambda\dot\gamma)^a\big)^{(1+a-n)/a}} &
          \quad&\text{Carreau-Yasuda}
\end{alignat*}
$$ (eq51-chap4)

Study the behavior for $\lim_{\dot\gamma\rightarrow0}\eta'/\dot\gamma$ and do you expect this behavior to lead to convergence problems of the Newto1n-Raphson scheme?

````

### Galerkin method for the linearized weak form

The discretized equations can be obtained by replacing 
the general spaces with the finite dimensional approximation spaces in the weak form {eq}`eq49-chap4`.

The approximation spaces for velocity and pressure can be taken the same as for the Stokes equations, at least if the monotonicity requirement Equation {eq}`eq13-chap4` is being fulfilled.cUsing Equations {eq}`eq18-chap4` together with Equations {eq}`eq21-chap4` the following set of discretized equations can be derived:

$$
\begin{gather*} 
\begin{split}
    \sum_{m=1}^{N_u}\big(\ten G_{km}(\vek u_h)+\ten K_{km}(\vek u_h)\big)\cdot \Delta \vek u_m+ \sum_{m=1}^{N_p} \vek b_{mk}\Delta p_m\\
  =-\vek g_{k}(\vek u_h)-\sum_{m=1}^{N_p}& \vek b_{mk}p_m+\vek f_k, \qquad k=1,\dots N_{u}
\end{split}\\
  \sum_{m=1}^{N_u}\vek b_{km}\cdot \Delta \vek u_m=-\sum_{m=1}^{N_u}\vek b_{km}\cdot \vek u_m, \qquad k=1,\dots N_{p}
\end{gather*}
$$ (eq52-chap4)

where $\vek b_{km}$ and $\vek f_k$ are the same as for Stokes (Equations {eq}`eq21-chap3`-{eq}`eq22-chap3`), the vectors $\vek g_{k}(\vek u_h)$ are defined in Equation {eq}`eq26-chap4`, the tensors $\ten K_{km}$ are defined in Equation {eq}`eq28-chap4` and the tensors $\ten G_{km}$ are defined by

$$
 \ten G_{km}(\vek u_h) = \int_\Omega 4\frac{\eta'(\dot\gamma_\text{e})}{\dot\gamma_\text{e}}
   \nao\phi_k\cdot\ten D(\vek u_h)\ten D(\vek u_h)\cdot\nao\phi_m\,d\Omega \qquad
  k,m=1,\dots N_{u}
$$ (eq53-chap4)

````{exercise}
:label: ex:4.8

Derive the expression for $\ten G_{k}(\vek u_h)$ in Equation {eq}`eq53-chap4`.

````

````{exercise}
:label: ex:4.9

Show that the following terms in the right-hand side of Equations {eq}`eq52-chap4` can be evaluated as:

$$
\begin{align*}
-\sum_{m=1}^{N_p} \vek b_{mk}p_m &= \int_\Omega p_h\nao \phi_k \, d\Omega \\
-\sum_{m=1}^{N_u}\vek b_{km}\cdot \vek u_m & = \int_\Omega \psi_k\nao\cdot\vek u_h \, d\Omega
\end{align*}
$$ (eq54-chap4)

````

The system Equations {eq}`eq52-chap4` can now be written as 

$$
\begin{align*}
[\mat{G}(\col{u})+ \mat{K}(\col{u})]\Delta \col{u} +\mat{B}^T\Delta\col{p} &= -\col{g}(\col{u}) -\mat{B}^T\col{p}+ \col{f} \\
   \mat{B}\Delta \col{u}  &= -\mat{B}\col{u}
\end{align*}
$$ (eq55-chap4)

or in full matrix form

$$
\begin{pmatrix}
  \mat{G}(\col{u})+\mat{K}(\col{u}) & \mat{B}^T\\
  \mat{B} & \mat{0}
\end{pmatrix}
\begin{pmatrix}
 \Delta \col{u} \\
 \Delta \col{p}
\end{pmatrix}=
\begin{pmatrix}
 -\col{g}(\col{u}) -\mat{B}^T\col{p}+ \col{f} \\
  -\mat{B}\col{u}
\end{pmatrix}
$$ (eq56-chap4)

where the matrix $\mat{K}(\col{u})$ has already been defined in Equation {eq}`eq30-chap4` and the matrix $\mat G(\col u)$ is defined similarly:

$$
   \mat{G}(\col{u})=  
  \begin{pmatrix}
    \mat{G}_{11} & \mat{G}_{12} & \cdots & \mat{G}_{1N_u} \\
    \mat{G}_{21} & \mat{G}_{22} & \cdots & \mat{G}_{2N_u} \\
      \vdots & \vdots & \ddots & \vdots \\
    \mat{G}_{N_u1} & \mat{G}_{N_u2} & \cdots & \mat{G}_{N_uN_u}
  \end{pmatrix}
$$ (eq57-chap4)

with $\mat{G_{km}}$ the $d\times d$ matrices of the components of $\ten G_{km}(\vek u_h)$ with respect to a Cartesian coordinate system given by the unit vectors $\vek e_i$, $i=1,\dots,d$. 

````{exercise}
:label: ex:4.10

Show that the matrix $\mat{G}(\col{u})$ is symmetric.

````

````{exercise}
:label: ex:4.11

 Show that the system matrix given in Equation {eq}`eq55-chap4` is symmetric.

````

The *Newton-Raphson iteration scheme* is now defined by a starting vector $\col{u^{(0)}}$ for the velocity and a sequence of solutions $\{\col{u^{(j)}},\col{p^{(j)}\}}$, $j=1,2,3,\dots$, given by

$$
 \col{u^{(j+1)}} =\col{u^{(j)}} + \Delta \col{u}, \qquad
 \col{p^{(j+1)}} =\col{p^{(j)}} + \Delta \col{p}, \qquad j=0,1,2,\dots
$$ (eq58-chap4)

with the iteration increments $\Delta \col u$ and $\Delta \col p$ obtained from


$$
\begin{align*}
[\mat{G}(\col{u^{(j)}})+ \mat{K}(\col{u^{(j)}})]\Delta \col{u} +\mat{B}^T\Delta\col{p} &= -\col{g}(\col{u^{(j)}}) -\mat{B}^T\col{p^{(j)}}+ \col{f}\\
   \mat{B}\Delta \col{u}  &= -\mat{B}\col{u^{(j)}}
\end{align*}
$$ (eq59-chap4)

For the starting vector we can use a Stokes solution as discussed in Remark~\ref{rem:startingvectoru0}.
A stopping criterium, such as $\Delta_u^{(j)}=\max_i|u_i^{(j)}-u_i^{(j-1)}|\leq\epsilon_u$ and $\Delta_p^{(j)}=\max_i|p_i^{(j)}-p_i^{(j-1)}|\leq\epsilon_p$ with $\epsilon_u$ and $\epsilon_p$ small positive values, needs to be defined to limit the number of iterations.

````{exercise}
:label: ex:4.12

 Show that the Newton-Raphson process can also be written as
 a starting vector $\col{u^{(0)}}$ for the velocity and a sequence of solutions $\{\col{u^{(j)}},\col{p^{(j)}\}}$, $j=1,2,3,\dots$, given by

$$
 \col{u^{(j+1)}} =\col{u^{(j)}} + \Delta \col{u}, \qquad j=0,1,2,\dots
$$ (eq60-chap4)

with the iteration increment $\Delta \col u$ and the pressure $\col{p^{(j+1)}}$ obtained from

$$
\begin{align*}
[\mat{G}(\col{u^{(j)}})+ \mat{K}(\col{u^{(j)}})]\Delta \col{u} +\mat{B}^T\col{p}^{(j+1)} &= -\col{g}(\col{u^{(j)}})+ \col{f}\\
   \mat{B}\Delta \col{u}  &= -\mat{B}\col{u^{(j)}}
\end{align*}
$$ (eq61-chap4)

````

````{exercise}
:label: ex:4.13

If you would use a stabilized element (see {numref}`Chap3.6`), which formulation would you prefer: solving for $(\Delta u,\Delta p)$ (Equations {eq}`eq59-chap4`) or $(\Delta u,p)$ (Equations {eq}`eq61-chap4`)? Explain why.

````

(Chap4.6.4)=
### Convergence of the Newton-Raphson iteration: lid-driven cavity flow

We now solve the driven cavity problem using Newton-Raphson iteration. The
iteration process for this problem using the Picard process is shown in {numref}`fig13-chap4`. Just replacing the Picard process with a Newton-Raphson process is insufficient and divergence is obtained. Apparently, the starting vector obtained using a Stokes problem is outside the convergence window of the Newton-Raphson process for this problem. Hence, we need to get closer to the actual solution of the problem to obtain convergence of the Newton-Raphons iteration process. We can use Picard iteration for that, since we have shown that this is a converging process for this problem. In {numref}`fig19-chap4` we show three different iteration processes:

1. only Picard iteration (the data from {numref}`fig13-chap4`),

2. a process of eight Picard iterations followed by Newton-Raphson,

3. a process of nine Picard iterations followed by Newton-Raphson.

The second one diverges. The third one, however, converges fast as might be expected from the Newton-Raphson theory. After the switch from Picard to Newton-Raphson there is an adaptation with relatively large differences between iterations, but still turns into a fast converging process for the third case. Note, that also here the pressure difference between iteration stops converging and remains at the same level similar to the Picard scheme due to precision loss. The velocity difference also stops converging, but this is close 


```{figure-md} fig19-chap4

<img src="media/chap4/fig19-chap4.png"  width="600px">

Difference between iterations for a Carreau fluid with $\eta_0=1000$, $n=0.2$ and $\lambda=(\eta_0/m)^{1/(1-n)}$. Shown are: Picard, eight Picard iterations followed by Newton-Raphson iteration and nine Picard iterations followed by Newton-Raphson iteration. 
```

```{figure-md} fig20-chap4

<img src="media/chap4/fig20-chap4.png"  width="600px">

zoom-in to the first 25 iterations.
```

The number of Picard iterations needed to obtain convergence for the Newton-Raphson process depends on the value of $n$. This is expected, since the Picard process convergence is slower for smaller $n$ and it takes more iterations to reach within a certain distance from the actual solution.  
In {numref}`tab1-chap4` we give the number of Picard iterations needed for various values of $n$. It is clear that the number increases significantly for the highly shear-thinning fluids. As soon as this number is reached, convergence is fast, however. For values $n\geq0.5$ Newton-Raphson converges without Picard iterations, however doing a few Picard iterations (say 2-4) anyway helps in lowering the actual needed number of total iterations.


```{list-table} Number of Picard iterations *before* the Newton Raphson iteration process needed to obtain convergence for a Carreau fluid with $\eta_0=1000$, $\lambda=1000^{1.25}=5623.41$ for various values of $n$.
:header-rows: 1
:name: tab1-chap4

* - $n$
  - 0.01
  - 0.02
  - 0.05
  - 0.1
  - 0.2
  - 0.3
  - 0.4
  - $n \ge 0.5$
* - Picard iterations
  - 110
  - 70
  - 48
  - 15
  - 9
  - 5
  - 3
  - 0

```



[^1]: The viscosity function $\eta(\dot\gamma)$ is measured in a steady shear flow  as a function of the shear rate $\dot\gamma$.

[^2]: We actually plot $p-1/A\int p\,d\Omega$.

[^3]: Again we plot $p-1/A\int p\,d\Omega$.

[^4]: h.o.t.$=$ higher-order terms.
