(Chap3)=
# Flow of a Newtonian fluid: the Stokes equations

In this chapter we consider the flow of incompressible Newtonian fluids where inertia can be neglected. These flows can be modelled by the Stokes equations, Equations {eq}`eq49-chap1`. Typical applications
where the Stokes equations can be useful are: flows of highly viscous fluids, such as in the mixing of molten polymers, and flows at very small length scales, such as suspensions.

The finite element discretization of the Stokes equations requires basis functions for both the velocity and pressure fields. This is a non-trivial problem, since the incompressibility constraint Equation {eq}`eq49-chap1` (second equation) does not involve the pressure field. Straightforward application of, for example, equal order basis functions leads to an unstable system. In fact, this so-called *mixed* approximation requires a compatibility condition to be safisfied. The objective of this chapter is to discuss the basic problems and present discretizations that are stable for the Stokes equations and can also be used for the extension to generalized Newtonian fluid models, the Navier-Stokes equations and to viscoelastic fluids. 

## Problem formulation

The momentum and mass balance are given by Equations {eq}`eq42-chap1` and {eq}`eq46-chap1` read

$$
\begin{align*}
 -\nao\cdot\ten\sigma^T&=\vek f&&\text{in }\Omega, \\
   \nao\cdot\vek u&=0&\qquad&\text{in }\Omega.
\end{align*}
$$ (eq1-chap3)

For a Newtonian fluid given by the constitutive equation Equation {eq}`eq47-chap1`, together with the
Cauchy stress expression Equation {eq}`eq44-chap1`, this reduces to the Stokes equations Equations {eq}`eq49-chap1`, which we repeat here:

$$
\begin{align*}
 -\nao\cdot(2\mu \ten D)+\nao p &=\vek f&&\text{in }\Omega \\
   \nao\cdot\vek u&=0&\qquad&\text{in }\Omega.
\end{align*}
$$ (eq2-chap3)

The boundary conditions (for now) are assumed to be

$$
\begin{align*}
    \vek u &= \vek u_\text{D}&&\qquad\text{ on }\Gamma_D,\\
    \ten\sigma\cdot\vek n=-p\vek n+2\mu \ten D\cdot\vek n&=\vek t_{\text{N}} &&\qquad\text{ on }\Gamma_{\text{N}}.
\end{align*}
$$ (eq3-chap3)

where the boundary $\Gamma=\Gamma_\text{D}\cup\Gamma_{\text{N}}$ and has been split into a
*Dirichlet* and a *Neumann* part (see {numref}`fig1-chap3`). The Neumann boundary condition is equivalent to an imposed traction of $\vek t_{\text{N}}$.


<a name="Remark:3.1"></a>
```{admonition} Remark 3.1
:class: note

Recall (see [Remark 2.1](#Remark:2.1)), that Dirichlet and Neumann boundary conditions are also known as essential and natural boundary conditions, respectively.
```
<a name="Remark:3.2"></a>
```{admonition} Remark 3.2
:class: note

Other types of boundary conditions are possible, for example  

- Dirichlet in $x$ and Neumann in $y$,
- Robin boundary conditions: $\vek u + \alpha\dpderiv{\vek u}{n}=\vek \beta$,
- Dirichlet in normal direction: $\vek n\cdot\vek u=\alpha$, Neumann in tangential direction,

however, these will not be considered in the theory, only in practical applications.
```



```{figure-md} fig1-chap3

<img src="media/chap3/fig1-chap3.png"  width="700px">

Domain $\Omega$ with boundary $\Gamma=\Gamma_\text{D}\cup\Gamma_{\text{N}}$, outwardly directed unit normal vector $\vek n$ on $\Gamma$, imposed velocity $\vek u_\text{D}$ on $\Gamma_\text{D}$ and imposed traction vector $\vek t_\text{N}$ on the Neumann boundary.
```

In order to obtain a unique solution for the velocity $\vek u$, the boundary $\Gamma_\text{D}$ should have a finite size. Otherwise, ``rigid translational and rotational modes'' would lead to a singular system. On the other hand, the quite common case of the whole boundary $\Gamma$ being Dirichlet ($\Gamma_\text{N}$ is absent) leads to a pressure field $p$ that is only unique up a constant. This constant (level) can be specified by an additional condition, such as specifying the pressure in a single point or imposing that

$$
   \int_\Omega p \,d\Omega = 0
$$

Furthermore, there is a *compatibility restriction* on the function $\vek u_\text{D}$ in Equation {eq}`eq3-chap3` (first equation). Integrating Equation {eq}`eq2-chap3` (second equation) over the domain $\Omega$ and applying the Gauss theorem Equation {eq}`eq4-A` gives:

$$
    0 = \int_\Omega \nao\cdot \vek u\,d\Omega = \int_\Gamma \vek n\cdot \vek u\,d\Gamma = 
   \int_\Gamma \vek n\cdot \vek u_\text{D}\,d\Gamma
$$ (eq4-chap3)

Specifying that the flowrate through the boundary $\Gamma$, as imposed by $\vek u_\text{D}$, needs to be zero. Or in other words: the inflow must be equal to outflow. There can be no storage inside the domain $\Omega$.
Obviously, no solution exists if $\vek u_\text{D}$ does not fulfill the compatibility condition Equation {eq4-chap3}.

````{exercise}
:label: ex:3.1

Argue, that also for the boundary condition in normal direction (item 3 of [Remark 3.2](#Remark:3.2)) the level of the pressure needs to be specified and the compatibility condition Equation {eq}`eq4-chap3` must be fulfilled.

````

We should note here, that a solution fulfilling the PDE is called a *classical* or *strong solution*. In order for this solution to be called classical, the right-hand side $\vek f$ needs to be a continuous function and the solution $\vek u$ needs to have continuous second-order derivatives. If this requirement can't be met, for example
$\vek f$ contains jumps in space or the boundary is not smooth enough, we needs to enlarge the solution space by fulfilling the PDE in a weak sense and turn to the so-called *weak formulation* of the PDE. In the next section we will present a weak formulation of the Stokes equations.

(Chap3.2)=
## Weak formulation

The Stokes system Equations {eq}`eq2-chap3` is a system in the unknown fields (functions) velocity $\vek u(\vek x,t)$ and pressure $p(\vek x,t)$. When deriving the weak form we need to be aware, that these functions require proper function spaces. As in Chapter {numref}`Chap2` we want to keep it simple and avoid the choice of proper function spaces as much as possible by saying that they need to be as big as possible as allowed by the theory. We will denote the corresponding function spaces by $\vek U$ and $P$, respectively. The solutions spaces $\vek U$ and $P$ are also called *trial spaces* for the velocity and pressure fields, respectively.

In order to derive the weak form of the Stokes equations, we could start directly from the Stokes system 
Equations {eq}`eq2-chap3` and follow a similar approach as in {numref}`chap2.3`. However, starting from the general system Equations {eq}`eq1-chap3` leads to a generic weak form that can be used for other fluid models as well. Multiplying the equations from {eq}`eq1-chap3` with test functions $\vek v$ and $q$, respectively, and integrating
over the domain $\Omega$ we obtain

$$
\begin{align*}
\int_\Omega \vek v\cdot(-\nao\cdot\ten\sigma^T-\vek f)\,d\Omega &=0\\
\int_\Omega q\nao\cdot\vek u\,d\Omega &=0
\end{align*}
$$ (eq5-chap3)

for all test functions $\vek v\in \vek V$ and $q\in Q$. We have denoted the function spaces for $\vek v$ and $q$ by $\vek V$ and $Q$, respectively. The spaces $\vek V$ and $Q$ are called *test spaces* of the velocity and pressure fields, respectively. 

<a name="Remark:3.3"></a>
```{admonition} Remark 3.3
:class: note

The procedure above for deriving a weak form is also known as the method of weighted residuals (see [Remark 2.9](#Remark:2.9)). 
The residual, i.e. the left-hand side minus the right-hand side,  is multiplied with a test function and integrated over the domain. Using the notation of [Remark 2.1](#Remark:2.8), Equation {eq}`eq5-chap3` can be written as

$$
    (\vek v, \vek r),\qquad\text{with }\vek r=-\nao\cdot\ten\sigma^T-\vek f
$$

for all test functions $\vek v\in \vek V$, where the residual $\vek r$ follows from Equation {eq}`eq1-chap3`. Another way of expression this is: the residual is required to be orthogonal to all functions from the test space $\vek V$.
```

It is quite common to impose the Dirichlet condition {eq}`eq3-chap3` in a strong way on the trial space $\vek U$, i.e. functions from this trial space should already fulfill the Dirichlet condition (see also Sec.~\ref{sec:weak_poisson} for the more simpler case using the diffusion equation). In that case, 
we must also restrict the test space $\vek V$ in the following way:

$$
\vek v = \vek 0\qquad\text{ on }\Gamma_D
$$ (eq6-chap3)

Assuming that a strong solution for $\vek u$ and $p$ exists that fulfills the PDE, the formulation in Equations {eq}`eq5-chap3` is fully equivalent to the original PDE. This is due to the fundamental lemma of variational calculus (see {numref}`chap2.3`). However, the (weak) formulation in Equations {eq}`eq5-chap3` allows solutions that go beyond the classical solution space. The reason is that the integral formulation only requires the integrals to exist, allowing discontinuous functions. For example, $\vek f$ can now have jumps in space and the second derivatives of $\vek u$ only needs to be (square) integrable. However, for application to finite elements these function requirements are still too high.

In order to lower the function requirements, we partially integrate the weak form of the momentum balance Equation {eq}`eq5-chap3` (first eqaution). For that, we use Equation {eq}`eq3-A` (9th equation) to find that

$$
\begin{align*}
\vek v\cdot(\nao\cdot\ten\sigma^T)&=(\nao\cdot\ten\sigma^T)\cdot\vek v\\
&=\nao\cdot(\ten\sigma^T\cdot\vek v)-\ten\sigma^T:(\nao\vek v)^T\\                                         
&=\nao\cdot(\ten\sigma^T\cdot\vek v)-(\nao\vek v)^T:\ten\sigma^T
\end{align*}
$$

and after substitution into {eq}`eq5-chap3` (first equation) we get

$$
  -\int_\Omega \nao\cdot(\ten\sigma^T\cdot\vek v)\,d\Omega +\int_\Omega (\nao\vek v)^T:\ten\sigma^T\,d\Omega=\int_\Omega \vek v\cdot\vek f\,d\Omega
$$ (eq7-chap)

Applying the Gauss theorem of {eq}`eq4-A` leads to

$$
-\int_\Gamma \vek v\cdot(\ten\sigma\cdot\vek n)\,d\Gamma +\int_\Omega (\nao\vek v)^T:\ten\sigma^T\,d\Omega=\int_\Omega \vek v\cdot\vek f\,d\Omega
$$ (eq8-chap3)

<a name="Remark:3.4"></a>
```{admonition} Remark 3.4
:class: note

The second term in the last equation could have been simplified, as follows

$$
  \int_\Omega (\nao\vek v)^T:\ten\sigma^T\,d\Omega=\int_\Omega \nao\vek v:\ten\sigma\,d\Omega
$$

however, the first version more easily can be used to put $\vek v$ ``in front of'' the equation when discretizing. 
```

````{exercise}
:label: ex:3.2

Show that, when using the symmetry property of $\ten\sigma$, the second term in {eq}`eq8-chap3` can be written as

$$
  \int_\Omega (\nao\vek v)^T:\ten\sigma^T\,d\Omega=\int_\Omega \ten D(\vek v):\ten\sigma\,d\Omega
$$

with $\ten D(\vek v)=(\nao\vek v+(\nao\vek v)^T)/2$.

````

As a final step we fill in the Neumann (traction) boundary condition of $\Gamma_\text{N}$ and the test function restriction {Eq}`eq6-chap3` on $\Gamma_{\text{D}}$ to arrive at the following weak form of the momentum balance

$$
  \int_\Omega (\nao\vek v)^T:\ten\sigma^T\,d\Omega=\int_{\Gamma_\text{N}} \vek v\cdot\vek t_\text{N}\,d\Gamma +\int_\Omega \vek v\cdot\vek f\,d\Omega 
$$ (eq9-chap3)

for all $\vek v\in \vek V$. This weak form of the momentum balance can be used for various constitutive laws, not just Newtonian fluids.

````{exercise}
:label: ex:3.3

 Show that, assuming a strong solution exists, Equation {Eq}`eq9-chap3` can be rewritten to

$$
  \int_{\Gamma_\text{N}} \vek v\cdot(\ten\sigma\cdot\vek n-\vek t_\text{N})\,d\Gamma+
  \int_\Omega \vek v\cdot(-\nao\cdot\ten\sigma^T-\vek f)\,d\Omega = 0
$$

for all $\vek v\in \vek V$. Argue, that this is equivalent to the original momentum balance Equation {eq}`eq1-chap3` (first equation) and the Neumann boundary conditions {eq}`eq3-chap3` (second line).

````

Substituting {eq}`eq44-chap2` and {eq}`eq47-chap2` into Equation {eq}`eq9-chap3`, combined with
the weak form of the continuity equation {eq}`eq5-chap3` (second equation), gives the following *weak form of the Stokes equations*: find $\vek u\in \vek U$, $p\in P$ such that

$$
\begin{align*}
 \int_\Omega \mu(\nao\vek v)^T:(\nao\vek u+(\nao\vek u)^T)\,d\Omega-\int_\Omega\nao\cdot\vek vp\,d\Omega
        &=\int_{\Gamma_\text{N}} \vek v\cdot\vek t_\text{N}\,d\Gamma +\int_\Omega \vek v\cdot\vek f\,d\Omega\\
  \int_\Omega q\nao\cdot\vek u\,d\Omega &=0
\end{align*}
$$ (eq10-chap3)

for all $\vek v\in \vek V$ and $q\in Q$.


<a name="Remark:3.5"></a>
```{admonition} Remark 3.5
:class: note

Limiting the trial and test functions to *divergence free* functions, i.e. to $\hat{\vek U}=\{\vek u\in \vek U\mid\nao\cdot\vek u=0\}$ and $\hat{\vek V}=\{\vek v\in \vek V\mid\nao\cdot\vek v=0\}$ reduces the weak form to: find $\vek u\in \hat{\vek U}$  such that

$$
 \int_\Omega \mu(\nao\vek v)^T:(\nao\vek u+(\nao\vek u)^T)\,d\Omega
        =\int_{\Gamma_\text{N}} \vek v\cdot\vek t_\text{N}\,d\Gamma +\int_\Omega \vek v\cdot\vek f\,d\Omega
$$ (eq11-chap3)

for all $\vek v\in \hat{\vek V}$.
Note, that the pressure $p$ (and the test function $q$) completely disappears. This weak form can be used to show that the velocity field $\vek u$ can be found from the *contrained minimization* problem {cite}`Bochev2009`

$$
 \min_{\vek u\in \hat{\vek U}}\left(\frac12\int_\Omega \mu(\nao\vek u+(\nao\vek u)^T):(\nao\vek u+(\nao\vek u)^T)\,d\Omega
  -\int_{\Gamma_\text{N}} \vek v\cdot\vek t_\text{N}\,d\Gamma -\int_\Omega \vek v\cdot\vek f\,d\Omega\right)
$$ (eq12-chap3)

From the full weak form Equations {eq}`eq10-chap3` it can be shown that the velocity $\vek u$ and pressure $p$ can be found from the *saddle point problem* {cite}`Bochev2009`

$$
\begin{multline*}
 \min_{\vek u\in \vek U}\max_{p\in P}\left(\frac12\int_\Omega \mu(\nao\vek u+(\nao\vek u)^T):(\nao\vek u+(\nao\vek u)^T)\,d\Omega-\int_\Omega\nao\cdot\vek up\,d\Omega\right.\\
        \left.-\int_{\Gamma_\text{N}} \vek v\cdot\vek t_\text{N}\,d\Gamma -\int_\Omega \vek v\cdot\vek f\,d\Omega\right)
\end{multline*}
$$ (eq13-chap3)

i.e. searching for a minimum in $\vek u$ and a maximum in $p$, hence we're looking for the stationary point in “a horse's saddle”. Comparing Equation {eq}`eq13-chap3`with Equation {eq}`eq12-chap3`, we see that the pressure field $p$ can be associated with a *Lagrange multiplier* field for imposing the incompressibility constraint.
```

<a name="Remark:3.6"></a>
```{admonition} Remark 3.6

Using the notation for inner products introduced in [Remark 2.8](#Remark:2.8),
the weak form of the Stokes equations can be written as: find $\vek u\in \vek U$, $p\in P$ such that

$$
\begin{align*}
 \mu\left((\nao\vek v)^T,\nao\vek u+(\nao\vek u)^T\right)-\left(\nao\cdot\vek v,p\right)
        &=\left(\vek v,\vek t_\text{N}\right)_{\Gamma_\text{N}}+\left(\vek v,\vek f\right)\\
   \left(q,\nao\cdot\vek u\right)&=0
\end{align*}
$$ (eq14-chap3)

for all $\vek v\in \vek V$ and $q\in Q$.

```

<a name="Remark:3.7"></a>
```{admonition} Remark 3.7

The pressure $p$ and the corresponding test function $q$ need to be in the $L^2(\Omega)$, whereas the velocity components need to be in the $H^1(\Omega)$, 
 since there are gradients of the velocity (see [Remark 2.11](#Remark:2.11) for the definition of $L^2(\Omega)$ and $H^1(\Omega)$).
In engineering practice it means that functions for the velocity components need to be continuous ($C^0(\Omega)$), since jumps would lead to Dirac delta functions for the derivatives. Using the above definition, we can define $P=Q=L^2(\Omega)$ and

$$
\begin{align*}
   \vek U&=\lbrace \vek u\in H^1(\Omega)^d\,\vert\, \vek u=\vek u_\text{D}\text{ on }\Gamma_\text{D}\rbrace\\
   \vek V&=\lbrace \vek v\in H^1(\Omega)^d\,\vert\, \vek v=\vek 0\text{ on }\Gamma_\text{D}\rbrace
\end{align*}
$$

where $H^1(\Omega)^d$ means that each component of the vector (in dimension $d=2$ or $3$) belongs to $H^1(\Omega)$.
Note, that $\vek U$ is only a linear (or vector) space if $\vek u_\text{D}=\vek 0$, since for $\vek u_\text{D}\neq\vek 0$ addition of two elements from $\vek U$ will result in a function not in $\vek U$. The test space $\vek V$ is a vector space, though.

```

(Chap3.3)=
## The Galerkin method

The objective of this section is to obtain a 
finite set of discretized equations, approximating the original Stokes differential equations with boundary conditions. For that, we will use the Galerkin discretization technique. See {numref}`Chap2.5` for the introduction of the Galerkin technique where it is applied to the steady diffusion equation.

The first step of the Galerkin technique consists of limiting the trial spaces $\vek U$ and $P$
to finite dimensional (sub)spaces using a finite set of basis functions. 
The second step consists of limiting the test spaces
$\vek V$ and $Q$ to finite dimensional (sub)spaces using the *same set of basis functions*.
This has been explained in the context of the steady diffusion equation in Section {numref}`Chap2.5`. The finite dimensional approximation spaces will be denoted by the subscript $h$, i.e. $\vek U_h$, $P_h$, $\vek V_h$ and $Q_h$, respectively.

For defining the approximation trial spaces $\vek U_h$ and $P_h$ we write the elements of these spaces as a linear combination of the *basis functions* $\phi_k(\vek x)$, $k=1,\dots,N_u$ and $\psi_k(\vek x)$, $k=1,\dots,N_p$, respectively 

$$
  \vek u_h(\vek x)=\sum_{k=1}^{N_{u}} \vek u_k\phi_k(\vek x), \qquad p_h(\vek x)=\sum_{k=1}^{N_{p}}p_k\psi_k(\vek x)
$$ (eq15-chap3)

Similarly the test spaces $\vek V_h$ and $Q_h$ are defined by writing elements of these spaces as a linear combination of the same basis functions as used for $\vek U_h$ and $P_h$:

$$
  \vek v_h(\vek x)=\sum_{k=1}^{N_{u}} \vek v_k\phi_k(\vek x), \qquad q_h(\vek x)=\sum_{k=1}^{N_{p}}q_k\psi_k(\vek x)
$$ (eq16-chap3)

Although the basis functions $\phi_k(\vek x)$ and $\psi_k(\vek x)$ can be quite general, we restrict ourselves here to the *finite element method* having basis functions as described in Section {numref}`Chap2.6`.

<a name="Remark:3.8"></a>
```{admonition} Remark 3.8

All components of the velocity field are defined using the same (scalar) basis functions $\phi_k$. Although most common finite elements can be written in this form, it excludes special basis functions such as based on the normal velocities on the element edges.

```



<a name="Remark:3.9"></a>
```{admonition} Remark 3.9

The number of velocity degrees of freedom is $n_u=d N_u$, where $d$ is the dimension of real space ($d=2$ or $d=3$). Writing the coefficients into components as $\vek u_k=u^k_i\vek e_i$ (using the Einstein summation convention, see Appendix {ref}`Appendix`), an element from the space $\vek U_h$ can be written as

$$
 \vek u_h=\sum_{k=1}^{N_{u}}\vek u_k\phi_k = \sum_{k=1}^{N_{u}} u^k_i\vek e_i\phi_k= 
 \sum_{k=1}^{N_{u}} u^k_i\vek\phi_i^k
$$ (eq17-chap3)

Hence, the $n_u=d N_u$ functions 

$$
 \vek \phi^k_i(\vek x)=\vek e_i\phi_k(\vek x),\quad i=1,\dots d,\quad k=1,\dots,N_u
$$ (eq18-chap3)

can be considered a basis of the $dN_u$ dimensional space $\vek U_h$ (and also $\vek V_h$) in a Cartesian coordinate system.
The number of pressure degrees of freedom $n_p=N_p$.

```

<a name="Remark:3.10"></a>
```{admonition} Remark 3.10

The spaces $\vek U_h$ and $P_h$ are subspaces of $\vek U$ and $P$ ($\vek U_h\subset \vek U$, $P_h\subset P$) and upon increasing the number of degrees of freedom, i.e. $n_u\rightarrow\infty$ and $n_p\rightarrow\infty$, the spaces should approach the real spaces $\vek U$ and $P$.

```


A subtlety still being ignored is that the Dirichlet conditions should be imposed on the space $\vek U$ and the subset requirement of [Remark 3.10](#Remark:3.10) requires this to be true even for the subset $\vek U_h$.
In Section {numref}`Chap2.7.3` it is extensively discussed how this can be achieved. Although the discussion is
in the context of the steady diffusion problem, the Stokes problem can be handled quite similarly.

The discretized Stokes equations can be obtained by replacing 
the general spaces with the finite dimensional approximation spaces in the weak form of the Stokes equations Equations {eq}`eq10-chap3`: 
find $\vek u_h\in \vek U_h$, $p_h\in P_h$ such that

$$
\begin{align*}
 \int_\Omega \mu(\nao\vek v_h)^T:(\nao\vek u_h+(\nao\vek u_h)^T)\,d\Omega-\int_\Omega\nao\cdot\vek v_hp_h\,d\Omega
        &=\int_{\Gamma_\text{N}} \vek v_h\cdot\vek t_\text{N}\,d\Gamma +\int_\Omega \vek v_h\cdot\vek f\,d\Omega\\
  \int_\Omega q_h\nao\cdot\vek u_h\,d\Omega &=0
\end{align*}
$$ (eq19-chap3)

for all $\vek v_h\in \vek V_h$ and $q_h\in Q_h$.  

````{exercise}
:label: ex:3.4

For the case that $\Gamma_\text{D}=\Gamma$ (the whole boundary is a Dirichlet boundary), the pressure $p$ is only unique up to a constant, if no further conditions on the pressure are imposed. Show, that this can only be the case for the numerical solution of the pressure $p_h$ if the compatibility condition Equation {eq}`eq4-chap3` is fulfilled exactly.

````

<a name="Remark:3.11"></a>
```{admonition} Remark 3.11

Another subtlety has been ignored here if (part of) the boundary $\Gamma$ is curved, for example a circle. 
The usual approach is to define the domain using isoparametric mapping using the same basis functions as used for the approximating function space (see Section {numref}`Chap2.6.4`). This means, that the domain is also being approximated and we should have used $\Omega_h$ and $\Gamma_{\text{N},h}$ in the discretized weak form.
The numerical domain and boundary only converges to the real ones after mesh refinement. Like most of the finite element literature, we will ignore this subtlety to avoid too unwieldy formulations.
Note, that in isogeometric analysis {cite}`Cottrell2009` the basis functions consist of NURBS, as used by solid-modelling software. In this way, regular fixed domains can be represented exactly, even for the coarsest discretization.
```

Using Equations {eq}`eq15-chap3` and {eq}`eq16-chap3`  we find that

$$
\begin{alignat*}{2}
\nao\vek u_h&=\sum_{k=1}^{N_u}\nao\phi_k\vek u_k, &\qquad \nao\cdot\vek u_h&=\sum_{k=1}^{N_u}\nao\phi_k\cdot\vek u_k\\
\nao\vek v_h&=\sum_{k=1}^{N_u}\nao\phi_k\vek v_k, &\qquad \nao\cdot\vek v_h&=\sum_{k=1}^{N_u}\nao\phi_k\cdot\vek v_k
\end{alignat*}
$$ (eq20-chap3)

and the discretized weak form Equations {eq}`eq19-chap3` can be written as

$$
\begin{gather*}
  \sum_{k=1}^{N_u} \vek v_k\cdot (\sum_{m=1}^{N_u} \ten A_{km}\cdot \vek u_m + \sum_{m=1}^{N_p} \vek b_{mk}p_m
  -\vek f_k)=0 \\
  -\sum_{k=1}^{N_p} q_k\sum_{m=1}^{N_u}\vek b_{km}\cdot \vek u_m=0
\end{gather*}
$$ (eq21-chap3)

for all $\vek v_k$, $k=1,\dots N_{u}$ and $q_k$, $k=1,\dots N_{p}$ with

$$
\begin{align*}
\ten A_{km} &= \int_\Omega \mu [(\nao\phi_k\cdot\nao\phi_m)\ten I+\nao\phi_m\nao\phi_k]\,d\Omega, \qquad
k,m=1,\dots N_{u} \\
\vek b_{km}&=-\int_\Omega \psi_k\nao\phi_m\,d\Omega, \qquad k=1,\dots N_{p}, \qquad m=1,\dots N_{u} \\
\vek f_k&=\int_{\Gamma_\text{N}}\phi_k\vek t_\text{N}\,d\Gamma
+\int_\Omega \phi_k\vek f\,d\Omega, \qquad k=1,\dots N_{u}
\end{align*}
$$ (eq22-chap3)

Since Equations {eq}`eq21-chap3` must be fulfilled for any $\vek v_k$ and $q_k$, we obtain the following set of equations

$$
\begin{align*}
  \sum_{m=1}^{N_u} \ten A_{km}\cdot \vek u_m + \sum_{m=1}^{N_p} \vek b_{mk}p_m
  &=\vek f_k, \qquad k=1,\dots N_{u} \\
  \sum_{m=1}^{N_u}\vek b_{km}\cdot \vek u_m&=0, \qquad k=1,\dots N_{p}
\end{align*}
$$ (eq23-chap3)

<a name="Remark:3.12"></a>
```{admonition} Remark 3.12

In order to arrive at a symmetric system matrix we have multiplied the second equation with $-1$.

```

````{exercise}
:label: ex:3.5

Verify that we have obtained $N_u$ vector equations for $N_u$ velocity vector unknowns $\vek u_k$ and
$N_p$ scalar equations for $N_p$ unknowns $p_k$.

````

````{exercise}
:label: ex:3.6

Derive the expressions for $\ten A_{km}$, $\vek b_{km}$ and $\vek f_k$ as given by Equations {eq}`eq22-chap3`.

````

In order to obtain a matrix system that can be easily programmed in an actual computer code, we introduce a Cartesian coordinate system given by the unit vectors $\vek e_i$, $i=1,\dots,d$. The components with respect to chosen coordinate system of the tensors $\ten A_{km}$  are put into $d\times d$ matrices $\mat{A_{km}}$ and the vectors $\vek u_k$, $\vek b_{km}$, $\vek f_k$ into $d$-dimensional column vectors $\col{u_k}$, $\col{b_{km}}$ and $\col{f_k}$, respectively. Next, the following matrices and column vectors are constructed:

$$

\begin{gather*}
   \mat{A}=  
  \begin{pmatrix}
    \mat{A}_{11} & \mat{A}_{12} & \cdots & \mat{A}_{1N_u} \\
    \mat{A}_{21} & \mat{A}_{22} & \cdots & \mat{A}_{2N_u} \\
      \vdots & \vdots & \ddots & \vdots \\
    \mat{A}_{N_u1} & \mat{A}_{N_u2} & \cdots & \mat{A}_{N_uN_u}
  \end{pmatrix},\qquad
     \mat{B}^T=  
  \begin{pmatrix}
    \col{b}_{11} & \col{b}_{12} & \cdots & \col{b}_{1N_p} \\
    \col{b}_{21} & \col{b}_{22} & \cdots & \col{b}_{2N_p} \\
      \vdots & \vdots & \ddots & \vdots \\
    \col{b}_{N_u1} & \col{b}_{N_u2} & \cdots & \col{b}_{N_uN_p}
  \end{pmatrix} \\
       \col{u}=  
  \begin{pmatrix}
    \col{u}_{1} \\
    \col{u}_{2} \\
      \vdots  \\
    \col{u}_{N_u} 
  \end{pmatrix},\qquad
         \col{p}=  
  \begin{pmatrix}
     p_{1} \\
     p_{2} \\
      \vdots  \\
     p_{N_p} 
  \end{pmatrix},\qquad
         \col{f}=  
  \begin{pmatrix}
    \col{f}_{1} \\
    \col{f}_{2} \\
      \vdots  \\
    \col{f}_{N_u} 
  \end{pmatrix}
\end{gather*}

$$ (eq24-chap3)


````{exercise}
:label: ex:3.7

 Show, that $\mat{A}$ has dimensions $n_u\times n_u=dN_u\times dN_u$, $\mat{B}$ has dimensions $n_p\times n_u=N_p\times dN_u$ and the column vectors 
 $\col{u}$ and $\col{f}$ have a dimension of $n_u=dN_u$.

````

````{exercise}
:label: ex:3.8

Show, that the components of the tensors $\ten A_{km}$ and vectors $\vek b_{km}$ can be written as

$$
\begin{align*}
   A^{km}_{ij}&=\vek e_i\cdot\ten A_{km}\cdot\vek e_j=
   \int_\Omega \mu (\phi_{k,\ell}\phi_{m,\ell}\delta_{ij}+\phi_{m,i}\phi_{k,j})\,d\Omega\\
   b^{km}_i&=\vek e_i\cdot\vek b_{km}=\int_\Omega \psi_k\phi_{m,i}\,d\Omega
\end{align*}
$$ (eq25-chap3)

where the Einstein summation convention ($\ell=1,\dots,d$ is a dummy index) and the comma notation for partial derivatives (see Appendix {ref}`Appendix`) has been used.

````

````{exercise}
:label: ex:3.9

Show that matrix $\mat{A}$ is symmetric, i.e. $\mat{A^T}=\mat{A}$.

````

````{exercise}
:label: ex:3.10

Show that the matrix $\mat A$ is positive definite, i.e.

$$
  \col{u^T} \mat{A} \col{u} > 0\quad\text{for any }\col{u}\neq \col{0}
$$ (eq26-chap3)

assuming that appropriate boundary conditions have been imposed to suppress global rigid motion of the fluid in the domain $\Omega$.

````

````{exercise}
:label: ex:3.11

Argue why $-\mat{B^T}$ is called the *discrete gradient* (matrix) and $-\mat{B}$ the *discrete divergence* (matrix).

````

Using the definition of the matrices and column vectors in Equations {eq}`eq24-chap3`, the system of equations Equations {eq}`eq23-chap3` can be written in matrix form as

$$
\begin{pmatrix}
  \mat{A} & \mat{B^T}\\
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
$$ (eq27-chap3)

The system matrix is symmetric and the structure can be further characterized by a positive definite matrix $\mat A$ on the upper part of diagonal and a zero block on the lower part of the diagonal. This is a characteristic feature of a so-called *mixed approximation*, appearing when the weak form is equivalent to a saddle point problem. The structure of the system matrix presents us with a big problem. Considering the last $n_p$ columns of the matrix, i.e.

$$
 \begin{pmatrix}
  \mat{B^T} \\
  \mat{0}
\end{pmatrix}  
$$

we see that for $n_p>n_u$, due to the $n_p$ zero rows, the maximum rank of this matrix ($=$ the number of independent columns) is the number of rows of $\mat B^T$, which is $n_u$. This means that if $n_p>n_u$, the coefficient matrix in Equation {eq}`eq27-chap3` becomes singular! 

In order to analyze the uniqueness of the solutions for $\vek u$ and $p$ further, we consider a system with $\col{f}=\col{0}$:

$$
\begin{align*}
 \mat{A}\col{u}+\mat{B^T}\col{p} &= \col{0} \\
 \mat{B}\col{u} &= \col{0}
\end{align*}
$$ (eq28-chap3)

with homogeneous boundary conditions, i.e. $\vek u_\text{D}=\vek 0$ in Equation {eq}`eq3-chap3` (first equation) and $\vek t_\text{N}=\vek 0$ in {eq}`eq3-chap3` (second equation). The solution should be unique, i.e. $\col{u}=\col{0}$, $\col{p}=\col{0}$.

````{exercise}
:label: ex:3.12

Assume, that both $\{\vek u_1,p_1\}$ and $\{\vek u_2,p_2\}$ satisfy the system Equation {eq}`eq27-chap3`. Argue, that the difference solution $\{\vek u_2-\vek u_1,p_2-p_1\}$ satifies the homogeneous system Equations {eq}`eq28-chap3` with homogeneous boundary conditions, which shows that a zero solution for the homogeneous problem is equivalent to having a unique solution.

````

Multiplying Equation {eq}`eq28-chap3` (second equation) with $\col{p^T}$ gives:

$$
 \col{p^T}\mat{B} \col{u} =0
$$ (eq29-chap3)

which can be used to find after multiplying Equation {eq}`eq28-chap3` (first equation) with $\col{u^T}$:

$$
\begin{split}
  \col{u^T}(\mat{A}\col{u}+\mat{B^T}\col{p})&=\col{u^T}\mat{A}\col{u}+\col{u^T}\mat{B^T}\col{p}\\
  &=\col{u^T}\mat{A}\col{u}+\col{p^T}\mat{B}\col{u}\\
        &=\col{u^T}\mat{A}\col{u}=0
\end{split}
$$ (eq30-chap3)

Since $\mat{A}$ is a positive definite matrix, we find $\col{u}=\col{0}$. Hence, the velocity solution of the discrete Stokes system is unique. This is not necessarily true for the pressure, since substituting $\col u=\col0$ into Equation {eq}`eq28-chap3` (first equation) gives:

$$
   \mat{B^T}\col{p} =\col{0}
$$ (eq31-chap3)

which only leads to $\col{p}=\col{0}$ if the *nullspace* of $\mat{B^T}$ consists of the zero vector $\col{0}$ only. In other words: the columns of the $n_u\times n_p$ matrix $\mat{B^T}$ need to be independent column vectors. Obviously, this can only be achieved if $n_p\leq n_u$, i.e. the number of pressure unknowns cannot be larger than the number of velocity unknowns. Although $n_p\leq n_u$ is a necessary requirement for uniqueness, it is not a sufficient one. In the next section we will see some examples, where the dimension of the nullspace of $\mat B^T$ is not zero, even for $n_p\leq n_u$. This leads to possible pressure solutions $\col{p}\neq\col{0}$ and therefore to non-unique solutions for the pressure. These solutions are called *spurious pressure modes*.

A necessary and sufficient condition for unique and convergent (stable) mixed approximations (such as velocity/pressure in the Stokes equations) is the so-called *inf-sup condition*, also known as the *Babuška-Brezzi* (BB) or *Ladyzhenskaya-Babuška-Brezzi* (LBB) condition. An analysis, whether a particular choice for the finite element spaces for velocity and pressure satisfies the inf-sup condition is rather technical and beyond the scope of these lecture notes. A good entry for further reading and references on this subject is the book of Elman *et al.* {cite}`Elman2014`, in particular the first part of Chapter 3, which has a pragmatic approach to the subject.

<a name="Remark:3.13"></a>
```{admonition} Remark 3.13

If $\Gamma_\text{D}=\Gamma$ (the whole boundary is a Dirichlet boundary), the pressure is only unique up to a constant, if no further conditions for the pressure have been imposed. In that case the dimension of the nullspace is one, having a basis vector $\col{e^T}=(1,1,\dots,1)$ for stable elements. Here, it is assumed that the pressure basis functions have the Kronecker delta property, i.e. $\psi_j(\vek x_i)=\delta_{ij}$.

```

<a name="Remark:3.14"></a>
```{admonition} Remark 3.14

In [Remark 3.5](#Remark:3.5) it is mentioned that the weak form of the Stokes equations is equivalent to a saddle point problem as a result of a constraint (incompressibility) imposed using a Lagrange multiplier. Imposing constraints using Lagrange multipliers generally leads to saddle point problems and thus to mixed approximations. Another common approach leading to saddle point problems is denoting the derivative of a primary unknown as an independent variable in the weak form. This also leads to mixed approximations and the problems associated with it.

```
(Chap3.4)=
## Examples of unstable elements

Although it is not recommended to use unstable elements for practical problems, it is good to present some and 
learn about the problems associated with them. In any case, it is quite revealing that simple elements using just linear and/or constant basis functions are unstable.

A well-known unstable element is the linear-constant pair on a triangle ($P_1$-$P_0$, see {numref}`fig2-chap3`). The velocity has (continuous) linear basis functions $P_1$ and the pressure is a constant per element $P_0$ (see {numref}`Chap2.6` for the notation of basis functions).


```{figure-md} fig2-chap3

<img src="media/chap3/fig2-chap3.png"  width="400px">

$P_1$-$P_0$ element.
```

Consider a 2D flow problem on an $N\times N$ mesh of which an example for $N=6$ is shown in {numref}`fig3-chap3`.

```{figure-md} fig3-chap3

<img src="media/chap3/fig3-chap3.png"  width="500px">

Example of an $N\times N$ mesh of triangles ($N=6$).
```

On the boundary (coloured in blue), homogeneous Dirichlet conditions are imposed, i.e. $\vek u_\text{D}=\vek 0$. The flow is generated by a right-hand side vector $\col{f}\neq\col{0}$.

````{exercise}
:label: ex:3.13

Show, that for an $N\times N$ mesh of $P_1$-$P_0$ triangular elements the number of unknown velocity degrees of freedom $n_u=2(N-1)^2$ and the number of unknown pressure degrees of freedom $n_p=2N^2-1$.
````

From [Exercise 3.13](#ex:3.13) we find that $n_p-n_u=4N-3$. Hence, $n_p>n_u$ for any value of $N$ and the pressure field is not unique, since the dimension of the nullspace of $\mat B^T$ is non-zero. Although the velocity field is unique, it is not accurate either. Actually $\col{u}=\col{0}$, i.e. the flow fully *locks* due to large number of pressure degrees of freedom leading to overconstraining of the incompressibility condition, whatever the value of $\col{f}$. The matrix $\mat{B}$ in the discrete incompressibility condition in Equation {eq}`eq28-chap3` (second equation) has more rows ($n_p$) than columns ($n_u$). There are $n_u$ independent rows leading to the only possible solution $\col{u}=\col{0}$. This can also be shown as follows. Since the pressure test function $q_h$ is defined *per element*, Equation {eq}`eq19-chap3` (second equation) is valid *elementwise*:

$$
\int_{\Omega_\text{e}} q_h\nao\cdot\vek u_h\,d\Omega =0
$$

and since both $q_h$ and $\nao\cdot\vek u_h$ are constant per element, we find that

$$
   \nao\cdot\vek u_h=0, \quad \text{for all }\vek x\in \Omega
$$

i.e. the volume rate-of-change is zero in each point of the domain. In {numref}`fig4-chap3` the $4\times4$ lower-left corner of an $N\times N$ mesh is shown. Starting from the lower-left corner, the velocities on the wall of the two elements in the corner are zero. The velocity in the internal node (common to the two elements) can only be upwards in order for the red element to not change it's volume. Similarly the velocity in that node can only be horizontal in order for the blue element to not change it's volume. Therefore the only possibility is the velocity vector to be zero in that node. The same procedure then can be applied to all elements one-by-one and therefore all nodal velocities in the complete mesh must be zero.

```{figure-md} fig4-chap3

<img src="media/chap3/fig4-chap3.png"  width="700px">

Locking of $P_1$-$P_0$ triangles.
```

Another well-known unstable element is the bilinear-constant pair on a quadrilateral ($Q_1$-$P_0$, see {numref}`fig5-chap3`). The velocity has (continuous) bilinear basis functions $Q_1$ and the pressure is a constant per element $P_0$.

```{figure-md} fig5-chap3

<img src="media/chap3/fig5-chap3.png"  width="400px">

$Q_1$-$P_0$ element.
```

Consider a 2D flow problem on an $N\times N$ mesh of which an example for $N=6$ is shown in {numref}`fig6-chap3`. 

```{figure-md} fig6-chap3

<img src="media/chap3/fig6-chap3.png"  width="500px">

Example of an $N\times N$ mesh of quadrilaterals ($N=6$).
```

On the boundary (coloured in blue), homogeneous Dirichlet conditions are imposed, i.e. $\vek u_\text{D}=\vek 0$. The flow is generated by a right-hand side vector $\col{f}\neq\col{0}$.

````{exercise}
:label: ex:3.14

Show, that for an $N\times N$ mesh of $Q_1$-$P_0$ quadrilateral elements the number of unknown velocity degrees of freedom $n_u=2(N-1)^2$ and the number of unknown pressure degrees of freedom $n_p=N^2-1$.

````

From [Exercise 3.13](#ex:3.13) we find that $n_u-n_p=(N-1)(N-3)$. Hence, $n_u\geq n_p$ for $N\geq3$ and the counting argument for stability is fulfilled, except for a very low number of elements. However, this element suffers from a severe spurious pressure mode taking the form of $+1$ in one element and $-1$ in it's neighbours (see {cite}`Gunzburger1989`, {cite}`Elman2014` and the references therein). This mode is called a *checkerboard mode} for obvious reasons (see {numref}`fig7-chap3`). 

```{figure-md} fig7-chap3

<img src="media/chap3/fig7-chap3.png"  width="500px">

Spurious pressure “checkerboard” mode.
```

For simplicity of implementation and efficiency of computation, some authors prefer “equal-order elements”, i.e. the same (continuous) basis functions for velocity and pressure. On triangles, the most simple element is the $P_1$--$P_1$ element (see {numref}`fig8-chap3`). On quadrilaterals, the most simple element is the $Q_1$--$Q_1$ element (see {numref}`fig9-chap3`). Both of these elements are unstable and suffer from spurious pressure modes. In fact, the nullspace of $\mat B^T$ for a square grid with an even number of $Q_1$--$Q_1$ elements is eight-dimensional {cite}`Elman2014`. In {numref}`Chap3.6` we will present a method to stabilize these elements and therefore use equal-order interpolation in practical problems.

````{exercise}
:label: ex:3.15

Show, that for an $N\times N$ mesh of both $P_1$--$P_1$ triangles (as in {numref}`fig3-chap3`) and $Q_1$--$Q_1$ quadrilaterals (as in {numref}`fig6-chap3`)  $n_u=2(N-1)^2$ and $n_p=(N+1)^2-1$. Also show, that you need $N\geq 6$ to fulfill the counting argument for stability $n_u\geq n_p$.

````

```{figure-md} fig8-chap3

<img src="media/chap3/fig8-chap3.png"  width="400px">

$P_1$-$P_1$ element.
```

```{figure-md} fig9-chap3

<img src="media/chap3/fig9-chap3.png"  width="400px">

$Q_1$-$Q_1$ element.
```

In the next section we will present some examples of stable elements, including the widely used variant based on quadratic basis functions for the velocity.

(Chap3.5)=
## Stable elements

In order to create stable elements the velocity space needs to be enhanced compared to the pressure space. 
We distinguish elements based on continuous pressure and discontinuous pressure approximations. However, we will only consider velocity spaces that are continuous.

(Chap3.5.1)=
## The MINI element

The first stable element is based on the $P_1$-$P_1$ triangular element ({numref}`fig8-chap3`) where the linear velocity field is replaced by an *extended* linear velocity field, see {numref}`fig10-chap3`. The basis functions of the extended linear velocity field consists of the regular linear functions extended with the third-order bubble function $\phi_\text{bubble}=\lambda_1\lambda_2\lambda_3$, where $\lambda_i$, $i=1,2,3$ are the barycentric coordinates (see {numref}`Chap2.6`). 

```{figure-md} fig10-chap3

<img src="media/chap3/fig10-chap3.png"  width="400px">

$P_1^+$-$P_1$ (MINI) element.
```


````{exercise}
:label: ex:3.16

Verify, that $\phi_\text{bubble}=0$ on the boundary of the element.

````

````{exercise}
:label: ex:3.17


Assume, that we have a local numbering of the MINI element as shown below:

```{figure-md} fig1x-chap3

<img src="media/chap3/fig1x-chap3.png"  width="400px">

Local numbering of the MINI element.
```



where the fourth node is positioned at the center of gravity of the triangle.
Show, that if we demand the Kronecker-delta property, the four basis functions (for each spatial direction of the velocity field) are given by:

$$
\begin{alignat*}{2}
\phi_1 &= \lambda_1-9\lambda_1\lambda_2\lambda_3, & \qquad \phi_2 &= \lambda_2-9\lambda_1\lambda_2\lambda_3 \\ 
\phi_3 &= \lambda_3-9\lambda_1\lambda_2\lambda_3, &  \phi_4 &= 27\lambda_1\lambda_2\lambda_3  
\end{alignat*}
$$

````

````{exercise}
:label: ex:3.18

Show that the counting argument for stability $n_u\geq n_p$ is fullfilled for an $N\times N$ mesh of the MINI element (as in {numref}`fig3-chap3`) for any $N$. Show also that it is even valid for a single element having Dirichlet boundary conditions.

````

The extension to three dimensions of the MINI element (using $\phi_\text{bubble}=\lambda_1\lambda_2\lambda_3\lambda_4$) is also a stable element (see {numref}`fig11-chap3`).

```{figure-md} fig11-chap3

<img src="media/chap3/fig11-chap3.png"  width="500px">

Three-dimensional $P_1^+$-$P_1$ (MINI) element.
```

We refer to {cite}`Boffi2013` and the references therein for an analysis of the MINI element.

Although the MINI element has been quite popular, a disadvantage is that the stability comes at the expense of a large number of degrees of freedom for the velocity. Compared to the $P_1$-$P_1$ element, the number of velocity degrees in 2D is roughly three times larger for the same mesh (compare [Exercise 3.15](#ex:3.15) and [Exercise 3.18](#ex:3.18)) without increasing the accuracy. It is even worse in 3D, where the factor is $5--6$. It is possible to eliminate the degrees in the internal node by static condensation on element level (see for example {cite}`Rens2002`). However, this is a rather complicated procedure that cannot be performed at the basis function level only, but involves all the terms in system of equations. Therefore, the recently developed $P_1$-$P_1$ element stabilized by pressure projection (see {numref}`Chap3.6`) looks to be a more attractive approach for low-order elements. 

(Chap3.5.2)=
## Taylor-Hood elements

This family of stable elements is characterized by *continuous* pressures, where the velocity space is one order higher than the pressure. It is named after Taylor \& Hood {cite}`Taylor1973`. The first well-known example is on a triangle, where the velocity space is quadratic ($P_2$) and the pressure space linear ($P_1$) as shown in {numref}`fig12-chap3`.

```{figure-md} fig12-chap3

<img src="media/chap3/fig12-chap3.png"  width="400px">

$P_2$-$P_1$ element (Taylor-Hood).
```

On a quadrilateral, the velocity space is bi-quadratic ($Q_2$) and the pressure space is bi-linear ($Q_1$) as shown in {numref}`fig13-chap3`.

```{figure-md} fig13-chap3

<img src="media/chap3/fig13-chap3.png"  width="400px">

$Q_2$-$Q_1$ element (Taylor-Hood).
```

````{exercise}
:label: ex:3.19

Show, that for an $N\times N$ mesh of both $P_2$--$P_1$ triangles (as in {numref}`fig3-chap3`) and $Q_2$--$Q_1$ quadrilaterals (as in {numref}`fig6-chap3`)  $n_u=2(2N-1)^2$ and $n_p=(N+1)^2-1$. Also show, that you need $N\geq 2$ to fulfill the counting argument for stability $n_u\geq n_p$. 


````


````{exercise}
:label: ex:3.20

In {cite}`Elman2014` it is shown that you need at least a patch of three $P_2$--$P_1$ triangles or two $Q_2$-$Q_1$  quadrilaterals to obtain stability of the velocity-pressure approximation. Show that this coincides with the counting argument for stability $n_u\geq n_p$.

````

The extension of the Taylor-Hood triangle and quadrilateral to three dimensions is straightforward and as expected (see {numref}`fig14-chap3`).



```{figure-md} fig14-chap3

<img src="media/chap3/fig14-chap3.png"  width="500px">

Three-dimensional $P_2$-$P_1$ tetrahedral element (Taylor-Hood).
```


```{figure-md} fig15-chap3

<img src="media/chap3/fig15-chap3.png"  width="500px">

Three-dimensional $Q_2$-$Q_1$ hexahedral element (Taylor-Hood).
```

The Taylor-Hood family can also be extended to higher-order polynomials {cite}`Boffi2013`.

The Taylor-Hood elements are optimally accurate (see {numref}`Chap3.7.1`) with a relatively low number of degrees of freedom. However, the mass balance (incompressibility) is not fulfilled exactly on elementlevel:

$$
  \int_{\Omega_\text{e}} \nao\cdot\vek u_h\,d\Omega=\int_{\Gamma_\text{e}}\vek n\cdot\vek u_h\,d\Gamma\neq 0
$$ (eq32-chap3)

However, fast convergence to zero (quadratically) can be expected.


````{exercise}
:label: ex:3.21

Argue, that the mass balance is fulfilled exactly on the global (mesh) level for Taylor-Hood elements.

````

(Chap3.5.3)=
### Crouzeix-Raviart elements

Some authors prefer elements that fulfill the mass balance exactly on element level. For that, we need elements with *discontinuous* pressure approximations. This family of stable elements, like the Taylor-Hood family, have a velocity space one order higher than the pressure space and
is named after Crouzeix \& Raviart {cite}`Crouzeix1973`.
The first well-known example is on a triangle, where
the velocity space is *extended* quadratic ($P_2^+$) and the pressure space linear, but discontinuous across the element boundary ($P_1^\text{d}$) as shown in {numref}`fig16-chap3`. Note, that using a quadratic velocity $P_2$ is unstable and the extension with a bubble function $\phi_\text{bubble}$, similar to the MINI element, is required for stability. We see that the requirement of exact mass balance on element level for a triangle comes at the expense of roughly 50\% more velocity degrees of freedom and three times (!) more pressure degrees as compared to the Taylor-Hood element.

```{figure-md} fig16-chap3

<img src="media/chap3/fig16-chap3.png"  width="500px">

$P_2^+$-$P_1^\text{d}$ element (Crouzeix-Raviart).
```

On a quadrilateral the velocity space is equal to the Taylor-Hood element, but the pressure is now linear and discontinuous ($P_1^\text{d}$), as shown in {numref}`fig17-chap3`. We see that the requirement of exact mass balance on element level for a triangle comes at the expense of three times more pressure degrees as compared to the Taylor-Hood element. Note, that a bilinear discontinuous pressure $Q_1^\text{d}$ is unstable. 

```{figure-md} fig17-chap3

<img src="media/chap3/fig17-chap3.png"  width="500px">

$Q_2$-$P_1^\text{d}$ element (Crouzeix-Raviart).
```

The extension of the Crouzeix-Raviart triangle to three dimensions (tetrahedron) is possible but leads to a quite complicated element with four additional side bubble functions and one volume bubble function ({cite}`Boffi2013`). The pressure space remains linear $P_1^\text{d}$ though. The additional bubble functions lead to a much larger number of degrees of freedom for the velocity space. The extension of the Crouzeix-Raviart quadrilateral to three dimensions (hexahedron) is straightforward, i.e. the interpolation is also $Q_2P_1^\text{d}$.

The Crouzeix-Raviart family can also be extended to higher-order polynomials {cite}`Boffi2013`.

````{exercise}
:label: ex:3.22

Assume, that we have a local numbering of the $P_2^+$-$P_1^\text{d}$ element as shown below:

```{figure-md} fig2x-chap3

<img src="media/chap3/fig2x-chap3.png"  width="500px">

Local numbering of the $P_2^+$-$P_1^\text{d}$ element.
```

where the seventh node is positioned at the center of gravity of the triangle.
Show, that if we demand the Kronecker-delta property, the seven basis functions (for each spatial direction of the velocity field) are given by:

$$
\begin{alignat*}{2}
 \phi_1 &= \lambda_1(2\lambda_1-1)+3\lambda_1\lambda_2\lambda_3, & \qquad \phi_2 &= 4\lambda_1\lambda_2-12\lambda_1\lambda_2\lambda_3 \\ 
 \phi_3 &= \lambda_2(2\lambda_2-1)+3\lambda_1\lambda_2\lambda_3, & \qquad \phi_4 &= 4\lambda_2\lambda_3-12\lambda_1\lambda_2\lambda_3 \\ 
  \phi_5 &= \lambda_3(2\lambda_3-1)+3\lambda_1\lambda_2\lambda_3, & \qquad \phi_6 &= 4\lambda_3\lambda_1-12\lambda_1\lambda_2\lambda_3 \\ 
  \phi_7 &= 27\lambda_1\lambda_2\lambda_3  
\end{alignat*}
$$

````

````{exercise}
:label: ex:3.23

Show, that a single Crouzeix-Raviart element fullfills the counting argument for stability $n_u\geq n_p$. In fact, a single Crouzeix-Raviart element is indeed stable {cite}`Elman2014`.

````

````{exercise}
:label: ex:3.24

Show, that for the Crouzeix-Raviart elements the following equations are fulfilled exactly on element level:

$$
  \int_{\Omega_\text{e}} \nao\cdot\vek u_h\,d\Omega=0, \qquad\int_{\Omega_\text{e}} x\nao\cdot\vek u_h\,d\Omega=0, \qquad 
  \int_{\Omega_\text{e}} y\nao\cdot\vek u_h\,d\Omega=0
$$


````

<a name="Remark:3.15"></a>
```{admonition} Remark 3.15
:class: note

There is a problem with the definition of the $P_1$ pressure space on non-parallelogram shaped quadrilaterals. The $P_1$ space is “incomplete” on a quadrilateral (the $Q_1$ is complete). The $P_1$ pressure space can be defined on the reference square and mapped onto the real element (mapped $P_1$) or defined as a $P_1$ space on the real element using the global coordinates $x$ and $y$. In {cite}`Arnold2002` it has been shown that for *general meshes* of quadrilaterals, optimal convergence is lost for the mapped $P_1$ space. 
See also {cite}`Boffi2002`, where it is shown that for general meshes
the order of convergence for both velocity and pressure reduces to $\mathcal{O}(h)$ when using the mapped $P_1$ pressure space. It should be noted however, that meshes containing rectangles and/or parallelograms are not affected. Also, meshes obtained after successive refinement of a given mesh and meshes in general that are of nearly parallelogram type, do have optimal convergence (see {cite}`Matthies2001`). So, for practical problems the additional error may be small. 
The convergence problem for mapped $P_1$ spaces can be solved by using *global shape functions* for the pressure {cite}`Boffi2002`:

$$
    \psi_1 = 1, \quad \psi_2=x, \quad \psi_3=y
$$

```


````{exercise}
:label: ex:3.25

How would you define $\psi_1$, $\psi_2$ and $\psi_3$, such that the degrees of freedom for the pressure in an element become $p(\vek x_7)$, $\plderiv{p}{x}(\vek x_7)$ and $\plderiv{p}{y}(\vek x_7)$, where $\vek x_7$ is the position of the seventh node as defined in [Exercise 3.22](#ex:3.22).

````

### Other elements

The stable elements discussed above are the most popular and most useful, but by no means the only ones. There exist non-conforming elements (i.e. $\vek u_h\not\in (H^1(\Omega))^d$), divergence free elements (i.e. $\nao\cdot\vek u_h=0$), “incomplete velocity” (serendipity) elements, but also elements based on vorticity-streamfunction formulations. We refer to the literature for a further study of these elements ({cite}`Elman2014`,{cite}`Boffi2013`,{cite}`Gunzburger1989`).

(Chap3.6)=
### Stabilized elements

Various stabilization techniques exist for velocity/pressure elements that are
not LBB stable (such as equal order elements). The goal of these techniques is to relax the incompressibility constraint such that the spurious pressure modes are suppressed without affecting the theoretical error bounds {cite}`Elman2014`. The effect is that in practical terms the zero block in Equation {eq}`eq27-chap3` becomes filled with non-zeros.

Most of these stabilization techniques are based on least-square terms of the residual added to the regular weak form (Galerkin-Least-Squares (GLS) {cite}`Bochev2009`).
Drawbacks of that approach are:

- The residual contains second-order derivatives.
- If the momentum balance contains other terms, like in Navier-Stokes or viscoelastic flows, these terms need to be taken into account in the residual.
- Tweaking parameters need to be introduced for which it is not straightforward what values need to be taken. Usually, these depend on the relative magnitude of the various terms and the size of the elements.
- The system becomes non-symmetric, even for the symmetric Stokes system.

Dohrmann $\&$ Bochev {cite}`Dohrmann2004` introduced a simple stabilization technique based on local pressure projections (and thus does not use the residual). Only a least-squares pressure term is added to the weak form. The additional term can be computed using standard finite element techniques, adds only to the zero
diagonal block, does not need a tweaking parameter and is symmetric. It only applies to continuous pressure interpolations, such as in equal order velocity-pressure approximations.

For the Dohrmann-Bochev stabilization we need to introduce an additional
pressure finite element space $P'_h$, consisting of piecewise discontinuous
polynomials:

$$
     p_h'=\tsum_i \psi'_i p'_i
$$

in addition to the normal continuous pressure finite element space $P_h$

$$
     p_h=\tsum_k \psi_k p_k
$$

For example a $P_1$-$P_1$ element with a $P_0^\text{d}$-space for $p_h'$, as shown in {numref}`fig18-chap3`.

```{figure-md} fig18-chap3

<img src="media/chap3/fig18-chap3.png"  width="500px">

$P_1$-$P_1$ element with $P_0^\text{d}$ for $p_h'$.
```

Now we consider the projection $\Pi p_h\in P'_h$ of $p_h\in P_h$ on $P'_h$:

$$
    (q'_h,\Pi p_h-p_h) \text{ for all }q_h'\in P'_h
$$

where the $L^2$ inner product has been defined in [Remark 2.8](#Remark:2.8). 
Writing $\Pi p_h= \psi'_i p'_i$ (omitting the $\tsum$ over $i$), we find

$$
    p'_i=m^{-1}_{ij}(\psi_j',p_h)
$$

and thus

$$
   \Pi p_h = \psi'_i m^{-1}_{ij}(\psi_j',p_h)=
                \psi'_i m^{-1}_{ij}(\psi_j',\psi_k)p_k
$$

where $m_{ij}=(\psi'_i,\psi'_j)$ is the pressure mass matrix of $P'_h$. Note,
that $\psi'_i$ is discontinuous across element boundaries and the inverse
$m^{-1}_{ij}$ matrix can be computed on element level.

The Stokes problem is determined by the stationary value of the functional
(saddle point problem)

$$
   L(\vek u,p) = \frac12\big(\ten D(\vek u),\mu\ten D(\vek u)\big)-
     (\nabla\cdot\vek u,p) -(\vek u,\vek f)
$$

The Dohrmann-Bochev stabilization consists of subtracting a least-squares term
from the functional $L(\vek u,p)$:

$$
  L'(\vek u,p) = \frac12\big(\ten D(\vek u),\mu\ten D(\vek u)\big)-
     (\nabla\cdot\vek u,p) -(\vek u,\vek f)-
     \frac1{2\mu'}(p-\Pi p,p-\Pi p)
$$

where $\mu'$ is a viscosity parameter for which $\mu'=\mu$ seems to work
just fine. The least-squares term leads to an additional term in the weak
 form of

$$
     \dots -\frac1{\mu'}(q_h-\Pi q_h,p_h-\Pi p_h) = 0
$$

which adds to the continuity equation and fills the ``zero diagonal block''.
Note, that 

$$
    (\Pi q_h,p_h-\Pi p_h)=0 \quad \text{ or } \quad(\Pi q_h,\Pi p_h)=(\Pi q_h,p_h)
$$

because the difference with respect to the projection should be orthogonal to
the projection itself. A different proof is by substituting the expressions for
$\Pi p_h$ and $m_{ij}$:

$$
\begin{multline*}
  (\Pi q_h,\Pi p_h)=
   (\psi'_i m^{-1}_{ij}(\psi_j',q_h), \psi'_k m^{-1}_{km}(\psi_m',p_h))=\\
   (\psi'_i,\psi'_k) m^{-1}_{ij}(\psi_j',q_h)m^{-1}_{km}(\psi_m',p_h)=
        m_{ik} m^{-1}_{ij}(\psi_j',q_h)m^{-1}_{km}(\psi_m',p_h)=\\
        \delta_{kj}(\psi_j',q_h)m^{-1}_{km}(\psi_m',p_h)=
                   (\psi_k',q_h)m^{-1}_{km}(\psi_m',p_h)=(\Pi q_h,p_h)
\end{multline*}
$$

Finally the additional term becomes:

$$
   -\frac1{\mu'}[(q_h,p_h)-(q_h,\Pi p_h)]
$$

with additional matrix

$$
  -\frac1{\mu'}[(\psi_k,\psi_m)-(\psi_k,\psi'_i)m^{-1}_{ij}(\psi'_j,\psi_m)]
$$

which is *symmetric* and negative definite. The structure of the system
matrix becomes

$$
 \begin{pmatrix}
   \mat{A} &\mat{B}^T \\[1ex]
   \mat{B} &-\mat C
 \end{pmatrix}
 \begin{pmatrix}
    \col{u}\\[1ex]
    \col p
 \end{pmatrix}
$$

with $\mat C$ symmetric positive definite.

Some concluding remarks:

- For equal order elements the polynomial order for $\psi'$ is chosen one order less than $\psi$. 
Even lower order is possible (like a constant for quadratic elements), but then accuracy is lost.
- The $P_1/P_1$ pair (see {numref}`fig18-chap3`) with pressure projection stabilization is an attractive alternative for the mini element
$(P_1^+/P_1)$, with less velocity degrees of freedom and constant gradients
inside the element.
- In {cite}`Dohrmann2004` it is mentioned that the stabilization might be
beneficial for iterative solvers, but this has not been tested.

(Chap3.7)=
## Accuracy

In this section we test the accuracy of the elements for the Stokes equations using an analytical test problem.

(Chap3.7.1)=
### Analytical test problem

For testing accuracy we use the following analytical solution of a Stokes problem ($\mu=1$) on a square domain $\Omega=[-1,1]\times[-1,1]$:

$$
 u_x=20xy^3,\quad u_y=5x^4-5y^4, \quad p=60x^2y-10y^3+\text{constant}
$$ (eq33-chap3)

as derived from the streamfunction $\psi(x,y)=5xy^4-x^5$ {cite}`Elman2014`. The contours of streamfunction $\psi$ are shown in {numref}`fig19-chap3`. A plot of velocity vector field $\vek u$ is shown in {numref}`fig20-chap3` and contours of the pressure field $p$ in {numref}`fig21-chap3`. Note, that the maximum velocity magnitude is $20$ and the maximum pressure difference in the domain is $80$.


```{figure-md} fig19-chap3

<img src="media/chap3/fig19-chap3.png"  width="500px">

Streamfunction contours of the analytical solution Equation {eq}`eq33-chap3`.
```

```{figure-md} fig20-chap3

<img src="media/chap3/fig20-chap3.png"  width="500px">

Velocity vectors (scaled) of the analytical solution Equation {eq}`eq33-chap3`.
```

```{figure-md} fig21-chap3

<img src="media/chap3/fig21-chap3.png"  width="500px">

Pressure contours of the analytical solution Equation {eq}`eq33-chap3` ($\text{constant}=0$).
```

We solve the problem numerically using a mesh of $N\times N$ elements (total of $2N^2$ triangles or $N^2$ quadrilaterals), similar to the {numref}`fig3-chap3` and {numref}`fig6-chap3`. The element size $h=2/N$. On the boundary we impose Dirichlet boundary conditions (exact solution in the nodes) and we impose the pressure level by imposing zero pressure in a single node of the mesh.

For defining the interpolation error and the error of the numerical solution we use the area-scaled $L^2(\Omega)$-norm for scalars and vectors:

$$
 \|a\|=\sqrt{\frac1A\int_\Omega a^2\,d\Omega},\qquad\|\vek a\|=\sqrt{\frac1A\int_\Omega |\vek a|^2\,d\Omega}
$$ (eq34-chap3)

where $A=\int_\Omega d\Omega$. 
The errors for velocity and pressure finite element solution become:

$$
 e^h_u=\|\vek u_h -\vek u\|, \qquad e^h_p=\| p_h-\frac1A\int_\Omega p_h\,d\Omega -p \|
$$ (eq35-chap3)

where the exact solution for $p$ should have zero level, i.e. $\int_\Omega p\,d\Omega=0$.

In the following we define $k$ and $m$ to be the (full) polynomial degree of the velocity and pressure field interpolations, i.e. the degree is $0$ for $P_0$, $1$ for $P_1$, $P_1^+$ and $Q_1$, $2$ for $P_2$, $P_2^+$ and $Q_2$. Whether the interpolation is continuous or discontinuous does not affect the degree of interpolation. 
Defining $\hat{\vek u}_h$ and $\hat p_h$ as the discrete fields assuming exact values in the nodes, for 
 sufficiently smooth solutions for $\vek u$ and $p$ the interpolation error 
can then be written as

$$
\begin{align*}
  \|\hat{\vek u}_h-\vek u\| &\leq \hat C_u h^{k+1}\\
  \|\hat p_h-p\| &\leq \hat C_p h^{m+1}
\end{align*}
$$ (eq36-chap3)


where $\hat C_u$ and $\hat C_p$ are constants depending on the solution $\vek u$ and $p$, but not on the size of the elements $h$. If the errors of the finite element solutions $e_u^h$ and $e_p^h$, as defined in Equation {eq}`eq35-chap3`, are
converging at the same rates as the interpolation errors in Equations {eq}`eq36-chap3`, the velocity-pressure element is called *optimally convergent*.
In the theory of errors (see for example the book of Elman *et al.* {cite}`Elman2014`), the theoretical error estimate is a sum of the error in $\nao\vek u$ and $p$. The interpolation error of $\nao\vek u$ is $\mathcal{O}(h^k)$.
Therefore, elements which are “balanced” regarding the errors usually have $m=k-1$. This leads to an error for the pressure of $\mathcal{O}(h^k)$. In theory, increasing $m$ would not improve the error bound (but see the first order elements below). Taking $m=k-2$ or lower affects the error of the pressure, but also “pollutes” the velocity and reduces the convergence rate of the velocity as well (see the $Q_2P_0$ element below). 

(Chap3.2.7)=
### Unstable elements

The unstable $P_1P_0$ element (see {numref}`Chap3.4`) does not give accurate solutions for velocity and pressure (it is even singular in the latter). However, the $Q_1Q_0$ element does give accurate solutions for the velocity field, at least for the test problem at hand. It is even slightly more accurate than the stabilized $Q_1Q_1$ element (see {numref}`Chap3.6`). The pressure field suffers from a spurious mode, though. If the solver does not detect the singularity, the obtained pressure field looks like an art painting  (see {numref}`fig22-chap3`) and is obviously wrong. This has inspired some to filter out the single spurious pressure mode to obtain a smooth pressure field. As noted by Gunzburger \cite{Gunzburger1989}, this a dangerous adventure with a possibility of poor pressure approximations and it is advised to stay away from this element.   

```{figure-md} fig22-chap3

<img src="media/chap3/fig22-chap3.png"  width="500px">

Pressure field of the test problem for 32$\times$32 $Q_1Q_0$ elements.
```

(Chap3.7.3)=
### The MINI element

The MINI element ($P_1^+$-$P_1$) has polynomial degrees of $k=1$ and $m=1$. Optimal convergence would mean $e_u^h=\mathcal{O}(h^2)$ and
$e_p^h=\mathcal{O}(h^2)$. From the theory of errors we would expect an error bound $\mathcal{O}(h)$ for the pressure though. In {numref}`fig23-chap3` we show the errors as obtained in the test problem. The errors are shown as a function of $N$. Since $h=2/N$, the order of convergence is *minus* the slope of the lines. We clearly see that $e_u^h$ is according to the theory, but the slope of the pressure is better than the theory would predict (slope around 1.7).

```{figure-md} fig23-chap3

<img src="media/chap3/fig23-chap3.png"  width="500px">

Error norms $e_u^h$ and $e_p^h$ for the MINI element as a function of $N$.
```

(Chap3.7.4)=
### Taylor-Hood elements

The Taylor-Hood elements ($P_2$-$P_1$ and $Q_2$-$Q_1$) have polynomial degrees of $k=2$ and $m=1$. Optimal convergence would mean $e_u^h=\mathcal{O}(h^3)$ and
$e_p^h=\mathcal{O}(h^2)$. From the theory of errors we expect the same rate of convergence. In {numref}`fig24-chap3` we show the errors as obtained in the test problem. We clearly see that both $e_u^h$ and
$e_p^h$ are according to the theory. The errors of the triangular and quadrilateral elements are similar.

```{figure-md} fig24-chap3

<img src="media/chap3/fig24-chap3.png"  width="500px">

Error norms $e_u^h$ and $e_p^h$ for the Taylor-Hood elements as a function of $N$.
```

(Chap3.7.5)=
### Crouzeiz-Raviart elements

The Crouzeix-Raviart elements ($P^+_2$-$P_1^\text{d}$ and $Q_2$-$P_1^\text{d}$) have polynomial degrees of $k=2$ and $m=1$. Optimal convergence would mean $e_u^h=\mathcal{O}(h^3)$ and
$e_p^h=\mathcal{O}(h^2)$. From the theory of errors we expect the same rate of convergence. In {numref}`fig25-chap3` we show the errors as obtained in the test problem. We clearly see that both $e_u^h$ and
$e_p^h$ are according to the theory. The errors of the quadrilateral elements are somewhat lower than the
errors for the triangular elements. This is surprising, considering the larger number of the degrees of freedom in the mesh for the triangular elements.

```{figure-md} fig25-chap3

<img src="media/chap3/fig25-chap3.png"  width="500px">

Error norms $e_u^h$ and $e_p^h$ for the Crouzeix-Raviart elements as a function of $N$.
```

Replacing the $P_1^\text{d}$ pressure interpolation in the quadrilateral $Q_2$-$P_1^\text{d}$ with $P_0$ pressure (see {numref}`fig26-chap3`) 

```{figure-md} fig26-chap3

<img src="media/chap3/fig26-chap3.png"  width="500px">

$Q_2$-$P_0^\text{d}$ element.
```

interpolation leads to a stable element. However, the reduction of the order of the pressure space also lowers the order of convergence of the velocity space by one (see {numref}`fig27-chap3`) and is considered “false economy” ({cite}`Elman2014`).

```{figure-md} fig27-chap3

<img src="media/chap3/fig27-chap3.png"  width="500px">

Error norms $e_u^h$ and $e_p^h$ for the $Q_2P_0$ element as compared to the $Q_2P_1^\text{d}$ Crouzeix-Raviart elements as a function of $N$.
```

(Chap3.7.6)=
### Stabilized elements

The pressure stabilized elements ($P_1$-$P_1$ and $Q_1$-$Q_1$) have polynomial degrees of $k=1$ and $m=1$. Optimal convergence would mean $e_u^h=\mathcal{O}(h^2)$ and
$e_p^h=\mathcal{O}(h^2)$. However, the pressure space is projected on a space that is one order lower ($P_0$). Therefore, from the theory of errors it is expected that the rate of convergence for the pressure solution $e_p^h=\mathcal{O}(h)$ {cite}`Dohrmann2004`. In {numref}`fig28-chap3` we show the errors as obtained in the test problem. We also show again the results of the MINI element.
We clearly see that $e_u^h$ is according to the theory, but the slope of the pressure error is better than the theory would predict (slope around 1.7). We see also that the convergence *rate* of the stabilized elements is very close to the MINI element. The pressure solution seems to be somewhat more accurate for the stabilized elements, though (more than a factor of two for the test problem). This is consistent with the results in {cite}`Dohrmann2004`.
In {cite}`Dohrmann2004` it is also found that the “superconvergence” of the pressure is only presence in the lowest order elements ($P_1$-$P_1$ and $Q_1$-$Q_1$). The convergence rate for the pressure of the higher order elements (such as $P_2$-$P_2$ and $Q_2$-$Q_2$ and beyond), is according to the theory.


```{figure-md} fig28-chap3

<img src="media/chap3/fig28-chap3.png"  width="500px">

Error norms $e_u^h$ and $e_p^h$ for the elements stabilized by pressure projection as a function of $N$.
```

