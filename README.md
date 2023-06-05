# C++ finite elements

Project assignment for the course C und C++ in der Simulationsentwicklung.

Group members:

Martin Fasser and Charlotte Vavourakis

## Instructions

* Create excecutables and run experiment (first dense implemention only)

```
mdkir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release #extra flag to avoid debug mode for the google benchmarking (significantly slower)
cmake --build .
./PDEsolver # run the PDE solver
./PDE-dense-bench # run the benchmark
```

## Goal

Implement a solver for the heat equation.

General formula: $\frac{du}{dt} = \nabla u$ (special case, thermal diffusivity = 1)

For domain $\Omega$ in two dimensions, the Laplacian can be rewritten: $\frac{du}{dt} = \frac{d^2u}{dx^2} + \frac{d^2u}{dy^2}$

## Solution

* The problem is rewritten to a weak formulation by multiplying with a test function $\phi : \Omega \subseteq \R^2 \to \R$ and integrating over domain $\Omega$.

* To pick $\phi$, $\Omega$ is discretisized using triangles. Then need to find u so that weak formulation is satisfied for all $\phi(D(\Omega))$.
Given is a mesh file generated with Gmsh (square.msh), will parse this .vtk file so that can be plotted with ParaView. Code given in read-mesh-template.cpp.

	* divide $\Omega$ in $n_T$ triangles with 3-tuple vertices $v_1,..., v_{n_{v}}$
	* for $(v_i, v_j, v_k) : \phi_i(x,y) = a + bx + cy$ and $\phi_i = 0$ inside each triangle that does not contain $v_i$

* Final problem rewritten in matrix notation, with $U = [u_1, ...u_{n_{v}}]$ : 

$\delta_tU = BU$ 

* This can then be solved using forward Euler:

$U^{n+1} = U^n + dt B U^n$ for a given dt, $U(t_0,x,y)$ and $\Omega$. 
