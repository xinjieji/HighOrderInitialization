#  High Order Reinitialization

This package mainly implements the 2D reinitialization method proposed by [Saye(2014)](dx.doi.org/10.2140/camcos.2014.9.107). It also wraps some useful functions for       implicit surface representation and the evolution calculation.
  
  
## About implicit surfaces and reinitialization
Implicit surfaces are a powerful technique for representing complex geometries without the need for explicit parameterization. An implicit surface is usually defined by a level set of a continuous function $\phi: \mathbb{R}^d \rightarrow \mathbb{R}$, where $\phi$ maps a point in the $d$- dimension Euclidean space to a scalar value, with the curve/surface defined as the set of points where $\phi(\mathbf{x}) = 0$. Implicit surfaces are able to represent topologies without any extra functions, and they guarantee the existence of solutions in the class of viscosity partial differential equations. As a result, implicit surfaces bring mathematical and computational efficiency in a variety of problems.

In the applications, one of the main challenges with implicit surfaces is that their level sets, which define the boundaries of the surface, can become distorted or irregular over time due to numerical errors or other factors. This can lead to inaccurate surface representations and can cause problems for applications such as collision detection and fluid simulation. Reinitialization/Redistancing is a process that is used to correct these distortions by redefining the level set function of an implicit surface to be a signed distance function (SDF). An SDF is a scalar field that assigns to each point in space the distance to the nearest point on the surface of the object, with positive values inside the object and negative values outside. If the reinitialization method is applied as a pre-processing of the simulations, we could preserve the stability and accuracy of the implicit surfaces.

## Tutorial
