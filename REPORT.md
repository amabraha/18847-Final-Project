# 18847 Final Project - Finite Element Method

## Authors

Alan Abraham (<amabraha@andrew.cmu.edu>)

Myles Mwathe (<mmwathe@andrew.cmu.edu>)

Yao Xiao (<yaox3@andrew.cmu.edu>)

# 1. Integrate Triangle into FEGrid (Yao)

- Scan .poly file

- Call `triangulate()` from Triangle library (written in **C**) to refine triangulation with specified area constraint

  ```c
  extern “C”
  ```

- Load nodes and elements from output

  - same as original constructor

## Sample VisIt output of .poly file

![A red and yellow letter  Description automatically generated](https://lh7-rt.googleusercontent.com/slidesz/AGV_vUcNkZAnuMOtObEco4b61-kyEHwnTv1bQqS6ZCCziejpdcjy8MUX3wJepjv3CFOt5fUKREzTCGdmbTtJr_c934v6xwJ77mYSKS8o7hX6YDNsshuNxnS_c_vj4BSSl6Cl9PlSFaAkWw=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe)



# 2. Template FEPoissonOperator by data type (Myles)

## Key changes required

- Converted classes to template \<typename T>

- Modified member variables and method signatures to use type T

- Added explicit template instantiations for each supported type

- Updated constructors and operators to handle templated types

## Modified core classes

- `FEPoissonOperator`
- `SparseMatrix` 
- `JacobiSolver`

## Implementation approach

- Separated declarations (.H) from implementations (.hpp)

- Added necessary header includes (e.g., \<complex>)

## Testing

- Created test cases for each data type

- Verified solution correctness for annulus mesh

- Ensured consistent behavior across all supported types

# 3. Non-zero Dirichlet boundary conditions (Myles)

## **Goal**

Set solution directly at boundaries to known values:

- **Before:** Boundary assumed to be zero.
- **Now:** Allows custom boundary values $\phi_\Phi$.

## Implementation

- Solve on *all* nodes (interior + boundary).
- Adjust matrix after assembly for boundary nodes:
  - Set diagonal = 1, other entries = 0.
  - Set RHS = known boundary values.

## Effect

- System remains solvable (positive definite, diagonally dominant, not symmetric).

## Sample .poly file with boundary conditions

![img](https://lh7-rt.googleusercontent.com/slidesz/AGV_vUdZekRxVhPDjEQeo4F-vuhwTqQlzQ0S3JMCYIOA28XXBHO3OFP65z1fSx0NeiJuNM-QGyr73XGGZAG1_TgMiiHbEWTAzVnMI4NoTlrfg9fXvmUSuQ7CYzP2yx_ueGLjBxLixCjwRg=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe)

# 4. Demonstrate Piecewise Linear elements converge at 2nd order accuracy

## **Setup**

- Pick smooth $\Phi$, build RHS from $-\Delta\Phi$.
- Enforce $\Phi$ on all boundary nodes.

## Exact Solution

```c++
static auto Phi_exact = [](const array<double, DIM> &X) -> double
{
    // smooth function that vanishes on the unit square boundary
    return X[0] * X[0] + X[1] * X[1];
};

static auto source2D = [](const array<double, DIM> &X) -> double
{
    return -4.0;
};
```

*Notes: Why still quadratic? 2D Linear function converges too fast to reveal the order.*

## **Workflow**

Refine mesh → smaller max element area.

1. Assemble stiffness matrix & load vector.
2. Stamp Dirichlet rows.
3. Solve for $\Phi_h$.
4. Measure nodal max-error.

![img](https://lh7-rt.googleusercontent.com/slidesz/AGV_vUdHT9WVviBnjy8-qEY9TIi_lGYL0i1RViE0kezg6nfS43SMs56nxdJluLY8Wa18SN_fsl8dcmGmW0YFsCcW6qebbYnHX2ZfoTLPjkdec42Uc2PA5_yFL0o9GSqGDJYWoMIKtUL4=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe)

## Results

```
Estimated order p:
between lvl 0->1 : p ≈ 2
between lvl 1->2 : p ≈ 2
between lvl 2->3 : p ≈ 0.265492
between lvl 3->4 : p ≈ 2.68471
```

# 5. Implement a time-dependent FEM solver (Alan)

Solve the differential equation

![img](https://lh7-rt.googleusercontent.com/slidesz/AGV_vUde4aIm1-oRGz73PTya1Byvaa_r34GmZJtHPVN0meSoLo9_s1GWhA0odmEZD8uwEj3ttxd4qbdW0TaRpTIJnm_aLbuFQ8Juh-NLj-bQAx6Ss4XrJp0NvM1fnFcLKLJN1TWB4KrN4Q=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe)

Via backwards euler, so we get the recurrence

![img](https://lh7-rt.googleusercontent.com/slidesz/AGV_vUchiAOUYPsk9WsypSIATFSvKSNUyEME6SixpP5yn-z_85zKKzWbN1xqwlJVfVJwWEzkvdHifkgdnsb8edYDSaUF2LeTp3on5f-e7m-FljmnHoAYmUTQmVFXb_Fsnh9lf1ZDVUNirA=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe)

$L$ is the finite element approximation of the laplacian (with boundary conditions)

In the next sections, to run the time dependent FEM solver, run `make run_TD` in `exec/`. This will create a database of indexed visit files. The computed solution will be in `exec/vtk_output/solutionXXX.vtk` where `XXX` is the time index of the file. The reference $\Phi$ function that we use will be outputted in `exec/phi_output/solutionXXX.vtk`.

# 6. Verification of Time Dependent FEM Solver

## General idea

- Pick a $\Phi$

- Let

  ![img](https://lh7-rt.googleusercontent.com/slidesz/AGV_vUcGxjLzBg9QUia3906UZ2q_XvkTBpQ1wwCyf6fKTfWUF4jS6FyDwoC67oVU8eA3e_0eCqYKe7hZg6VJIvNGe0IQ9XCAEbbUn8a7CrPM3mSDwPl6N4GzLl41d74q_aSj8eJaXMxe9Q=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe)

The theoretical solution to the differential equation should converge to the $\Phi$ we picked regardless of the initial conditions (if we ignore boundary conditions).

*Notes: There is inherent error with a numerical time integrator chasing after the theoretical solution.*

## Example1: Time independent $\Phi$

Final relative error **0.009**

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUf82zJlupZxmJ7qQ3ARAlApiJrf_7GTMGjV1FJHxM_WhwR-DrXtIjvVxjnnx97skOaKK8CGo2xnL1JuyUNSRKFoVSdwTDq7ftSmLAyii384cjVTAMrmrSopC6MyOUp6yOhgPBpZdw=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="img" style="zoom:67%;" />

### Computed solution

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUfYBzsYOJaaEjqD7LET_nN3cdB49CPMv0aedlENSo0wCGb4teGyOfOmCJpkuY0dch1Q-tzyNKnP4pmfCHdxCgh78YjluZZjGhTr520iKv7N05yNCFcxzlXaPdRbjtCWKJ0rWUCICg=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="img" style="zoom: 33%;" />



### Reference $\Phi$

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUds0QTJmjFJglo3jbwEWs_15iorc40i3qOMGlDKFU-dRlzRZNdDNQhFTLfmWKBIbJfjuj6J5fACKaR2lDz_YKKa0OGHuWf-w2Xai4DJYZYPXhpDtUSbevtYvvDBVlZzR6noDonO=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="img" style="zoom: 33%;" />

## Example2: Simple time dependent $\Phi$

Final relative error **0.1385**

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUcwsybEUFnNctyasqTbEhbxO1PGgdRJp5jJHTHEBRQnz0tq--PZi2zaXNHYr6GvXRFfqp-t3bRFcO0cyuofOkXP4EUDojrcaGT1sBTyv3c1BWSFSFrj51fNdAp8qPSsNu44KrBMmQ=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="img" style="zoom:67%;" />

### Computed solution

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUdmST2GKgVqCjTgVeWxI3gMmDfLocOCqeaDoRcj7tPdx_M-5XGwNxZhvK06qI4lf4gvFO4x1lqoz8SMJBIjsckjeBZVufSXT0ZKV5BhLGrVwX2_oGZD4fij8isYH5Ta4OJq0NkJSQ=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="img" style="zoom: 50%;" />

### Reference $\Phi$

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUfY8RRSd576iZCrPfGpEeAW-IUnNmm5XIFbbKfRtQcDrI6NlyEpgHPhqtkwmRh38tcjSpSS0gFECkCaFRCo7BZwGgQgT_Hw_r32OffX839DsnECyqNQxXgdejLTbT-_U1dyEVpoDw=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="img" style="zoom: 50%;" />

## Example3: More complex time dependent $\Phi$

Final relative error **0.1057**

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUcqhFZWKnvysyhKvmUqOFnyMtkPh1l4c1UF9rp2M-jnV8nNBqKcmmRGnJNkbZadrAmWHtWnzNt5sm7dxom1d_jn262cNKn9wS8fBoSKU518Cm7BogQgzAV6vod8UyGd3vT-5otjOA=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="img" style="zoom: 67%;" />

### Computed solution

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUdTUwBYANthIlp-lixqOsPAkF-Pn2IwiqCjcEaEnkNN-Ny9nR2fq61dZIp0YTrMnQI0XFAvC79fP1Brx1dIlYamqgEuGUfRrZOZkWqBzGiBkyyuiD_3DrKoN7ANnCymOQdgfA4gUw=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="img" style="zoom: 50%;" />

### Reference $\Phi$

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUe95jeM1N0aw9FBCBhkrviJZSHoaFv6ZbUribQLpyrpiAQxdJEvn_azqfRz6NYdL5fTZrXv05z2a5u-wdiraBnb-dxJjUeXCKc6RkkrA3nQ-ckyp1fdzP-2TgamOeZcuiJvuIm2ag=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="img" style="zoom: 50%;" />

# 7. Time-dependent animations

To create a time-dependent animation, run `make run_TD` in `exec/`. Then open visit and open the database in `exec/vtk_output/solutionXXX.vtk`. Then add the mesh and pseudocolor you want to see and use the next arrow to play the animation. For example, with the last example of $\Phi(x,y,t)=2x+y+3\sin t$, some snapshots are shown below:

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUdmAzJEj64y4xI_6E2aAug4azUTAYbfCCu_oYi_dTxg6TvPL4vd99Eh4iM67GU-siS9rCKIif8xt584n0whtFuXF_F-txGalpeRTikJileM4UI_sJN2EkP73WtyCV1Lav5kn1Tp=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="img" style="zoom:50%;" />

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUfff7OaoPb7h-1s0k8Hd57VdMTCLlSOsWkk66W68XYDx2fJgCOdIgtNArQW7zvVkb7OIezK-rGMltRw1W_23K8stpy40ndJUW0bHySyKGZ1XRcoQt8b0xUmYpOyHgZLpy50VtH0YQ=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="img" style="zoom:50%;" />

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUf3GHd91sADeLMc3VaJ9ugOIMw35qjN80TTiYIKWKM5jVo3Sv2KDOIGLk4Vn9pd5Lp6eui3X_8xEYiI9Z1M2g7tpG8SOa-1FsKDCkuS18nJirTqdK1j8C74cC5qfFyLEmG7h5SCHQ=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="img" style="zoom:50%;" />



# 8. Extrude the 2D elements into 3D (Yao)

- Extrude each triangle vertically in a prism
- Split prism into 3 tetrahedrons

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUctC_8DSBZxSQAvA94sQJTwt9j_XjEqXfEIQvCJ-jA-X_B2oubIgCYrLD1_0C1LPd2qgcWJSuTZ6SkY02u2YUjN0XVDJFryxfMgVWqgj8eUL7Pdr_lKJhTWyzsCDMR9accf2UEz=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="Prisms are split into three tetrahedra. Prisms can be split in six ways, depending on the direction of the triangulation of the quadrilateral faces. The six splits can be characterized with labels (FRR, RFR, RRF, RFF, FRF, and FFR) depending on the triangulation of the quadrilateral faces." style="zoom:50%;" />



# 9. Solve Poisson Equation in 3D (Yao)

2 key differences of another dimension

## **FEGrid::gradient()**

Solve linear system in 3D (inverse matrix: compute determinant, cofactors)

We have

$$
dx \times \Delta\phi= \begin{pmatrix}
 -1\\
 -1\\
 -1
\end{pmatrix}
$$

Then solve

$$
\Delta\phi= (dx)^{-1}\times \begin{pmatrix}
 -1\\
 -1\\
 -1
\end{pmatrix}
$$

## **FEGrid::elementArea()**

Compute tetrahedron volume instead (value: determinant / 6).

## Example extrusion plot of square mesh

<img src="https://lh7-rt.googleusercontent.com/slidesz/AGV_vUfsgB9V-8wQ5OfCXuRvJjnU3iUcjCO60t_OT8c6j1KYd39dlWon3BWr2XLC2qwFdBuIXLY0y_cRo7i1BypBFvfAL6dUgsgQ90IC7aw_lp3wepLqXsOD0hHFKBBiWrM7j7hOw4tm=s2048?key=BrZTxrthEqqZHSlVwKRtJ_Oe" alt="img" style="zoom:50%;" />

# 10. Verify 3D steady solutions

The process closely resembles the 2D version.

## Exact Solution

```c++
static auto Phi_exact = [](const array<double, DIM> &X) -> double
{
    // smooth function that vanishes on the unit square boundary
    return X[0] + X[1] + X[2];
};

static auto source3D = [](const array<double, DIM> &X) -> double
{
    return 0.0;
};
```

## Results

```
Estimated order p:
between lvl 0→1 : p ≈ -0.0939246
between lvl 1→2 : p ≈ 1.39756
between lvl 2→3 : p ≈ 0.307041
between lvl 3→4 : p ≈ 0.413396
```

# References

Dompierre, Julien & Labbé, Paul & Vallet, Marie-Gabrielle & Camarero, Ricardo. (1999). How to Subdivide Pyramids, Prisms, and Hexahedra into Tetrahedra.. 195-204. 