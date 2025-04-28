# 18847 Final Project

#### Alan Abraham
#### Myles Mwathe
#### Yao Xiao


1. Integrate Triangle into FEGrid (Yao)
2. Template FEPoissonOperator by data type
3. Non-zero Dirichlet boundary conditions
4. Demonstrate Piecewise Linear elements converge at 2nd order accuracy
5. Implement a time-dependent FEM solver (Alan)
6. Verify correct time-dependent behavior
7. Create time-dependent animations of interesting source terms and boundary conditions

For time dependent, implemented a new class in `src/fe_td/FETimeDependent.H`. This class will approximate the solution to `dphi/dt = f - Lphi` by using backward euler. It calls the step function for each iteration of the time integration method.
The testcase is in `exec/FEMain_TD.cpp`. To run it call `make run_TD` in `exec/`. It will create a `.vtk` database in `exec/vtk_output/` with the approximated solution. It also creates a `.vtk` database in `exec/phi_output/` with the source data used to derive the RHS `f`. These should be roughly equal and the error is printed at the end of the test.

8. Extrude the 2D elements into an extrusion into 3D
9. Solve Poisson Equation in 3D
10. Verify 3D steady solutions
