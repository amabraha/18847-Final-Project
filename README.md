# 18847 Final Project
 
 #### Alan Abraham
 #### Myles Mwathe
 #### Yao Xiao
 
 
 1. Integrate Triangle into FEGrid
- To test, run `make poly` in `exec/`.

 2. Template FEPoissonOperator by data type
- To test, run `make run_float`, `make run_double`, `make run_complex_float`, or `make run_complex_double` in `exec/`

 3. Non-zero Dirichlet boundary conditions
- To test, run `make run2d` in `exec/`

 4. Demonstrate Piecewise Linear elements converge at 2nd order accuracy
- To test, run `make run2d` in `exec/`
 5. Implement a time-dependent FEM solver
- To test, run `make run_TD` in `exec/`
 6. Verify correct time-dependent behavior
- To test, run `make run_TD` in `exec/`
 7. Create time-dependent animations of interesting source terms and boundary conditions
- To test, run `make run_TD` in `exec/`
 8. Extrude the 2D elements into an extrusion into 3D
- To test, run `make clean` and then `make run3d DIM=3`
 9. Solve Poisson Equation in 3D
- To test, run `make clean` and then `make run3d DIM=3`
 10. Verify 3D steady solutions
- To test, run `make clean` and then `make run3d DIM=3`