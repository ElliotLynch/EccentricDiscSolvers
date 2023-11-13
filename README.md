# EccentricDiscSolvers
Python Code For Studying Eccentric Waves in 3D fluid Discs

Theory Covered in Papers:

 * Hamiltonian hydrodynamics of eccentric discs; Ogilvie G. I, Lynch E. M. 2019; 10.1093/mnras/sty3436
   - Introduces Hamiltonian formalism for fluid discs.
 * Linear and non-linear eccentric mode evolution in unstratified MHD discs; Lynch E. M., Dewberry J. W. 2023; 10.1093/mnras/stad2678  
   - Extends theory to include magnetic fields in cylindrical coordinates

# Current Fuctionality

 * 2D and 3D globally isothermal eccentric mode solvers
   - 2D includes magnetic fields
 * 2D adiabatic eccentric mode solvers
 * gamma=1, 3D adiabatic eccentric mode solver
   - This differs from a locally isothermal solver in the entropy is constant on orbits
 * Includes vertical structure solver for adiabatic and isothermal hydro-discs
 * class to compute geometric terms (e.g. metric coefficients/Christoffel symbols) in orbital coordinates
   - Includes routing for coordinate transform
 * class to compute fluid properties in an eccentric disc

# TODO

 * Update examples to include comparison with Ogilvie & Lynch 2019
 * Improve documentation
 * Add magnetic fields to 2D adiabatic solver
 * Create a solver for tabulated Hamiltonians
   - This is required for a practical 3D adiabatic solver
