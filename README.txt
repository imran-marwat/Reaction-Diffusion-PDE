Numerical Solution of the Reaction-Diffusion equation - Imran Marwat

The Reaction-Diffusion equation is a 2nd order differential equation that models the change in space and time of the concentration of one or more chemical substances.
Based on the finite difference method, a solution is found for this equation. An added consideration is made for when the equation contains an advection term (see https://en.wikipedia.org/wiki/Reaction%E2%80%93diffusion_system#Applications_and_universality).

Adv_Class.py contains all the methods for simulating this equation where phi, the concentration of some chemical, is represented as a 2D ‘order_array’. Periodic boundary conditions are imposed on the lattice.

-Using the finite difference method, an algorithm is created to find the value of the order array, at a given point in the lattice, for the (n+1)th timestep based on the lattice values in the nth timestep. The solution in the time domain ultimately converges to a stable solution.

Adv_Script.py contains instructions for running the methods above for various tasks.

-Available tasks are:

<task> = "viz": animates the solution of the equation. User needs to close animation window to end program.

<task> = “data”: simulates the system until a user defined tolerance limit is reached i.e. solution converges. Plots the steady state array as well as a graph of the concentration vs distance from the centre of the lattice.

To Run: python3 Adv_Script.py <dimension> <kappa> <sigma> <tolerance> <advection velocity> <task>
Recommended: dimension = 50, kappa = 0.01, sigma = 10, tolerance = 0.001, vo = 0.01, 0.1 or 0.5



