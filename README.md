# 2D Advection Simulation Program

This project simulates 2D advection using a Gaussian distribution within a domain that mimics the atmospheric boundary layer. The simulation calculates how a cloud of material, represented by the Gaussian, moves horizontally over time based on a velocity profile that varies logarithmically with height.

## Features

- **Customizable Gaussian Initial Conditions**: Initializes a Gaussian distribution representing the material cloud.
- **Boundary Conditions**: Implements boundary conditions for all sides of the computational domain.
- **Logarithmic Velocity Profile**: Applies a height-dependent velocity profile to simulate the horizontal wind speed in the atmospheric boundary layer.

$$
v_x(z) = \frac{u_*}{\kappa} ln \left(\frac{z}{z0}\right)
$$

- **Time Evolution**: Uses finite difference approximations and the CFL condition to compute the advection over a set number of time steps.

$$
\frac{\partial{u}}{\partial{t}}= - \left(v_x \frac{\partial{u}}{\partial{x}} + v_y \frac{\partial{u}}{\partial{y}}\right)
$$

## Program Output

The program generates two output files:

1. `initial.dat`: Contains the initial values of the scalar field `u(x, y)` at the start of the simulation.
2. `final.dat`: Contains the final values of `u(x, y)` after the simulation completes.

Both files contain three columns: `x`, `y`, and `u`, representing the spatial coordinates and the scalar field's value at those coordinates.

## Compilation and Execution

### Prerequisites

- **GCC Compiler**: Ensure that GCC is installed on your system. The program is written in C and tested with `gcc`.

### Compilation

To compile the program, run the following command:

```bash
gcc -o advection2D -std=c99 advection2D.c -lm
```

### Running the Program

After compiling, you can run the executable as follows:

```bash
./advection2D
```

This will generate two output files: `initial.dat` and `final.dat`, containing the initial and final values of the scalar field.

## Implementation Details

- **Grid Properties**: The program defines a grid of ($1000 	\times 1000$) points over a domain of ($30 m \times 30 m$).
- **Boundary Conditions**: Material is introduced into the domain from the left boundary, while other boundaries are assigned a value of zero.
- **Logarithmic Wind Profile**: The horizontal wind speed is computed using a logarithmic profile, with constants:

  - $u^\ast = 0.1 m/s$
  - Roughness length $z_0 = 1.0 m$
  - Von Kármán constant $\kappa = 0.41$
- **Time Stepping**: The program uses the CFL condition for stable time stepping, with a CFL number of 0.9 and a total of 1000 time steps.

## Notes

- The program initializes an empty computational domain and applies a boundary condition that introduces material into the domain over time.
- Ensure sufficient disk space to store the output files, especially for larger grid sizes.
- The program serves as a basic framework for further modifications, such as incorporating different velocity profiles, boundary conditions, or higher grid resolutions.
