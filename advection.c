/*******************************************************************************
2D advection example program which advects a Gaussian u(x,y) at a fixed velocity



Outputs: initial.dat - initial values of u(x,y)
         final.dat   - final values of u(x,y)

         The output files have three columns: x, y, u

         Compile with: gcc -o advection2D -std=c99 advection2D.c -lm

Notes: The time step is calculated using the CFL condition

********************************************************************************/

/*********************************************************************
                     Include header files
**********************************************************************/

#include <math.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/*********************************************************************
                      Main function
**********************************************************************/

int main()
{
    /* Grid properties */
    const int NX = 1000;     // Number of x points
    const int NY = 1000;     // Number of y points
    const float xmin = 0.0;  // Minimum x value
    const float xmax = 30.0; // Maximum x value
    const float ymin = 0.0;  // Minimum y value
    const float ymax = 30.0; // Maximum y value

    /* Parameters for the Gaussian initial conditions */
    const float x0 = 0.1;                  // Centre(x)
    const float y0 = 0.1;                  // Centre(y)
    const float sigmax = 0.0;              // Width(x)
    const float sigmay = 0.0;              // Width(y)
    const float sigmax2 = sigmax * sigmax; // Width(x) squared
    const float sigmay2 = sigmay * sigmay; // Width(y) squared

    /* Boundary conditions */
    const float bval_left = 0.0;  // Left boudnary value
    const float bval_right = 0.0; // Right boundary value
    const float bval_lower = 0.0; // Lower boundary
    const float bval_upper = 0.0; // Upper bounary

    /* Constants for updating Left boundary value*/
    const float y_0 = 15.0;
    const float t_0 = 3.0;
    const float sigma_y = 5.0;
    const float sigma_t = 1.0;
    const float sigma_y_sqr = sigma_y * sigma_y;
    const float sigma_t_sqr = sigma_t * sigma_t;
    float y_2; // y squared (used for boundary calculation)
    float t_2; // t squared (used for boundary calculation)

    /* Time stepping parameters */
    const float CFL = 0.9;   // CFL number
    const int nsteps = 1000; // Number of time steps

    /* Velocity */
    // const float velx = 1.0; // Velocity in x direction
    float velx = 1.0;       // Velocity in x direction (varies based on height)
    const float vely = 0.0; // Velocity in y direction

    /* Parameters for logarithmic profile used for updating (velx) */
    const float u_star = 0.1; // u* constant in (m/s)
    const float z_0 = 1.0;    // z0 constant in (m)
    const float K = 0.41;     // Von Karman's constant

    /* Arrays to store variables. These have NX+2 elements
        to allow boundary values to be stored at both ends */
    float x[NX + 2];            // x-axis values
    float y[NX + 2];            // y-axis values
    float u[NX + 2][NY + 2];    // Array of u values
    float dudt[NX + 2][NY + 2]; // Rate of change of u

    /* Array to store vertically averaged distribution of u(x,y) */
    float vert_avg[NX]; // value of average u around x-axis

    float x2; // x squared (used to calculate initial conditions)
    float y2; // y squared (used to calculate initial conditions)

    /* Calculate distance between points */
    float dx = (xmax - xmin) / ((float)NX);
    float dy = (ymax - ymin) / ((float)NY);

    /* Calculate the time step using the CFL condition */
    /* The fabs function gives the absolute value in case the velocity is -ve */
    float dt = CFL / ((fabs(velx) / dx) + (fabs(vely) / dy));

    /*** Report information about the calculation ***/
    printf("Grid spacing dx     = %g\n", dx);
    printf("Grid spacing dy     = %g\n", dy);
    printf("CFL number          = %g\n", CFL);
    printf("Time step           = %g\n", dt);
    printf("No. of time steps   = %d\n", nsteps);
    printf("End time            = %g\n", dt * (float)nsteps);
    printf("Distance advected x = %g\n", velx * dt * (float)nsteps);
    printf("Distance advected y = %g\n", vely * dt * (float)nsteps);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(x, dx, NX)
#endif
    /*** Place x points in the middle of the cell ***/
    /* LOOP 1 */
    for (int i = 0; i < NX + 2; i++)
    {
        x[i] = ((float)i - 0.5) * dx;
    }

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(y, dy, NY)
#endif
    /*** Place y points in the middle of the cell ***/
    /* LOOP 2 */
    for (int j = 0; j < NY + 2; j++)
    {
        y[j] = ((float)j - 0.5) * dy;
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) collapse(2) private(x2, y2)
#endif
    /*** Set up Gaussian initial conditions ***/
    /* LOOP 3 */
    for (int i = 0; i < NX + 2; i++)
    {
        for (int j = 0; j < NY + 2; j++)
        {
            x2 = (x[i] - x0) * (x[i] - x0);
            y2 = (y[j] - y0) * (y[j] - y0);
            u[i][j] = exp(-1.0 * ((x2 / (2.0 * sigmax2)) + (y2 / (2.0 * sigmay2))));
        }
    }

    /*** Write array of initial u values out to file ***/
    FILE *initialfile;
    initialfile = fopen("initial.dat", "w");

    /* LOOP 4 */
    /***************************************************************************
    This loop can not be parallelized since it is an I/O operation. In such operations, each thread is trying to
    access the same resource and write into it. With respect to running time and efficiency, it makes more sense to do
    this part in a serial manner.
    ***************************************************************************/
    for (int i = 0; i < NX + 2; i++)
    {
        for (int j = 0; j < NY + 2; j++)
        {
            fprintf(initialfile, "%g %g %g\n", x[i], y[j], u[i][j]);
        }
    }
    fclose(initialfile);

    /*** Update solution by looping over time steps ***/
    /* LOOP 5 */
    /***************************************************************************
    This loop is the time loop and is doing an iterative operation. The value of (u) is being updated during each
    iteration of this loop so the order of the operation is important. Therefore there is no need to parallelize
    this loop.
    ***************************************************************************/
    for (int m = 0; m < nsteps; m++)
    {

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(y_2, t_2)
#endif
        /*** Apply boundary conditions at u[0][:] and u[NX+1][:] ***/
        /* LOOP 6 */
        for (int j = 0; j < NY + 2; j++)
        {
            y_2 = (y[j] - y_0) * (y[j] - y_0);
            t_2 = ((m * dt) - t_0) * ((m * dt) - t_0);

            u[0][j] = exp(-1.0 * ((y_2 / (2 * sigma_y_sqr)) + (t_2 / (2 * sigma_t_sqr))));
            u[NX + 1][j] = bval_right;
        }

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
        /*** Apply boundary conditions at u[:][0] and u[:][NY+1] ***/
        /* LOOP 7 */
        for (int i = 0; i < NX + 2; i++)
        {
            u[i][0] = bval_lower;
            u[i][NY + 1] = bval_upper;
        }

#ifdef _OPENMP
#pragma omp parallel for collapse(2) default(shared) private(velx)
#endif
        /*** Calculate the rate of change of u using leftward difference ***/
        /* Loop over points in the domain but not boundary values */
        /* LOOP 8 */
        for (int i = 1; i < NX + 1; i++)
        {
            for (int j = 1; j < NY + 1; j++)
            {
                /** Updating the velocity according to log profile equation **/
                if (y[j] > z_0)
                    velx = (u_star / K) * log(y[j] / z_0);
                else
                    velx = 0.0;

                dudt[i][j] = -velx * (u[i][j] - u[i - 1][j]) / dx - vely * (u[i][j] - u[i][j - 1]) / dy;
            }
        }
#ifdef _OPENMP
#pragma omp parallel for collapse(2) default(none) shared(NX, NY, dudt, dt, u)
#endif
        /*** Update u from t to t+dt ***/
        /* Loop over points in the domain but not boundary values */
        /* LOOP 9 */
        for (int i = 1; i < NX + 1; i++)
        {
            for (int j = 1; j < NY + 1; j++)
            {
                u[i][j] = u[i][j] + dudt[i][j] * dt;
            }
        }
    } // time loop

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    /* Calculating the vertically averaged distribution of u(x,y) */
    for (int i = 1; i < NX + 1; i++)
    {
        float sum = 0.0;

#ifdef _OPENMP
#pragma omp parallel for default(shared) reduction(+ : sum)
#endif
        for (int j = 1; j < NY + 1; j++)
        {
            sum += u[i][j];
        }

        // printf("%d  sum: %g\n", i, sum);
        vert_avg[i] = (float)(sum / NY);
        sum = 0.0;
    }

    /*** Write array of final u values out to file ***/
    FILE *avgfile;
    avgfile = fopen("avg.dat", "w");
    for (int i = 0; i < NX + 2; i++)
    {
        fprintf(avgfile, "%g %g\n", x[i], vert_avg[i]);
    }
    fclose(avgfile);

    /*** Write array of final u values out to file ***/
    FILE *finalfile;
    finalfile = fopen("final.dat", "w");

    /* LOOP 10 */

    /***************************************************************************
    This loop is doing the same operation as "LOOP 4" therefore there is no need for parallelising it.
    ***************************************************************************/
    for (int i = 0; i < NX + 2; i++)
    {
        for (int j = 0; j < NY + 2; j++)
        {
            fprintf(finalfile, "%g %g %g\n", x[i], y[j], u[i][j]);
        }
    }
    fclose(finalfile);

    return 0;
}

/* End of file ******************************************************/
