#include "header.h"

void Fields::solvePressure()
{
    double error = 1;
    double resP, rP;
    double flux_no_pressure_correction;

    while(error > 1e-3)
    {
        resP = 0.0;

        for (int k=1; k<=nz_; k++)
        {
            for (int j=1; j<=ny_; j++)
            {
                for (int i=1; i<=nx_; i++)
                {
                    flux_no_pressure_correction = 
                    resP  = 0.5 * ( (field_.u_pred[i+1][j][k] - field_.u_pred[i-1][j] )/dx_ 
                    + ( field_.v_pred[i][j+1][k] - field_.v_pred[i][j-1][k] )/dy_ ) 
                    + ( field_.w_pred[i][j][k+1] - field_.w_pred[i][j][k-1] )/dz_ ) - dt_ * ( p_star[i+1][j] - 2*p_star[i][j] + p_star[i-1][j])/(dx*dx) - dt * ( p_star[i][j+1] - 2*p_star[i][j] + p_star[i][j-1])/(dy*dy);
                res = res + R * R;                                                     // Adding square of residues for RMS
                p_star[i][j] = - 1.5*R/( 2*dt/(dx*dx) + 2*dt/(dy*dy) ) + p_star[i][j];           // correction for next iteration
                }
            }
        }

    }
}