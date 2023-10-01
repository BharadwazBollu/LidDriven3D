#include "header.h"

void Fields::solvePressure()
{
    double error = 1;
    double resP, rP;
    double flux_no_pressure_correction, diffsuionP;

    std::cout << "Entered Pressure" << std::endl;

    while(error > 1e-3)
    {
        resP = 0.0;

        for (int k=1; k<=nz_; k++)
        {
            for (int j=1; j<=ny_; j++)
            {
                for (int i=1; i<=nx_; i++)
                {
                    flux_no_pressure_correction = 0.5 * ( field_.u_pred[i+1][j][k] - field_.u_pred[i-1][j][k] ) * areaX_
                    + 0.5 * ( field_.v_pred[i][j+1][k] - field_.v_pred[i][j-1][k] ) * areaY_
                    + 0.5 * ( field_.w_pred[i][j][k+1] - field_.w_pred[i][j][k-1] ) * areaZ_;

                    diffsuionP = ( field_.p[i+1][j][k] -2*field_.p[i][j][k] + field_.p[i-1][j][k] ) * diff_disc_coffX_
                    + ( field_.p[i][j+1][k] -2*field_.p[i][j][k] + field_.p[i][j-1][k] ) * diff_disc_coffY_
                    + ( field_.p[i][j][k+1] -2*field_.p[i][j][k] + field_.p[i][j][k-1] ) * diff_disc_coffZ_;

                    rP  = density_/dt_ * flux_no_pressure_correction - diffsuionP;
                    resP = resP + rP * rP;                                                   
                    field_.p[i][j][k] = -rP/diff_cen_coff + field_.p[i][j][k]; 
                }
            }
        }

        error = sqrt( resP/ (nx_ * ny_ *nz_) );
        // std::cout << " Pressure residual " << error << std::endl;

        // BC for pressure

        // East & West
        for (int k=1; k<=nz_; k++)
        {
            for (int j=1 ; j<=ny_; j++)     
            {
                field_.p[0][j][k] = field_.p[1][j][k];
                field_.p[nx_+1][j][k] = field_.p[nx_][j][k];
            }
        }

        // North & South
        for (int k=1; k<=nz_; k++)
        {
            for (int i=1; i<=nx_; i++)    
            {
                field_.p[i][0][k] = field_.p[i][1][k];
                field_.p[i][ny_+1][k] = field_.p[i][ny_][k];
            }
        }

        // Front & Back
        for (int j=1; j<=ny_; j++)
        {
            for (int i=1; i<=nx_; i++) 
            {
                field_.p[i][j][0] = field_.p[i][j][1];
                field_.p[i][j][nz_+1] = field_.p[i][j][nz_];
            }
        }

    }

    //std::cout << " Done Pressure " << std::endl;

}