#include "header.h"

void LidDrivenCavity::solvePressure()
{
    error_ = 1;
    p_local_iteartion_ = 0;

    while(error_ > pressure_tolerance_)
    {
        p_rms_ = 0.0;
        // loop for ineterior cells
        for (int k=1; k<=nz_; k++)
        {
            for (int j=1; j<=ny_; j++)
            {
                for (int i=1; i<=nx_; i++)
                {
                    flux_no_pressure_correction_ = 0.5 * ( field_.u_pred[i+1][j][k] - field_.u_pred[i-1][j][k] ) * x_area_
                    + 0.5 * ( field_.v_pred[i][j+1][k] - field_.v_pred[i][j-1][k] ) * y_area_
                    + 0.5 * ( field_.w_pred[i][j][k+1] - field_.w_pred[i][j][k-1] ) * z_area_;

                    p_diffusion_ = ( field_.p[i+1][j][k] -2*field_.p[i][j][k] + field_.p[i-1][j][k] ) * x_diffusion_coefficient_
                    + ( field_.p[i][j+1][k] -2*field_.p[i][j][k] + field_.p[i][j-1][k] ) * y_diffusion_coefficient_
                    + ( field_.p[i][j][k+1] -2*field_.p[i][j][k] + field_.p[i][j][k-1] ) * z_diffusion_coefficient_;

                    p_residual_  = density_/dt_ * flux_no_pressure_correction_ - p_diffusion_;
                    p_rms_ = p_rms_ + p_residual_ * p_residual_;                                                   
                    field_.p[i][j][k] = -p_residual_/diff_cen_coff + field_.p[i][j][k]; 
                }
            }
        }

        // Boundary conditions for pressure
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
        // error for gauss seidal
        p_rms_ = sqrt( p_rms_/cell_count_ );
        error_ = p_rms_;
        // counting the number of iterations
        p_local_iteartion_++;
    }

    std::cout << "Gauss seidal solver:  Solving for pressure, residual = " 
    << p_rms_ << " No iterations " << p_local_iteartion_ << std::endl;

}