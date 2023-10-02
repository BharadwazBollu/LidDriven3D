#include "header.h"

void LidDrivenCavity::solvePrediction()
{
    error_ = 1;
    u_local_iteration_ = 0;
    v_local_iteration_ = 0;
    w_local_iteration_ = 0;

    while(error_ > prediction_tolerance_)
    {
        u_rms_ = 0.0;
        v_rms_ = 0.0;
        w_rms_ = 0.0;
        // loop for all interior cells
        for (int k=1; k<=nz_; k++)
        {
            for (int j=1; j<=ny_; j++)
            {
                for (int i=1; i<=nx_; i++)
                {
                    // Calculating fluxes
                    fluxE_ = 0.5 * ( field_.u_curr[i][j][k] + field_.u_curr[i+1][j][k] ) * x_area_;
                    fluxW_ = -0.5 * ( field_.u_curr[i][j][k] + field_.u_curr[i-1][j][k] ) * x_area_;
                    fluxN_ = 0.5 * ( field_.v_curr[i][j][k] + field_.v_curr[i][j+1][k] ) * y_area_;
                    fluxS_ = -0.5 * ( field_.v_curr[i][j][k] + field_.v_curr[i][j-1][k] ) * y_area_;
                    fluxF_ = 0.5 * ( field_.w_curr[i][j][k] + field_.w_curr[i][j][k+1] ) * z_area_;
                    fluxB_ = -0.5 * ( field_.w_curr[i][j][k] + field_.w_curr[i][j][k-1] ) * z_area_;
                    
                    // Upwind scheme for convection
                    // East
                    conv_diag_ = 0.0;
                    if ( fluxE_ >= 0.0 )
                    {
                        uE_ = field_.u_curr[i][j][k];
                        vE_ = field_.v_curr[i][j][k];
                        wE_ = field_.w_curr[i][j][k];
                        conv_diag_ = conv_diag_ + fluxE_;
                    }
                    else
                    {
                        uE_ = field_.u_curr[i+1][j][k];
                        vE_ = field_.v_curr[i+1][j][k];
                        wE_ = field_.w_curr[i+1][j][k];
                    }
                    // West
                    if ( fluxW_ >= 0.0 )
                    {
                        uW_ = field_.u_curr[i][j][k];
                        vW_ = field_.v_curr[i][j][k];
                        wW_ = field_.w_curr[i][j][k];
                        conv_diag_ = conv_diag_ + fluxW_;
                    }
                    else
                    {
                        uW_ = field_.u_curr[i-1][j][k];
                        vW_ = field_.v_curr[i-1][j][k];
                        wW_ = field_.w_curr[i-1][j][k];
                    }
                    // North
                    if ( fluxN_ >= 0.0 )
                    {
                        uN_ = field_.u_curr[i][j][k];
                        vN_ = field_.v_curr[i][j][k];
                        wN_ = field_.w_curr[i][j][k];
                        conv_diag_ = conv_diag_ + fluxN_;
                    }
                    else
                    {
                        uN_ = field_.u_curr[i][j+1][k];
                        vN_ = field_.v_curr[i][j+1][k];
                        wN_ = field_.w_curr[i][j+1][k];
                    }
                    // South
                    if ( fluxS_ >= 0.0 )
                    {
                        uS_ = field_.u_curr[i][j][k];
                        vS_ = field_.v_curr[i][j][k];
                        wS_ = field_.w_curr[i][j][k];
                        conv_diag_ = conv_diag_ + fluxS_;
                    }
                    else
                    {
                        uS_ = field_.u_curr[i][j-1][k];
                        vS_ = field_.v_curr[i][j-1][k];
                        wS_ = field_.w_curr[i][j-1][k];
                    }
                    // Front
                    if ( fluxF_ >= 0.0 )
                    {
                        uF_ = field_.u_curr[i][j][k];
                        vF_ = field_.v_curr[i][j][k];
                        wF_ = field_.w_curr[i][j][k];
                        conv_diag_ = conv_diag_ + fluxF_;
                    }
                    else
                    {
                        uF_ = field_.u_curr[i][j][k+1];
                        vF_ = field_.v_curr[i][j][k+1];
                        wF_ = field_.w_curr[i][j][k+1];
                    }
                    // Back
                    if ( fluxB_ >= 0.0 )
                    {
                        uB_ = field_.u_curr[i][j][k];
                        vB_ = field_.v_curr[i][j][k];
                        wB_ = field_.w_curr[i][j][k];
                        conv_diag_ = conv_diag_ + fluxB_;
                    }
                    else
                    {
                        uB_ = field_.u_curr[i][j][k-1];
                        vB_ = field_.v_curr[i][j][k-1];
                        wB_ = field_.w_curr[i][j][k-1];
                    }
                    // total convection in x, y, z
                    x_convection_ = fluxE_ * uE_ + fluxW_ * uW_ + fluxN_ * uN_ 
                    + fluxS_ * uS_ + fluxF_ * uF_ + fluxB_ * uB_;
                    y_convection_ = fluxE_ * vE_ + fluxW_ * vW_ + fluxN_ * vN_ 
                    + fluxS_ * vS_ + fluxF_ * vF_ + fluxB_ * vB_;
                    z_convection_ = fluxE_ * wE_ + fluxW_ * wW_ + fluxN_ * wN_ 
                    + fluxS_ * wS_ + fluxF_ * wF_ + fluxB_ * wB_;
                    // total diffusion in x, y, z
                    x_diffusion_ = ( field_.u_pred[i+1][j][k] -2*field_.u_pred[i][j][k] + field_.u_pred[i-1][j][k] ) * x_diffusion_coefficient_
                    + ( field_.u_pred[i][j+1][k] -2*field_.u_pred[i][j][k] + field_.u_pred[i][j-1][k] ) * y_diffusion_coefficient_
                    + ( field_.u_pred[i][j][k+1] -2*field_.u_pred[i][j][k] + field_.u_pred[i][j][k-1] ) * z_diffusion_coefficient_;
                    y_diffusion_ = ( field_.v_pred[i+1][j][k] -2*field_.v_pred[i][j][k] + field_.v_pred[i-1][j][k] ) * x_diffusion_coefficient_
                    + ( field_.v_pred[i][j+1][k] -2*field_.v_pred[i][j][k] + field_.v_pred[i][j-1][k] ) * y_diffusion_coefficient_
                    + ( field_.v_pred[i][j][k+1] -2*field_.v_pred[i][j][k] + field_.v_pred[i][j][k-1] ) * z_diffusion_coefficient_;
                    z_diffusion_ = ( field_.w_pred[i+1][j][k] -2*field_.w_pred[i][j][k] + field_.w_pred[i-1][j][k] ) * x_diffusion_coefficient_
                    + ( field_.w_pred[i][j+1][k] -2*field_.w_pred[i][j][k] + field_.w_pred[i][j-1][k] ) * y_diffusion_coefficient_
                    + ( field_.w_pred[i][j][k+1] -2*field_.w_pred[i][j][k] + field_.w_pred[i][j][k-1] ) * z_diffusion_coefficient_;
                    // central coeffient diagonal element
                    central_coefficient_ = Vp_ * density_/dt_ + conv_diag_ * density_ + mu_ * diff_cen_coff;
                    // residual for x velocity
                    u_residual_ = Vp_ * density_ * ( field_.u_curr[i][j][k] - field_.u_pred[i][j][k] )/dt_
                    - x_convection_ * density_ + mu_ * x_diffusion_ ;
                    u_rms_ = u_rms_ + u_residual_ * u_residual_ ;
                    field_.u_pred[i][j][k] = u_residual_/central_coefficient_ + field_.u_pred[i][j][k];
                    // residual for y velocity
                    v_residual_ = Vp_ * density_ * ( field_.v_curr[i][j][k] - field_.v_pred[i][j][k] )/dt_
                    - y_convection_ * density_ + mu_ * y_diffusion_ ;
                    v_rms_ = v_rms_ + v_residual_ * v_residual_ ;
                    field_.v_pred[i][j][k] = v_residual_/central_coefficient_ + field_.v_pred[i][j][k];
                    // residual for z velocity
                    w_residual_ = Vp_ * density_ * ( field_.w_curr[i][j][k] - field_.w_pred[i][j][k] )/dt_
                    - z_convection_ * density_ + mu_ * y_diffusion_;
                    w_rms_ = w_rms_ + w_residual_ * w_residual_;
                    field_.w_pred[i][j][k] = w_residual_/central_coefficient_ + field_.w_pred[i][j][k];
                }
            }
        }
        u_rms_ = sqrt( u_rms_/cell_count_ );
        v_rms_ = sqrt( v_rms_/cell_count_ );
        w_rms_ = sqrt( w_rms_/cell_count_ );
        // counting the number of iterations 
        if (u_rms_ > prediction_tolerance_)
        {
            u_local_iteration_++;
        }
        if (v_rms_ > prediction_tolerance_)
        {
            v_local_iteration_++;
        }
        if (w_rms_ > prediction_tolerance_)
        {
            w_local_iteration_++;
        }
        // error for the gauss seidal
        error_ = std::max(u_rms_, std::max(v_rms_, w_rms_));
   }
   
    std::cout << "Gauss seidal solver:  Solving for Ux prediction, residual = " 
    << u_rms_ << " No iterations " << u_local_iteration_ << std::endl;
    std::cout << "Gauss seidal solver:  Solving for Uy prediction, residual = " 
    << v_rms_ << " No iterations " << v_local_iteration_ << std::endl;
    std::cout << "Gauss seidal solver:  Solving for Uz prediction, residual = " 
    << w_rms_ << " No iterations " << w_local_iteration_ << std::endl;
}