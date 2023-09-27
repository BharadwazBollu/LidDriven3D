#include "header.h"

void Fields::solvePrediction()
{
    double error = 1;
    double resU, resV, resW;
    double rU, rV, rW;

    while(error > 1e-3)
    {
        resU = 0.0;
        resV = 0.0;
        resW = 0.0;

        for (int k=1; k<=nz_; k++)
        {
            for (int j=1; j<=ny_; j++)
            {
                for (int i=1; i<=nx_; i++)
                {
                    // Calculating fluxes
                    fluxE_ = 0.5 * ( field_.u_curr[i][j][k] + field_.u_curr[i+1][j][k] ) * areaX_;
                    fluxW_ = 0.5 * ( field_.u_curr[i][j][k] + field_.u_curr[i-1][j][k] ) * areaX_;
                    fluxN_ = 0.5 * ( field_.v_curr[i][j][k] + field_.v_curr[i][j+1][k] ) * areaY_;
                    fluxS_ = 0.5 * ( field_.v_curr[i][j][k] + field_.v_curr[i][j-1][k] ) * areaY_;
                    fluxF_ = 0.5 * ( field_.w_curr[i][j][k] + field_.w_curr[i][j][k+1] ) * areaZ_;
                    fluxB_ = 0.5 * ( field_.w_curr[i][j][k] + field_.w_curr[i][j][k-1] ) * areaZ_;

                    // Upwind scheme

                    conv_diag_ = 0.0;
                    if ( fluxE_ >= 0.0 )
                    {
                        uE_ = field_.u_curr[i][j][k];
                        vE_ = field_.u_curr[i][j][k];
                        wE_ = field_.u_curr[i][j][k];
                        conv_diag_ = conv_diag_ + fluxE_;
                    }
                    else
                    {
                        uE_ = field_.u_curr[i+1][j][k];
                        vE_ = field_.u_curr[i+1][j][k];
                        wE_ = field_.u_curr[i+1][j][k];
                    }

                    if ( fluxW_ >= 0.0 )
                    {
                        uW_ = field_.u_curr[i][j][k];
                        vW_ = field_.u_curr[i][j][k];
                        wW_ = field_.u_curr[i][j][k];
                        conv_diag_ = conv_diag_ + fluxW_;
                    }
                    else
                    {
                        uW_ = field_.u_curr[i-1][j][k];
                        vW_ = field_.u_curr[i-1][j][k];
                        wW_ = field_.u_curr[i-1][j][k];
                    }


                    if ( fluxN_ >= 0.0 )
                    {
                        uN_ = field_.u_curr[i][j][k];
                        vN_ = field_.u_curr[i][j][k];
                        wN_ = field_.u_curr[i][j][k];
                        conv_diag_ = conv_diag_ + fluxN_;
                    }
                    else
                    {
                        uN_ = field_.u_curr[i+1][j][k];
                        vN_ = field_.u_curr[i+1][j][k];
                        wN_ = field_.u_curr[i+1][j][k];
                    }

                    if ( fluxS_ >= 0.0 )
                    {
                        uS_ = field_.u_curr[i][j][k];
                        vS_ = field_.u_curr[i][j][k];
                        wS_ = field_.u_curr[i][j][k];
                        conv_diag_ = conv_diag_ + fluxS_;
                    }
                    else
                    {
                        uS_ = field_.u_curr[i-1][j][k];
                        vS_ = field_.u_curr[i-1][j][k];
                        wS_ = field_.u_curr[i-1][j][k];
                    }

                    if ( fluxF_ >= 0.0 )
                    {
                        uF_ = field_.u_curr[i][j][k];
                        vF_ = field_.u_curr[i][j][k];
                        wF_ = field_.u_curr[i][j][k];
                        conv_diag_ = conv_diag_ + fluxF_;
                    }
                    else
                    {
                        uF_ = field_.u_curr[i][j][k+1];
                        vF_ = field_.u_curr[i][j][k+1];
                        wF_ = field_.u_curr[i][j][k+1];
                    }

                    if ( fluxB_ >= 0.0 )
                    {
                        uB_ = field_.u_curr[i][j][k];
                        vB_ = field_.u_curr[i][j][k];
                        wB_ = field_.u_curr[i][j][k];
                        conv_diag_ = conv_diag_ + fluxB_;
                    }
                    else
                    {
                        uB_ = field_.u_curr[i][j][k-1];
                        vB_ = field_.u_curr[i][j][k-1];
                        wB_ = field_.u_curr[i][j][k-1];
                    }

                    convectionX_ = fluxE_ * uE_ + fluxW_ * uW_ + fluxN_ * uN_ 
                    + fluxS_ * uS_ + fluxF_ * uF_ + fluxB_ * uB_;
                    convectionY_ = fluxE_ * vE_ + fluxW_ * vW_ + fluxN_ * vN_ 
                    + fluxS_ * vS_ + fluxF_ * vF_ + fluxB_ * vB_;
                    convectionZ_ = fluxE_ * wE_ + fluxW_ * wW_ + fluxN_ * wN_ 
                    + fluxS_ * wS_ + fluxF_ * wF_ + fluxB_ * wB_;

                    diffusionX_ = ( field_.u_pred[i+1][j][k] -2*field_.u_pred[i][j][k] + field_.u_pred[i-1][j][k] ) * diff_disc_coffX_
                    + ( field_.u_pred[i][j+1][k] -2*field_.u_pred[i][j][k] + field_.u_pred[i][j-1][k] ) * diff_disc_coffY_
                    + ( field_.u_pred[i][j][k+1] -2*field_.u_pred[i][j][k] + field_.u_pred[i][j][k-1] ) * diff_disc_coffZ_;

                    diffusionY_ = ( field_.v_pred[i+1][j][k] -2*field_.v_pred[i][j][k] + field_.v_pred[i-1][j][k] ) * diff_disc_coffX_
                    + ( field_.v_pred[i][j+1][k] -2*field_.v_pred[i][j][k] + field_.v_pred[i][j-1][k] ) * diff_disc_coffY_
                    + ( field_.v_pred[i][j][k+1] -2*field_.v_pred[i][j][k] + field_.v_pred[i][j][k-1] ) * diff_disc_coffZ_;

                    diffusionZ_ = ( field_.w_pred[i+1][j][k] -2*field_.w_pred[i][j][k] + field_.w_pred[i-1][j][k] ) * diff_disc_coffX_
                    + ( field_.w_pred[i][j+1][k] -2*field_.w_pred[i][j][k] + field_.w_pred[i][j-1][k] ) * diff_disc_coffY_
                    + ( field_.w_pred[i][j][k+1] -2*field_.w_pred[i][j][k] + field_.w_pred[i][j][k-1] ) * diff_disc_coffZ_;


                    double tmp;
                    tmp = Vp_ * density_/dt_ + density_ * conv_diag_ + nu_ * diff_cen_coff;

                    rU = Vp_ * density_/dt_ * ( field_.u_curr[i][j][k] - field_.u_pred[i][j][k] ) 
                    - convectionX_ * density_ + nu_ * diffusionX_ ;
                    resU = resU + rU * rU ;
                    field_.u_pred[i][j][k] = rU/tmp + field_.u_pred[i][j][k];

                    rV = Vp_ * density_/dt_ * ( field_.v_curr[i][j][k] - field_.v_pred[i][j][k] ) 
                    - convectionY_ * density_ + nu_ * diffusionY_ ;
                    resV = resV + rV * rV ;
                    field_.v_pred[i][j][k] = rW/tmp + field_.v_pred[i][j][k];

                    rW = Vp_ * density_/dt_ *( field_.w_curr[i][j][k] - field_.w_pred[i][j][k] ) 
                    - convectionZ_ * density_ + nu_ * diffusionZ_ ;
                    resW = resW + rW * rW;
                    field_.w_pred[i][j][k] = rW/tmp + field_.w_pred[i][j][k];

                }
            }
        }

        error = sqrt( (resU + resV + resW)/ (nx_ * ny_ *nz_) );
   }

}