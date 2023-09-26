#include "header.h"

void hello(){
    std::cout << "hello" ;
}

Fields::Fields(int nx, int ny, int nz):nx_(nx),ny_(ny),nz_(nz)
{

    fields.u.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    fields.v.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    fields.w.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    fields.u_pred.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    fields.v_pred.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    fields.w_pred.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    fields.u_curr.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    fields.v_curr.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    fields.w_curr.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    fields.p.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    
}


void Fields::initialise()
{
    // Top BC 
    int k = nz_+1;
    for(int j=0; j<=ny_; j++)
    {
        for(int i=0; i<=nx_; i++){
            fields.u[i][j][k] = 1.0;
            fields.u_curr[i][j][k] = 1.0;
            fields.u_pred[i][j][k] = 1.0;
        }
    }
}

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

        double dx = 1.0/nx_;
        double dy = 1.0/nx_;
        double dz = 1.0/nz_;
        double Vp = dx*dy*dz;
        double dt = 0.001;
        double nu = 0.01;

        double area_x = dy * dx;
        double area_y = dx * dz;
        double area_z = dx * dy;

        double fe, fw, fn, fs, ff, fb;
        double ue, uw, un, us, uf, ub;
        double ve, vw, vn, vs, vf, vb;
        double we, ww, wn, ws, wf, wb;

        double cfx, cfy, cfz;   // convective in x, y, z
        double dfx, dfy, dfz;   // diffusion in x, y, z

        double dfcx, dfcy, dfcz;    // diffusion denomiators coefficeint
        dfcx = 1 / ( dx * dx );
        dfcy = 1 / ( dy * dy );
        dfcz = 1 / ( dz * dz );

        double cenCoff;

        double temp = Vp * nu * 2 * ( dfcx + dfcy + dfcz );

        for (int k=1; k<=nz_; k++)
        {
            for (int j=1; j<=ny_; j++)
            {
                for (int i=1; i<=nx_; i++)
                {
                    // Calculating fluxes
                    fe = 0.5 * ( fields.u_curr[i][j][k] + fields.u_curr[i+1][j][k] ) * area_x;
                    fw = 0.5 * ( fields.u_curr[i][j][k] + fields.u_curr[i-1][j][k] ) * area_x;
                    fn = 0.5 * ( fields.v_curr[i][j][k] + fields.v_curr[i][j+1][k] ) * area_y;
                    fs = 0.5 * ( fields.v_curr[i][j][k] + fields.v_curr[i][j-1][k] ) * area_y;
                    ff = 0.5 * ( fields.w_curr[i][j][k] + fields.w_curr[i][j][k+1] ) * area_z;
                    fb = 0.5 * ( fields.w_curr[i][j][k] + fields.w_curr[i][j][k-1] ) * area_z;

                    // Upwind scheme

                    cenCoff = 0.0;
                    if ( fe >= 0.0 )
                    {
                        ue = fields.u_curr[i][j][k];
                        ve = fields.u_curr[i][j][k];
                        we = fields.u_curr[i][j][k];
                        cenCoff = cenCoff + fe;
                    }
                    else
                    {
                        ue = fields.u_curr[i+1][j][k];
                        ve = fields.u_curr[i+1][j][k];
                        we = fields.u_curr[i+1][j][k];
                    }

                    if ( fw >= 0.0 )
                    {
                        uw = fields.u_curr[i][j][k];
                        vw = fields.u_curr[i][j][k];
                        ww = fields.u_curr[i][j][k];
                        cenCoff = cenCoff + fw;
                    }
                    else
                    {
                        uw = fields.u_curr[i-1][j][k];
                        vw = fields.u_curr[i-1][j][k];
                        ww = fields.u_curr[i-1][j][k];
                    }


                    if ( fn >= 0.0 )
                    {
                        un = fields.u_curr[i][j][k];
                        vn = fields.u_curr[i][j][k];
                        wn = fields.u_curr[i][j][k];
                        cenCoff = cenCoff + fn;
                    }
                    else
                    {
                        un = fields.u_curr[i+1][j][k];
                        vn = fields.u_curr[i+1][j][k];
                        wn = fields.u_curr[i+1][j][k];
                    }

                    if ( fs >= 0.0 )
                    {
                        us = fields.u_curr[i][j][k];
                        vs = fields.u_curr[i][j][k];
                        ws = fields.u_curr[i][j][k];
                        cenCoff = cenCoff + fs;
                    }
                    else
                    {
                        us = fields.u_curr[i-1][j][k];
                        vs = fields.u_curr[i-1][j][k];
                        ws = fields.u_curr[i-1][j][k];
                    }

                    if ( ff >= 0.0 )
                    {
                        uf = fields.u_curr[i][j][k];
                        vf = fields.u_curr[i][j][k];
                        wf = fields.u_curr[i][j][k];
                        cenCoff = cenCoff + ff;
                    }
                    else
                    {
                        uf = fields.u_curr[i][j][k+1];
                        vf = fields.u_curr[i][j][k+1];
                        wf = fields.u_curr[i][j][k+1];
                    }

                    if ( fb >= 0.0 )
                    {
                        ub = fields.u_curr[i][j][k];
                        vb = fields.u_curr[i][j][k];
                        wb = fields.u_curr[i][j][k];
                        cenCoff = cenCoff + fb;
                    }
                    else
                    {
                        ub = fields.u_curr[i][j][k-1];
                        vb = fields.u_curr[i][j][k-1];
                        wb = fields.u_curr[i][j][k-1];
                    }

                    cfx = fe * ue + fw * uw + fn * un + fs * us + ff * uf + fb * ub;
                    cfy = fe * ve + fw * vw + fn * vn + fs * vs + ff * vf + fb * vb;
                    cfz = fe * we + fw * ww + fn * wn + fs * ws + ff * wf + fb * wb;

                    dfx = ( fields.u_pred[i+1][j][k] -2*fields.u_pred[i][j][k] + fields.u_pred[i-1][j][k] ) * dfcx + ( fields.u_pred[i][j+1][k] -2*fields.u_pred[i][j][k] 
                    + fields.u_pred[i][j-1][k] ) * dfcy + ( fields.u_pred[i][j][k+1] -2*fields.u_pred[i][j][k] + fields.u_pred[i][j][k-1] ) * dfcz ;
                    dfy = ( fields.v_pred[i+1][j][k] -2*fields.v_pred[i][j][k] + fields.v_pred[i-1][j][k] ) * dfcx + ( fields.v_pred[i][j+1][k] -2*fields.v_pred[i][j][k] 
                    + fields.v_pred[i][j-1][k] ) * dfcy + ( fields.v_pred[i][j][k+1] -2*fields.v_pred[i][j][k] + fields.v_pred[i][j][k-1] ) * dfcz ;
                    dfz = ( fields.w_pred[i+1][j][k] -2*fields.w_pred[i][j][k] + fields.w_pred[i-1][j][k] ) * dfcx + ( fields.w_pred[i][j+1][k] -2*fields.w_pred[i][j][k] 
                    + fields.w_pred[i][j-1][k] ) * dfcy + ( fields.w_pred[i][j][k+1] -2*fields.w_pred[i][j][k] + fields.w_pred[i][j][k-1] ) * dfcz ;


                    rU = Vp/dt * ( fields.u_curr[i][j][k] - fields.u_pred[i][j][k] ) - cfx + Vp*nu * dfx ;
                    resU = resU + rU * rU ;
                    fields.u_pred[i][j][k] = rU/( Vp/dt + cenCoff + temp ) + fields.u_pred[i][j][k];

                    rV = Vp/dt * ( fields.v_curr[i][j][k] - fields.v_pred[i][j][k] ) - cfy + Vp*nu * dfz ;
                    resV = resV + rV * rV ;
                    fields.v_pred[i][j][k] = rW/( Vp/dt + cenCoff + temp ) + fields.v_pred[i][j][k];

                    rW = Vp/dt * ( fields.w_curr[i][j][k] - fields.w_pred[i][j][k] ) - cfz + Vp*nu * dfz ;
                    resW = resW + rW * rW;
                    fields.w_pred[i][j][k] = rW/( Vp/dt + cenCoff + temp ) + fields.w_pred[i][j][k];




                }
            }
        }

        error = sqrt( (resU + resV + resW)/ (nx_ * ny_ *nz_) )

   }



}
