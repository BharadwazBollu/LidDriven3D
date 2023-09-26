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
    while(error < 1e-3)
    {
        double dx = 1.0/nx_;
        double dy = 1.0/nx_;
        double dz = 1.0/nz_;
        double Vp = dx*dy*dz;

        double area_x = dy * dx;
        double area_y = dx * dz;
        double area_z = dx * dy;

        double fe, fw, fn, fs, ff, fb;
        double uff, vff, fduf, fdvf;
        double ue, uw, un, us;
        double ve, vw, vn, vs;

        for (int k=1; k<=nz_; k++)
        {
            for (int j=1; j<=ny_; j++)
            {
                for (int i=1; i<=nx_; i++)
                {
                    fe = 0.5 * ( fields.u_curr[i][j][k] + fields.u_curr[i+1][j][k] ) * area_x;
                    fw = 0.5 * ( fields.u_curr[i][j][k] + fields.u_curr[i-1][j][k] ) * area_x;
                    fn = 0.5 * ( fields.v_curr[i][j][k] + fields.v_curr[i][j+1][k] ) * area_y;
                    fs = 0.5 * ( fields.v_curr[i][j][k] + fields.v_curr[i][j-1][k] ) * area_y;
                    ff = 0.5 * ( fields.w_curr[i][j][k] + fields.w_curr[i][j][k+1] ) * area_z;
                    fb = 0.5 * ( fields.w_curr[i][j][k] + fields.w_curr[i][j][k-1] ) * area_z;
                }
            }
        }




   }



}
