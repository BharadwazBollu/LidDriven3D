#include "header.h"

void hello(){
    std::cout << "hello" ;
}

Fields::Fields(int nx, int ny, int nz, double nu, double density, double dt)
:nx_(nx), ny_(ny), nz_(nz), nu_(nu), density_(density), dt_(dt)
{

    field_.u.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    field_.v.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    field_.w.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    field_.u_pred.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    field_.v_pred.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    field_.w_pred.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    field_.u_curr.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    field_.v_curr.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    field_.w_curr.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    field_.p.resize(nz_+2, std::vector<std::vector<double>>(ny_+2, std::vector<double>(nx_+2, 0.0)));
    
}


void Fields::initialise()
{
    // Top BC 
    int k = nz_+1;
    for(int j=0; j<=ny_; j++)
    {
        for(int i=0; i<=nx_; i++){
            field_.u[i][j][k] = 1.0;
            field_.u_curr[i][j][k] = 1.0;
            field_.u_pred[i][j][k] = 1.0;
        }
    }
}


