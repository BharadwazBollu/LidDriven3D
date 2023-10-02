#include "header.h"

void LidDrivenCavity::initialize()
{
    // Top boundary condition (moving lid)
    int k = nz_+1;
    for(int j=0; j<=ny_+1; j++)
    {
        for(int i=0; i<=nx_+1; i++){
            field_.u[i][j][k] = lid_velocity_;
            field_.u_curr[i][j][k] = lid_velocity_;
            field_.u_pred[i][j][k] = lid_velocity_;
        }
    }
}