#include "header.h"

double Fields::calculateError()
{
    double error = 0.0;

    for (int k=1; k<=nz_; k++)
    {
        for (int j=1; j<=ny_; j++)
        {
            for (int i=1; i<=nx_; i++)
            {
                error =  error + pow((field_.u[i][j][k] - field_.u_curr[i][j][k]),2)
                + pow((field_.v[i][j][k] - field_.v_curr[i][j][k]),2)
                + pow((field_.w[i][j][k] - field_.w_curr[i][j][k]),2);
                field_.u_curr[i][j][k] = field_.u[i][j][k];
                field_.v_curr[i][j][k] = field_.v[i][j][k];
                field_.w_curr[i][j][k] = field_.w[i][j][k];
            }
        }
    }

    return error;
}