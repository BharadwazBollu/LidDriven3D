#include "header.h"

LidDrivenCavity::LidDrivenCavity(int nx, int ny, int nz, double lid_velocity,
double mu, double density, double dt, double x_length, double y_length, 
double z_length, double prediction_tolerance, double pressure_tolerance,
double correction_tolerance)
:nx_(nx), ny_(ny), nz_(nz), lid_velocity_(lid_velocity),
mu_(mu), density_(density),  dt_(dt), x_length_(x_length), 
y_length_(y_length), z_length_(z_length), 
prediction_tolerance_(prediction_tolerance),
pressure_tolerance_(pressure_tolerance),
correction_tolerance_(correction_tolerance)
{
    // Giving size for field variable as per input including ghost cells
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



