#ifndef header_h
#define header_h

#include <algorithm>
#include <iostream>
#include <math.h>
#include <vector>

class LidDrivenCavity{
    public:

        LidDrivenCavity(int nx, int ny, int nz, double lid_velocity,
        double mu, double density, double dt, double x_length, double y_length, 
        double z_length, double prediction_tolerance, double pressure_tolerance,
        double correction_tolerance, int write_interval);

        void initialize();
        void solvePrediction();
        void solvePressure();
        void solveCorrection();
        double calculateError();
        void exportToVTK(int iteration);
        
    private:
        // Number of cells in x, y, z
        int nx_, ny_, nz_;
        // Top moving lid velocity
        double lid_velocity_;
        // dynamic viscosity, density, time step
        double mu_, density_, dt_;
        // length in x, y, z
        double x_length_, y_length_, z_length_;
        // tolerances for prediction, pressure, correction
        double prediction_tolerance_, pressure_tolerance_, correction_tolerance_;
        // write interval to export results
        int write_interval_;
        // local iteration count for u, v, w, p
        int u_local_iteration_, v_local_iteration_, w_local_iteration_, p_local_iteartion_;
        // root means square for u, v, w, p
        double u_rms_, v_rms_, w_rms_, p_rms_;
        // residual for u, v, w, p
        double u_residual_, v_residual_, w_residual_, p_residual_;
        // error for gauss seidal
        double error_;
        // central coefficient for prediction and correction 
        double central_coefficient_;
        // interior total count
        int cell_count_ = nx_ * ny_ * nz_;
        // Cell dimensions in x, y, z
        double dx_ = x_length_/nx_;
        double dy_ = y_length_/ny_;
        double dz_ = z_length_/nz_;
        // volume of the cell
        double Vp_ = dx_ * dy_ * dz_;
        // areas in x, y, z
        double x_area_ = dy_ * dz_;
        double y_area_ = dx_ * dz_;
        double z_area_ = dx_ * dy_;
        // Flux in all directions
        double fluxE_, fluxW_, fluxN_, fluxS_, fluxF_, fluxB_;
        // Upwind velocities declartion
        double uE_, uW_, uN_, uS_, uF_, uB_;
        double vE_, vW_, vN_, vS_, vF_, vB_;
        double wE_, wW_, wN_, wS_, wF_, wB_;
        // Total cell convection in x, y, z
        double x_convection_, y_convection_, z_convection_;
        // Convection diagonal coefficient for cell
        double conv_diag_;
        // Total cell diffusion in x, y, z, p
        double x_diffusion_, y_diffusion_, z_diffusion_, p_diffusion_;
        // Diffusion descretization constants ( to save repereated computation )
        double x_diffusion_coefficient_ = x_area_ / dx_;
        double y_diffusion_coefficient_ = y_area_ / dy_;
        double z_diffusion_coefficient_ = z_area_ / dz_;
        // Diffusion central coefficient
        double diff_cen_coff = 2 * ( x_diffusion_coefficient_ + y_diffusion_coefficient_ + z_diffusion_coefficient_ );
        // total flux per cell with no pressure correction
        double flux_no_pressure_correction_;
        // declaring field varaibles 
        struct data
        {
            std::vector<std::vector<std::vector<double>>> u;
            std::vector<std::vector<std::vector<double>>> v;
            std::vector<std::vector<std::vector<double>>> w;
            std::vector<std::vector<std::vector<double>>> u_pred;
            std::vector<std::vector<std::vector<double>>> v_pred;
            std::vector<std::vector<std::vector<double>>> w_pred;
            std::vector<std::vector<std::vector<double>>> u_curr;
            std::vector<std::vector<std::vector<double>>> v_curr;
            std::vector<std::vector<std::vector<double>>> w_curr;
            std::vector<std::vector<std::vector<double>>> p;
        };
        data field_;
};
// git reset --hard HEAD
#endif
