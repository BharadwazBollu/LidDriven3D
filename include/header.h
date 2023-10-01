#ifndef header_h
#define header_h

#include <math.h>
#include <vector>
#include <iostream>

class Fields{
    public:

        Fields(int nx, int ny, int nz, double nu, double density, double dt);

        void initialise();
        void solvePrediction();
        void solvePressure();
        void solveCorrection();
        double calculateError();
        void exportToVTK(int iter);
        
    private:
        int nx_, ny_, nz_;
        double nu_, density_, dt_;

        // Cell dimensions in x, y, z
        double dx_ = 1.0/nx_;
        double dy_ = 1.0/ny_;
        double dz_ = 1.0/nz_;

        // Volume of the cell
        double Vp_ = dx_ * dy_ * dz_;

        // Areas in x, y, z
        double areaX_ = dy_ * dz_;
        double areaY_ = dx_ * dz_;
        double areaZ_ = dx_ * dy_;

        // Flux in all directions
        double fluxE_, fluxW_, fluxN_, fluxS_, fluxF_, fluxB_;

        // Upwind velocities declartion
        double uE_, uW_, uN_, uS_, uF_, uB_;
        double vE_, vW_, vN_, vS_, vF_, vB_;
        double wE_, wW_, wN_, wS_, wF_, wB_;

        // Total cell Convection in x, y, z
        double convectionX_, convectionY_, convectionZ_;

        // Convection diagonal constant for cell
        double conv_diag_;

        // Total cell Diffusion in x, y, z
        double diffusionX_, diffusionY_, diffusionZ_;

        // Diffusion descretization denominators constants ( to save repereated computation )
        double diff_disc_coffX_ = areaX_ / dx_;
        double diff_disc_coffY_ = areaY_ / dy_;
        double diff_disc_coffZ_ = areaZ_ / dz_;

        // Diffusion central coefficient
        double diff_cen_coff = 2 * ( diff_disc_coffX_ + diff_disc_coffY_ + diff_disc_coffZ_ );

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
