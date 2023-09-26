#ifndef header_h
#define header_h

#include <math.h>
#include <vector>
#include <iostream>

void hello();

class Fields{
    public:

        Fields(int nx, int ny, int nz);

        void initialise();
        void solvePrediction();
        void solvePressure();
        void solveCorrection();
        
    private:
        int nx_, ny_, nz_;

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

        data fields;


};

#endif