#include "header.h"

int main(){

int iteration = 1;
int write_interval = 500;
double err = 1.0;
LidDrivenCavity fields(50, 50, 50, 1.0, 0.01, 1.0, 0.001, 1.0, 1.0, 1.0, 1e-8, 1e-4, 1e-8);
fields.initialize();

while ( err > 1e-3)
{
    fields.solvePrediction();
    fields.solvePressure();
    fields.solveCorrection();

    err = fields.calculateError();

    std::cout << " Iteration " << iteration << " error = " << err << std::endl;
    iteration++;
    
    if ( iteration % write_interval == 0)
    {
        fields.exportToVTK(iteration);
    }
}

fields.exportToVTK(iteration);

return 0;
}