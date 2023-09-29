#include "header.h"

int main(){

int iteration = 1;
double err = 1.0;
Fields fields(20, 20, 20, 0.01, 1.0, 0.001);
fields.initialise();

while ( err > 1e-3)
{

    fields.solvePrediction();
    fields.solvePressure();
    fields.solveCorrection();

    err = fields.calculateError();

    std::cout << " Iteration " << iteration << " error = " << err << std::endl;
    iteration++;

    if ( iteration % 1 == 0)
    {
        fields.exportToVTK(iteration);
    }

}




return 0;
}