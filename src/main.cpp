#include "header.h"

int main(){

int iteration = 1;
double err = 1.0;
Fields fields(30, 30, 30, 0.01, 1.0, 0.001);
fields.initialise();
//fields.exportToVTK(iteration);

while ( err > 1e-3)
{

    
    fields.solvePrediction();
    fields.solvePressure();
    fields.solveCorrection();

    err = fields.calculateError();

    std::cout << " Iteration " << iteration << " error = " << err << std::endl;
    iteration++;
    

    if ( iteration % 100 == 0)
    {
        fields.exportToVTK(iteration);
    }

}

return 0;
}