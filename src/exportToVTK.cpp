#include "header.h"
#include <fstream>

void Fields::exportToVTK(int iter)
{
    std::ofstream fid;

    fid.open("LidDriven.vtk");

    fid << "# vtk DataFile Version 2.0" << std::endl;
    fid << "Cube example" << std::endl;
    fid << "ASCII" << std::endl;
    fid << "DATASET STRUCTURED_POINTS" << std::endl;
    fid << "DIMENSIONS " << nx_ + 2 << " " << ny_ + 2 << " " << nz_ + 2 << std::endl;
    fid << "ASPECT_RATIO " << 1 << " " << 1 << " " << 1 << std::endl;
    fid << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
    fid << "POINT_DATA " << (nx_ + 2) * (ny_ + 2) * (nz_ + 2) << std::endl;
    fid << "SCALARS " << "U" << " double 1 " << std::endl;

    fid << "LOOKUP_TABLE my_table" << std::endl;
    for (int k = 0; k <= nz_+1; k++) {
        for (int j = 0; j <= ny_+1; j++) {
            for (int i = 0; i <= nx_+1; i++) {
                fid << field_.u[i][j][k] << std::endl;
            }
        }
    }

    fid << "SCALARS " << "V" << " double 1 " << std::endl;

    fid << "LOOKUP_TABLE my_table" << std::endl;
    for (int k = 0; k <= nz_+1; k++) {
        for (int j = 0; j <= ny_+1; j++) {
            for (int i = 0; i <= nx_+1; i++) {
                fid << field_.v[i][j][k] << std::endl;
            }
        }
    }

    fid << "SCALARS " << "W" << " double 1 " << std::endl;

    fid << "LOOKUP_TABLE my_table" << std::endl;
    for (int k = 0; k <= nz_+1; k++) {
        for (int j = 0; j <= ny_+1; j++) {
            for (int i = 0; i <= nx_+1; i++) {
                fid << field_.w[i][j][k] << std::endl;
            }
        }
    }

    fid.close();

}