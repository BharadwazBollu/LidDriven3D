#include "header.h"
#include <fstream>

void LidDrivenCavity::exportToVTK(int iter)
{
    std::ofstream fid;

    std::string filename = "LidDriven_" + std::to_string(iter) +  ".vtk";

    fid.open(filename);

    fid << "# vtk DataFile Version 2.0" << std::endl;
    fid << "Lid Driven Cavity 3D" << std::endl;
    fid << "ASCII" << std::endl;
    fid << "DATASET STRUCTURED_GRID" << std::endl;
    fid << "DIMENSIONS " << nx_ + 2 << " " << ny_ + 2 << " " << nz_ + 2 << std::endl;
    fid << "POINTS " << ( (nx_ + 2) * (ny_ + 2) * (nz_ + 2) ) << " double" << std::endl;

    // Write the point coordinates
    double z_ = 0.0;
    for (int k = 0; k <= nz_+1; k++) {
        if ( k==1){
        z_ = dz_/2;}
        if ( k==nz_+1){
        z_ = 1.0;}
        double y_ = 0.0;
        for (int j = 0; j <= ny_+1; j++) {
            if ( j==1){
            y_ = dy_/2;}
            if ( j==ny_+1){
            y_ = 1.0;}
            double x_ = 0.0;
            for (int i = 0; i <= nx_+1; i++) {
                if ( i==1){
                x_ = dx_/2;}
                if ( i==nx_+1){
                x_ = 1.0;}
                fid << x_ << " " << y_ << " " << z_ << std::endl;
                x_ = x_ + dx_;
            }
        y_ = y_ + dy_;    
        }
    z_ = z_ + dz_;    
    }

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

    fid << "SCALARS " << "P" << " double 1 " << std::endl;

    fid << "LOOKUP_TABLE my_table" << std::endl;
    for (int k = 0; k <= nz_+1; k++) {
        for (int j = 0; j <= ny_+1; j++) {
            for (int i = 0; i <= nx_+1; i++) {
                fid << field_.p[i][j][k] << std::endl;
            }
        }
    }

    fid.close();

}