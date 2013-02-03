/* 
 * File:   main.cpp
 * Author: jf
 *
 * Created on January 15, 2013, 1:30 PM
 */

#include <cstdlib>

#include "FFT_MA.hpp"

using namespace std;

void print_row(double arg[], int cols){
    for (int j=0;j<cols;j++)
        cout << arg[j]<< " ";
    cout << "\n";
}

void print_matrix(double** arg, int rows, int cols)
{
        for (int i=0; i<rows; i++)
            print_row(arg[i], cols);
            
}

/*
 * 
 */
int main(int argc, char** argv) {
    int nx = 120, ny = 80, nz = 1, Nrel = 100;
    double cell = 1.0, CorLenMax = 15.0, CorLenMed = 10.0, CorLenMin = 5.0;
    double theta = 0.0, beta = 0.0, alpha = 45.0;
    double mean = 0.0, var = 1.0;
    std::string name = "Spherical";
    
    FFT_MA rf(Nrel, nx, ny, nz, cell,CorLenMax, CorLenMed, CorLenMin, 
              theta, beta, alpha, mean, var, name);
    
// Change default output file to a user-defined:
//    string S  = "test.bin";
//    rf.setWriteFile(S);
    
    
//  test 1 - random fields    
// Change default seed file to a user-defined:
    int seed[2] = {1234,5678};
    rf.setSeed(seed);    
        
    double** randomfields = rf.fields(Nrel, nx*ny*nz);
    rf.randomfield(randomfields);
    
    rf.writeFields(randomfields,Nrel);
    
//    print_matrix(randomfields,Nrel, nx*ny*nz);
    
    rf.freeFields(randomfields,Nrel);
    
//  test 2 - co-simulated random fields
    
    // remember to allocate for aliasing effects, i.e. 20 to each direction
    array4<double> Z0,Z1,Z2;
    Z0.Allocate(Nrel,nx+20,ny+20,nz+20);
    Z1.Allocate(Nrel,nx+20,ny+20,nz+20);
    Z2.Allocate(Nrel,nx+20,ny+20,nz+20);
    
    rf.setSeed(seed);
    rf.normal(Z0);
    rf.normal(Z1);
    // first get the normal fields
        
    rf.setNoise(Z0);
        
    randomfields = rf.fields(Nrel, nx*ny*nz);
    rf.randomfield(randomfields);
    
    string S0 = "Z0.bin";
    rf.setWriteFile(S0);
    
    rf.writeFields(randomfields,Nrel);
    
    rf.setNoise(Z1);
    
    // the second normal field
        
    randomfields = rf.fields(Nrel, nx*ny*nz);
    rf.randomfield(randomfields);
    
    string S1 = "Z1.bin";
    rf.setWriteFile(S1);
    
    rf.writeFields(randomfields,Nrel);
    
    
    // now use the known noise to generate co-simulated fields
    // we choose rho to 0.8
        
    double rho = 0.8;
    
    for (unsigned int i=0;i<Z2.Size();i++)
        Z2(i) = rho*Z1(i) + sqrt(1.0 - pow(rho,2))*Z0(i);
    
    rf.setNoise(Z2);
    rf.randomfield(randomfields);
    
    string S2 = "cosim.bin";
    rf.setWriteFile(S2);
    
    rf.writeFields(randomfields,Nrel);
    
    rf.freeFields(randomfields,Nrel);
    
    Z0.Deallocate();Z1.Deallocate();Z2.Deallocate();
        
    return 0;
}

