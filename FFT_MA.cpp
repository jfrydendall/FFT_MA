/* 
 * File:   FFT_MA.cpp
 * Author: jf
 * 
 * Created on January 15, 2013, 1:26 PM
 */

#include "FFT_MA.hpp"

FFT_MA::FFT_MA(int Nrel, int Nx, int Ny, int Nz, double cell,
        double CorLenMax, double CorLenMed, double CorLenMin, double theta,
        double beta, double alpha, double mean, double var, std::string CorrFunc) {
    m_nrel = Nrel;
    m_corLenMax = CorLenMax/cell;
    m_corLenMed = CorLenMed/cell;
    m_corLenMin = CorLenMin/cell;
    m_cell = cell;
    m_mean = mean;
    m_std = std::sqrt(var);
    m_Nx = Nx;
    m_Ny = Ny;
    m_Nz = Nz;
    m_l = 2 * 10;
    m_nxx = Nx + m_l;
    m_nyy = Ny + m_l;
    m_nzz = Nz + m_l;
    m_nzp = m_nzz / 2 + 1;
    m_scale = 1.0 / (m_nxx * m_nyy * m_nzz);
    m_align = sizeof (Complex);
    m_size = m_nxx*m_nyy*m_nzz;
    name = CorrFunc;
    
    if (m_corLenMin > m_corLenMax) {
        std::cout << "Warning:" << std::endl;
        std::cout << "Choose m_corLenMin<=m_corLenMax because m_corLenMax is the direction of maximum continuity." << std::endl;
        std::cout << "Use the input 'ang' to change the direction of maximum continuity" << std::endl;
    }

    m_gamma = m_corLenMin / m_corLenMax; // anistrophy factor < 1 (=1 for isotropy)
    m_delta = m_corLenMin / m_corLenMed;
    m_theta = pi / 180 * theta; // Transform angle into radians
    m_beta = pi / 180 * beta; // Transform angle into radians
    m_alpha = pi / 180 * alpha; // Transform angle into radians
    
    
    pi = atan(1)*4;

    // internal default values - these values can be modified with the set functions
    m_writeFile = "FFT_MA.bin";
    m_setNoiseFlag = false;
    m_setCorrelationFlag = false;
    fftw::maxthreads = omp_get_max_threads();
    //set seed
    set_seed(123,456);
    
}

void FFT_MA::setNoise(array4<double>& Z){
    this->m_Z = Z;
}

void FFT_MA::setCorrelation(array3<double>& C){
    this->m_C = C;
}

void FFT_MA::setSeed(int seed[2]){
    set_seed(seed[0],seed[1]);
}

void FFT_MA::setWriteFile(string writefile){
    this->m_writeFile = writefile;
}


void FFT_MA::distance(array3<double>& r) {
    double Xrot, Yrot, Zrot, xx, yy, zz;
    
    for (int i = 0; i<this->m_nxx; i++) {
        for (int j = 0; j<this->m_nyy; j++) {
            for (int k= 0; k<this->m_nzz; k++) {
                xx = i + this->m_cell / 2.0  - floor(this->m_nxx/2);
                yy = j + this->m_cell / 2.0  - floor(this->m_nyy/2);
                zz = k + this->m_cell / 2.0  - floor(this->m_nzz/2);
                
                Xrot = (cos(this->m_theta) * cos(this->m_alpha) - sin(this->m_theta) * sin(this->m_beta) * sin(this->m_alpha)) * xx +
                        (cos(this->m_theta) * sin(this->m_alpha) + sin(this->m_theta) * sin(this->m_beta) * cos(this->m_alpha)) * yy +
                        (-sin(this->m_theta) * cos(this->m_beta)) * zz;

                Yrot = (-cos(this->m_beta) * sin(this->m_alpha)) * xx +
                        (cos(this->m_beta) * cos(this->m_alpha)) * yy +
                        (sin(this->m_beta)) * zz;

                Zrot = (sin(this->m_theta) * cos(this->m_alpha) + cos(this->m_theta) * sin(this->m_beta) * sin(this->m_alpha)) * xx +
                        (sin(this->m_theta) * sin(this->m_alpha) - cos(this->m_theta) * sin(this->m_beta) * cos(this->m_alpha)) * yy +
                        (cos(this->m_theta) * cos(this->m_beta)) * zz;
                
                r(i, j, k) = sqrt(pow(m_delta * Xrot, 2.0) + pow(m_gamma * Yrot, 2.0) + pow(Zrot, 2.0));
            }
        }
    }
}

void FFT_MA::correlation_functions(array3<double>& C) {
    
    array3<double> r(this->m_nxx, this->m_nyy, this->m_nzz, this->m_align);
    FFT_MA::distance(r);
    if (this->name == "Spherical") {
        for (int i = 0; i < this->m_size; i++) {
            if (r(i)>this->m_corLenMin)
                r(i) = this->m_corLenMin;
            C(i) = 1.0 - (1.5 * (r(i) / this->m_corLenMin) - 0.5 * pow(r(i) / this->m_corLenMin, 3));
        }

    } else if (this->name == "Exponential") {
        for (int i = 0; i < this->m_size; i++)
            C(i) = exp(-3.0 * r(i) / this->m_corLenMin);
    }
    else if (this->name == "Gaussian") {
        for (int i = 0; i < this->m_size; i++){
            C(i) = exp(-3.0 * pow(r(i), 2.0) / pow(this->m_corLenMin, 2.0));
        }
    }
    r.Deallocate();
}

void FFT_MA::normal(array4<double>& Z){

    for (int i = 0; i < this->m_size*this->m_nrel; i++)
        Z(i) = rnorm(0,1);
}

void FFT_MA::randomfield(double** fields) {

    array3<Complex> fC(this->m_nxx, this->m_nyy, this->m_nzp, this->m_align);
    array3<Complex> fZ(this->m_nxx, this->m_nyy, this->m_nzp, this->m_align);
    array3<Complex> fR(this->m_nxx, this->m_nyy, this->m_nzp, this->m_align);

    array3<double> C(this->m_nxx, this->m_nyy, this->m_nzz, this->m_align);
    array4<double> Z(this->m_nrel,this->m_nxx, this->m_nyy, this->m_nzz, this->m_align);

    array3<double> phi(this->m_nxx, this->m_nyy, this->m_nzz, this->m_align);

    rcfft3d Forward3C(this->m_nzz, C, fC);
    rcfft3d Forward3Z(this->m_nzz, Z[0], fZ);
    crfft3d Backward3(this->m_nzz, fR, phi);

    int size = this->m_nxx * this->m_nyy * this->m_nzp;

    if(m_setCorrelationFlag)
        C = this->m_C;
    else
        this->correlation_functions(C);
    
    if(m_setNoiseFlag)
        Z = this->m_Z;
    else
        this->normal(Z);
    
    Forward3C.fft(C, fC);
    double t1 = omp_get_wtime();
    for (int rel = 0; rel < m_nrel; rel++) {
        cout<<"realization:  "<<rel<<endl;
        Forward3Z.fft(Z[rel], fZ);
        for (int i = 0; i < size; i++)
            fR(i) = this->m_std * sqrt(fC(i)) * fZ(i);

        Backward3.fftNormalized(fR, phi);

        // slice array down to size : nx, ny, nz
        
        this->slice(phi,fields[rel]);
    }
    double t2 = omp_get_wtime();
    cout<<"Total time of loop:  "<<t2-t1<<endl;
    C.Deallocate(); fC.Deallocate();
    fZ.Deallocate();fR.Deallocate();
    phi.Deallocate(); Z.Deallocate();

}

double** FFT_MA::fields(int rows, int cols){
    double** matrix;
    matrix = new double* [rows];
    for (int i = 0;i<rows;i++)
        matrix[i] = new double [cols];
    return matrix;
}

void FFT_MA::freeFields(double** matrix, int rows) {
    for (int i = 0; i < rows; i++)
        delete [] matrix[i];
    delete [] matrix;
}

void FFT_MA::writeFields(double** matrix, int rows){
    FILE* pFile;
    int size = this->m_Nx * this->m_Ny * this->m_Nz;
    pFile = fopen(this->m_writeFile.c_str(), "wb");
    for (int rel = 0; rel < rows; ++rel){
        //Some calculations to fill a[]
        fwrite(matrix[rel], 1, size*sizeof(double), pFile);
    }
    fclose(pFile);
}

void FFT_MA::slice(array3<double>& A,double* values) {
    // slice array down to size : nx, ny, nz
    int nyz = this->m_Ny*m_Nz;
    int nz = this->m_Nz;
    
    for (int i = 0; i<this->m_Nx; i++) {
        for (int j = 0; j<this->m_Ny; j++) {
            for (int k = 0; k<this->m_Nz; k++) {
                values[i * nyz + j * nz + k] = A(i, j, k) + this->m_mean;
            }
        }
    }
}


FFT_MA::~FFT_MA() {
}

