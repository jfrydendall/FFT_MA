/* 
 * File:   FFT_MA.hpp
 * Author: jf
 *
 * Created on January 15, 2013, 1:26 PM
 */

#ifndef FFT_MA_HPP
#define	FFT_MA_HPP


#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <stdio.h>
#include "Array.h"
#include "fftw++.h"
#define MATHLIB_STANDALONE
#include <Rmath.h>

using namespace std;
using namespace Array;
using namespace fftwpp;

class FFT_MA {
    
public:
    FFT_MA(int, int,int,int,double,double,double,double,double,double,double,double,double,std::string);
    void distance(array3<double>&);
    void correlation_functions(array3<double>&);
    void normal(array4<double>&);
    void randomfield(double**);
    void writeFields(double**,int);
    double** fields(int,int);
    void freeFields(double**,int);
    void setNoise(array4<double>&);
    void setCorrelation(array3<double>&);
    void setSeed(int []);
    void setWriteFile(string);
    virtual ~FFT_MA();
private:
    int m_size;
    int m_nxx, m_nyy, m_nzz, m_nzp;
    size_t m_align;
    int  m_l;
    double m_delta, m_gamma, m_scale, m_std;
    int m_nrel;
    int m_Nx, m_Ny, m_Nz;
    double m_cell, m_corLenMax, m_corLenMed, m_corLenMin, m_theta, m_beta, m_alpha;
    double m_mean, m_var;
    array3<double> m_Z,m_C;
    bool m_setNoiseFlag, m_setCorrelationFlag;
    std::string name, m_writeFile;
    void slice(array3<double>&, double*);
    double pi;
};

#endif	/* FFT_MA_HPP */

