#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"spectrum.h"

const double c = 2.99792458E5;
const double pi=3.14159265359;
  
int sub_absorber(   const double waves[], double cont[],
                    double N, double b, double z, 
                    double f, double Gamma, double restWave, 
                    bool subFlag, size_t size){

    /*
    subtract absorption from continuum level given some absorber
    caller is responsible for initializing flux

    */

    double tau_threshold = 0.001;

    double cwave = (1.0+z)*restWave;
    double vdopp = b/cwave*1.0e13;
    double alpha = Gamma/(4.0*pi*vdopp)/(1.0+z);    
    double factor = pow(10.0, N) * 2.647E-2 * f/(sqrt(pi)*vdopp) * 1.0/(1.0+z);
    int i=0;
    int cpix = get_wave_index(waves, size, cwave);

    if (cpix >= size-3){
        fprintf(stderr,"cext/get_absorption.c: center pix too big: %d out of %d pixels for %5.1lf\n",(int)cpix,(int)size,cwave);
        return 0;
    }
    if (cpix < 1){
        fprintf(stderr,"cext/get_absorption.c: center pix too small: %d out of %d pixelsfor %5.1lf\n",(int)cpix,(int)size,cwave);
        return 0;
    }


    double tau;
    double wave_high = (waves[cpix] + waves[cpix+1])/2.;
    double wave_low  = (waves[cpix] + waves[cpix-1])/2.;

    int nsamp = 10;
    double delta_wave = (wave_high - wave_low)/nsamp;

    double tau_sum = 0.0;

    for (i=0; i<nsamp;  ++i) {
        double wave_i = wave_low + i*delta_wave;
        double vbar = c/b * (wave_i/cwave - 1.0);
        tau_sum += factor*voigt(vbar, alpha);
    }

    tau = tau_sum / nsamp;

    //if (subFlag) cont[cpix] -= tau;
    //else         cont[cpix] += tau;
    if (subFlag) cont[cpix] /= exp(tau);
    else         cont[cpix] *= exp(tau);


    int lpix = cpix;
    while (tau > tau_threshold && lpix > 1) {
        lpix -= 1;
        double wave_high = (waves[lpix] + waves[lpix+1])/2;
        double wave_low  = (waves[lpix] + waves[lpix-1])/2;

        double vbar = c/b * (wave_low/cwave - 1.0);
        double tau_low = factor*voigt(vbar, alpha);

        vbar = c/b * (wave_high/cwave - 1.0);
        double tau_high = factor*voigt(vbar, alpha);

        tau = (tau_low + tau_high)/2;

        //if (subFlag) cont[lpix] -= tau;
        //else         cont[lpix] += tau;

        if (subFlag) cont[lpix] /= exp(tau);
        else         cont[lpix] *= exp(tau);
    }

    //
    // Now walk the other way ...
    int hpix = cpix;

    double vbar = c/b * (waves[hpix]/cwave - 1.0);
    tau = factor*voigt(vbar, alpha);

    while (tau > tau_threshold && hpix < size-3) {
        hpix += 1;

        double wave_high = (waves[hpix] + waves[hpix+1])/2;
        double wave_low  = (waves[hpix] + waves[hpix-1])/2;

        double vbar = c/b * (wave_low/cwave - 1.0);
        double tau_low = factor*voigt(vbar, alpha);

        vbar = c/b * (wave_high/cwave - 1.0);
        double tau_high = factor*voigt(vbar, alpha);

        tau = (tau_low + tau_high)/2;


        //if (subFlag) cont[hpix] -= tau;
        //else         cont[hpix] += tau;

        if (subFlag) cont[hpix] /= exp(tau);
        else         cont[hpix] *= exp(tau);
    }

    return 0;
}

double* get_absorption(double cont[], const double waves[], absorber abs[], 
                        size_t num_absorbers, size_t len_arr){
    /*inputs are arrays of various absorber attributes*/

    int i;

    double* absorption = (double*)malloc(len_arr*sizeof(double));
    for(i=0;i<len_arr;++i){
        absorption[i]=cont[i]; /*initialie absorber flux as continuum level*/
    }

    

    for (i=0;i<num_absorbers;i++){ // subtract absorption
        sub_absorber(  waves, absorption, 
                            abs[i].N,abs[i].b, abs[i].z, 
                            abs[i].f, abs[i].gamma, abs[i].rest, 
                            true, len_arr);
    }
    return absorption;
}
