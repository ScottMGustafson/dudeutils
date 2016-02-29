#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"spectrum.h"

double* spectrum(   const double waves[], const double flux[], const double err[], double continuum[],
                    double xs[], double ys[],
                    double Nval[], double bval[], double redshift[], 
                    double restwv[], double gmma[], double osc_strn[],
                    int starts[], int ends[],
                    int len_cont_points, int len_arr, int len_abs, int len_pairs){
    /*
    int i=0;
    continuum_point cont_points[len_cont_points];
    absorber absorbers[len_abs];
    index_pair pairs[len_pairs];

    //double continuum[len_arr];


    for(i=0;i<len_pairs;++i){
        pairs[i]=(index_pair){.start=starts[i], .end=ends[i]};
    }
    

    
    memcpy(absorption, continuum, sizeof(continuum));

    for(i=0;i<len_abs;i++){
        absorbers[i] = (absorber) { .N=Nval[i], .b=bval[i], .z=redshift[i], 
                                    .rest=restwv[i], .gamma=gmma[i], .f=osc_strn[i]};
    }
    double* cont = (double*)get_continuum( cont_points, waves, len_cont_points, len_arr); 
    //double* absorption =(double*)get_absorption(cont, waves, absorbers,len_abs,len_arr);
    double* absorption =(double*)get_absorption(continuum, waves, absorbers,len_abs,len_arr);

    free(cont);
    cont*=NULL;
    //spec spectrum= (spec) {.waves=waves, .flux=flux, 
    //                        .absorption=absorbed_flux, .continuum=continuum}; 
    //double chi2 = chi2(absorption, flux, err, pairs, len_pairs);
    */
    return continuum;
}

double* return_absorption(double continuum[], const double waves[], double Nval[], double bval[], double redshift[], 
                    double restwv[], double gmma[], double osc_strn[], size_t len_abs, size_t len_arr){
    int i;
    absorber absorbers[len_abs];

    for(i=0;i<len_abs;i++){
        absorbers[i] = (absorber) { .N=Nval[i], .b=bval[i], .z=redshift[i], 
                                    .rest=restwv[i], .gamma=gmma[i], .f=osc_strn[i]};
    }

    double* absorption =(double*)get_absorption(continuum, waves, absorbers,len_abs,len_arr);
    //memcpy(absorption, continuum, sizeof(continuum));


    return absorption;
}

double* return_continuum(double xs[], double ys[], const double waves[], size_t len_cont_points, size_t len_arr){
    continuum_point cont_points[len_cont_points];
    int i;
    for(i=0;i<len_cont_points;++i){
        cont_points[i]=(continuum_point){.x=xs[i], .y=ys[i]};
    }
    double* cont =(double*)get_continuum(cont_points, waves, len_cont_points, len_arr);
    return cont;
}

double get_chi2(double absorption[], const double flux[], const double err[], const double waves[], double starts[], double ends[], size_t len_pairs, size_t len_arr){
    int i;
    int beg, end;

    index_pair pairs[len_pairs];

    for(i=0;i<len_pairs;++i){
        beg = get_wave_index(waves, len_arr, starts[i]);
        end = get_wave_index(waves, len_arr, ends[i]);
        pairs[i]=(index_pair){.start=beg, .end=end};
    }
    return chi2(absorption, flux, err, pairs, len_pairs);
}

