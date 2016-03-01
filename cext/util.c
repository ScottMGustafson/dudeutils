#include<stdio.h>
#include<math.h>
#include"spectrum.h"

double pix_to_wave(double crval, double cdelt, double crpix, double pix, bool loglin){
    if(loglin){
        return pow(10.,(crval+cdelt*(pix-crpix)));
    }else{
        return crval+cdelt*(pix-crpix);
    }
}

int get_wave_index(const double wave_arr[], size_t size, double wave){
    /*assumes that double[] waves is sorted alread*/
    int i=0;  
    double vdisp=2.14  //velocity dispersion for HIRES in kms
    double c=2.99792458E5 //c in km/s

    //double tol=0.02;

    //with a loglin scale:
    double tol=wave*vdisp/c;   //one pixel in wavlength

    while ( i < size && fabs(wave_arr[i]-wave)>tol ){
        //printf("diff=%lf \n",fabs(wave_arr[i]-wave));
        ++i;
    }

    return ( i == size ? -1 : i );
}

double chi2(double mod[], const double obs[], const double err[], index_pair pairs[], size_t len_pairs) {
    int i,k;
    double result = 0.0;
    double diff=0.0;
    for (k=0;k<len_pairs;++k){
        
        for (i = pairs[k].start; i < pairs[k].end; ++i) {
            diff = (obs[i] - mod[i])/err[i];
            result += diff*diff;
        }
    }
    return result;
}
