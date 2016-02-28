#include<stdio.h>
#include"spectrum.h"
#include<string.h>
#include<stdlib.h>

double* get_continuum(  continuum_point cont_points[], 
                    const double waves[], 
                    size_t pts_size, size_t arr_size){
    /* continuum points must be pre-sorted by x*/

    double* continuum=(double*)malloc(arr_size*sizeof(double));
    //memcpy(continuum, waves);

    compute(1, pts_size-3, cont_points, continuum, waves, arr_size);
    return continuum;
}

void compute(int fpoint, int lpoint, continuum_point cont_points[], 
            double continuum[], const double waves[], 
            size_t arr_size){
    /*calculate continuum for given start and end point*/

    int bSplineNDiv, i,j,k;

    bSplineNDiv=50;
    //bSplineNDiv=200;

    for (i = fpoint; i <=lpoint; ++i) {  //iterate over cont points
          for (j = 0; j < bSplineNDiv; ++j) {   //iterate over bSplineNDiv
                double t0 = ((double) j) / bSplineNDiv;
                double t1 = ((double) (j + 1)) / bSplineNDiv;
                
                double x0 = 0.0;
                double y0 = 0.0;
                double x1 = 0.0;
                double y1 = 0.0;
                for (k = -1; k <= 2; ++k) {  //iterate over bases
                      continuum_point point = cont_points[i + k];  //get continuum points around current point (1 to the left, 2 to the right)
                      x0 += splineBasis(k, t0) * point.x;
                      x1 += splineBasis(k, t1) * point.x;
                      y0 += splineBasis(k, t0) * point.y;
                      y1 += splineBasis(k, t1) * point.y;
                }
                setContinLine(x0, x1, y0, y1, waves, continuum, arr_size);
          }
    }
}

void effectRegion(int i, continuum_point cont_points[]) {
    double t0 = 0;
    double t1 = 1;
    int k;
    double x0 = 0.0;
    double y0 = 0.0;
    double x1 = 0.0;
    double y1 = 0.0;
    for (k = -1; k <= 2; ++k) {
          continuum_point point = cont_points[i + k];
          x0 += splineBasis(k, t0) * point.x;
          x1 += splineBasis(k, t1) * point.x;
          y0 += splineBasis(k, t0) * point.y;
          y1 += splineBasis(k, t1) * point.y;
    }  
 
    return;
}

void update_region(int index, continuum_point cont_points[], 
            double continuum[], const double waves[], 
            size_t arr_size, size_t conts_pts_sz){
    
    if (index>((int)conts_pts_sz -4) || index<2){
        return;
    }
    compute(index-2, index+1, cont_points,continuum,waves,arr_size);
}

double splineBasis(int i, double t) {
    /*get the correct spline basis function*/
    switch (i) {
    case -1:
      return (((3.-t)*t - 3)*t + 1)/6.0;
    case 0:
      return (((3*t - 6)*t)*t + 4)/6.0;
    case 1:
      return (((3.-3*t)*t + 3)*t + 1)/6.0;
    case 2:
      return (t*t*t)/6.0;
    }
    

    fprintf(stderr,"spectrum/splineBasis.c was inappropriately called.  calculation has failed");
    return -1000000000.;
} 

void setContinLine(double x1, double x2, double y1, double y2, 
        const double waves[], double cont[], size_t wave_sz) {     

    /*set the continuum line*/

    int spix = get_wave_index(waves, wave_sz, x1)+1;
    int epix = get_wave_index(waves, wave_sz, x2);
    int j;
    double slope = (y2 - y1)/(x2 - x1);
    double intercept = y1 - x1*slope;
    
    for (j = spix; j <= epix; ++j) {
       cont[j] = waves[j]*slope + intercept;
    }
}

