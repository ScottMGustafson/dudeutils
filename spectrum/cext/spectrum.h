#ifndef HEADER_GUARD
#define HEADER_GUARD

/*typedef struct Spectrum{
    //size_t size;
    double wave;
    double continuum;
    double absorption;
    double flux;
    //double error[size];
} spec;*/

typedef enum { false, true } bool;

typedef struct Absorber{
    double N;
    double b;
    double z;
    double rest;
    double gamma;
    double f;
} absorber;

typedef struct ContinuumPoint{
    double x;
    double y;
} continuum_point;

typedef struct IndexPair{
    int start;
    int end;
} index_pair;

/*spec get_spectrum(double waves[], double flux[], 
        continuum_point cont_points[], absorber absorbers[], 
        size_t num_abs, size_t size, size_t pts_size);*/

double* spectrum(   const double waves[], const double flux[], const double err[], double continuum[],
                    double xs[], double ys[],
                    double Nval[], double bval[], double redshift[], 
                    double restwv[], double gmma[], double osc_strn[],
                    int starts[], int ends[],
                    int len_cont_points, int len_arr, int len_abs, int len_pairs);

double* return_absorption(double continuum[], const double waves[], double Nval[], double bval[], double redshift[], 
                    double restwv[], double gmma[], double osc_strn[], size_t len_abs, size_t len_arr);

double* return_continuum(double xs[], double ys[], const double waves[], size_t len_cont_points, size_t len_arr);

double get_chi2(double absorption[], const double flux[], const double err[], const double waves[],double starts[], double ends[], size_t len_pairs, size_t len_arr);
/************************************************************/

double voigt(double v, double a);

/*************************************************************/

double pix_to_wave(double crval, double cdelt, double crpix, double pix, bool loglin);

int get_wave_index(const double wave_arr[], size_t size, double wave);

double chi2(double mod[], const double obs[], const double err[], index_pair pairs[], size_t len_pairs);

/**************************************************************/

int sub_absorber(   const double waves[], double cont[],
                    double N, double b, double z, 
                    double f, double Gamma, double restWave, 
                    bool subFlag, size_t size);

double* get_absorption(double cont[], const double waves[], absorber abs[], 
                        size_t num_absorbers, size_t len_arr);


/************************************************************/

double* get_continuum(  continuum_point cont_points[], 
                    const double waves[], 
                    size_t pts_size, size_t arr_size);

void compute(int fpoint, int lpoint, continuum_point cont_points[], 
            double continuum[], const double waves[], 
            size_t arr_size);

void effectRegion(int i, continuum_point cont_points[]);

void update_region(int index, continuum_point cont_points[], 
            double continuum[], const double waves[], 
            size_t arr_size, size_t conts_pts_sz);

double splineBasis(int i, double t);

void setContinLine(double x1, double x2, double y1, double y2, 
        const double waves[], double cont[], size_t wave_sz);

#endif

