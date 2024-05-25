#include<cstdio>


double fabs(double x);
void vector1_equality_vector2(int dimension,double* vector1,const double* vector2);
double Max(double V1,double V2);
double Min(double V1,double V2);

double RKDOPR_8_7(
    int dimension,
    double t,
    const double* x_start,
    double step,
    double* x_end,
    double ( *(*F) )(double ,const double*),
    double* tmp_memory
);

double IntersectionSearch_with_increases_up(
int dimension,
double* sistem_parametr,
double ( *(*F) )(double ,const double*),
double* tmp_memory
);

double IntersectionSearch_with_increases_bot(
int dimension,
double* sistem_parametr,
double ( *(*F) )(double ,const double*),
double* tmp_memory
);

double IntersectionSearch_with_decreases_up(
int dimension,
double* sistem_parametr,
double ( *(*F) )(double ,const double*),
double* tmp_memory
);

double IntersectionSearch_with_decreases_bot(
int dimension,
double* sistem_parametr,
double ( *(*F) )(double ,const double*),
double* tmp_memory
);


double IntersectionCheck_with_top_up(
int dimension,
double* sistem_parametr,
double ( *(*F) )(double ,const double*),
double* tmp_memory
);

double IntersectionCheck_with_top_bot(
int dimension,
double* sistem_parametr,
double ( *(*F) )(double ,const double*),
double* tmp_memory
);


double IntersectionCheck_with_bottom_up(
int dimension,
double* sistem_parametr,
double ( *(*F) )(double ,const double*),
double* tmp_memory
);

double IntersectionCheck_with_bottom_bot(
int dimension,
double* sistem_parametr,
double ( *(*F) )(double ,const double*),
double* tmp_memory
);


void shift(
int dimension,
double* sistem_parametr
);


void shift_print(
int dimension,
double* sistem_parametr
);


double F_x_main(double t,const double* x);

double F_y_main(double t, const double* x);

double F_px_main(double t, const double* x);

double F_py_main(double t, const double* x);

double RKDOPR_8_7_interval_with_control(
const int dimension,
double* sistem_parametr,
double* xstart,
double* tmp_memory
);


double fun_alpha_nevaska(
    double alpha1,
    double alpha2,
    double* x_alpha,
    double ( *(*F) )(double,const double*),
    double* sistem_parametr,
    double* tmp_memory
);

int derivative_matrix_up(
    const double* x_plus_delta1,
    const double* x_minus_delta1,
    const double* x_plus_delta2,
    const double* x_minus_delta2,
    double delta,
    double* matrix
);

void print_matrix(const double* matrix);

double resh_SLAU(double* h,const double* matrix,const double* b);

void cappa_vector(double* cappa_v,const double* matrix);

double S_function(const double* x,const double* cappa_v);

double vector_alpha_sampling(
    double* new_vector_alpha,
    const double* vector_alpha,
    const double* h,
    const double* cappa_v,
    double ( *(*F) )(double,const double*),
    double* sistem_parametr,
    double* tmp_memory
);

void print_vector(double* vector);

int RKDOPR_8_7_interval_with_control_main(
    double px0,
    double py0
);
