#include<cstdio>
#include<iostream>
#include<cmath>

#include "fun.h"

double tol_local_err_RK=pow(10.,-12);
double tol_intersection=pow(10.,-12);
int p=8;
double fac=0.9;
double facmax=1.3;
double facmin=0.7;
double epsilon=pow(10.,-10);
int N_MAX_vector_alpha=20;
double step_MAX=1.;
double T=10.;



static double   c2=1./18.,
                c3=1./12.,
                c4=0.125,
                c5=0.3125,
                c6=0.375,
                c7=0.1475,
                c8=0.465,
                c9=5490023248./9719169821.,
                c10=0.65,
                c11=1201146811./1299019798.,
                a21=1./18.,
                a31=1./48.,
                a32=0.0625, 
                a41=0.03125,
                a43=0.09375,
                a51=0.3125,
                a53=-1.171875,
                a54=1.171875,
                a61=0.0375,
                a64=0.1875,
                a65=0.15,
                a71=29443841./614563906.,
                a74=77736538./692538347.,
                a75=-28693883./1125000000.,
                a76=23124283./1800000000.,
                a81=16016141./946692911.,
                a84=61564180./158732637.,
                a85=22789713./633445777.,
                a86=545815736./2771057229.,
                a87=-180193667./1043307555.,
                a91=39632708./573591083.,
                a94=-433636366./683701615.,
                a95=-421739975./2616292301.,
                a96=100302831./723423059.,
                a97=790204164./839813087.,
                a98=800635310./3783071287.,
                a101=246121993./1340847787.,
                a104=-37695042795./15268766246.,
                a105=-309121744./1061227803.,
                a106=-12992083./490766935.,
                a107=6005943493./2108947869.,
                a108=393006217./1396673457.,
                a109=123872331./1001029789.,
                a111=-1028468189./846180014.,
                a114=8478235783./508512852.,
                a115=1311729495./1432422823.,
                a116=-10304129995./1701304382.,
                a117=-48777925059./3047939560.,
                a118=15336726248./1032824649.,
                a119=-45442868181./3398467696.,
                a1110=3065993473./597172653.,
                a121=185892177./718116043.,
                a124=-3185094517./667107341.,
                a125=-477755414./1098053517.,
                a126=-703635378./230739211.,
                a127=5731566787./1027545527.,
                a128=5232866602./850066563.,
                a129=-4093664535./808688257.,
                a1210=3962137247./1805957418.,
                a1211=65686358./487910083.,
                a131=403863854./491063109.,
                a134=-5068492393./434740067.,
                a135=-411421997./543043805.,
                a136=652783627./914296604.,
                a137=11173962825./925320556.,
                a138=-13158990841./6184727034.,
                a139=3936647629./1978049680.,
                a1310=-160528059./685178525.,
                a1311=248638103./1413531060.,
                b1=14005451./335480064.,
                b6=-59238493./1068277825.,
                b7=181606767./758867731.,
                b8=561292985./797845732.,
                b9=-1041891430./1371343529.,
                b10=760417239./1151165299.,
                b11=118820643./751138087.,
                b12=-528747749./2220607170.,
                b13=0.25,
                bh1=13451932./455176623.,
                bh6=-808719846./976000145.,
                bh7=1757004468./5645159321.,
                bh8=656045339./265891186.,
                bh9=-3867574721./1518517206.,
                bh10=465885868./322736535.,
                bh11=53011238./667516719.,
                bh12=2./45.;


double fabs(double x)
{
    if(x>0)return x;
    else return -x;
}

void vector1_equality_vector2(int dimension,double* vector1,const double* vector2){
    for(int i=0;i<dimension;i++){
        vector1[i]=vector2[i];
    }
}

double Max(double V1,double V2)
{
    if(V1>V2)return V1;
    else return V2;
}

double Min(double V1,double V2)
{
    if(V1>V2)return V2;
    else return V1;
}

double RKDOPR_8_7(
    int dimension,
    double t,
    const double* x_start,
    double step,
    double* x_end,
    double ( *(*F) )(double ,const double*),
    double* tmp_memory
){
    double* coefficient=tmp_memory;
    double* v_arg=coefficient+13*dimension;
    double* x_end_h=v_arg+dimension;

    double err;

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i];
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+0]=(*(F[i]))(t,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a21)*coefficient[13*i+0]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+1]=(*(F[i]))(t+(c2)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a31)*coefficient[13*i+0]+(a32)*coefficient[13*i+1]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+2]=(*(F[i]))(t+(c3)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a41)*coefficient[13*i+0]+(a43)*coefficient[13*i+2]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+3]=(*(F[i]))(t+(c4)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a51)*coefficient[13*i+0]+(a53)*coefficient[13*i+2]+(a54)*coefficient[13*i+3]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+4]=(*(F[i]))(t+(c5)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a61)*coefficient[13*i+0]+(a64)*coefficient[13*i+3]+(a65)*coefficient[13*i+4]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+5]=(*(F[i]))(t+(c6)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a71)*coefficient[13*i+0]+(a74)*coefficient[13*i+3]+(a75)*coefficient[13*i+4]+(a76)*coefficient[13*i+5]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+6]=(*(F[i]))(t+(c7)*step,v_arg);
    }


    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a81)*coefficient[13*i+0]+(a84)*coefficient[13*i+3]+(a85)*coefficient[13*i+4]+(a86)*coefficient[13*i+5]+(a87)*coefficient[13*i+6]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+7]=(*(F[i]))(t+(c8)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a91)*coefficient[13*i+0]+(a94)*coefficient[13*i+3]+(a95)*coefficient[13*i+4]+(a96)*coefficient[13*i+5]+(a97)*coefficient[13*i+6]+(a98)*coefficient[13*i+7]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+8]=(*(F[i]))(t+(c9)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a101)*coefficient[13*i+0]+(a104)*coefficient[13*i+3]+(a105)*coefficient[13*i+4]+(a106)*coefficient[13*i+5]+(a107)*coefficient[13*i+6]+(a108)*coefficient[13*i+7]+(a109)*coefficient[13*i+8]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+9]=(*(F[i]))(t+(c10)*step,v_arg);
    }


    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a111)*coefficient[13*i+0]+(a114)*coefficient[13*i+3]+(a115)*coefficient[13*i+4]+(a116)*coefficient[13*i+5]+(a117)*coefficient[13*i+6]+(a118)*coefficient[13*i+7]+(a119)*coefficient[13*i+8]+(a1110)*coefficient[13*i+9]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+10]=(*(F[i]))(t+(c11)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a121)*coefficient[13*i+0]+(a124)*coefficient[13*i+3]+(a125)*coefficient[13*i+4]+(a126)*coefficient[13*i+5]+(a127)*coefficient[13*i+6]+(a128)*coefficient[13*i+7]+(a129)*coefficient[13*i+8]+(a1210)*coefficient[13*i+9]+(a1211)*coefficient[13*i+10]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+11]=(*(F[i]))(t+(1.)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a131)*coefficient[13*i+0]+(a134)*coefficient[13*i+3]+(a135)*coefficient[13*i+4]+(a136)*coefficient[13*i+5]+(a137)*coefficient[13*i+6]+(a138)*coefficient[13*i+7]+(a139)*coefficient[13*i+8]+(a1310)*coefficient[13*i+9]+(a1311)*coefficient[13*i+10]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+12]=(*(F[i]))(t+(1.)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)
    {
    x_end[i]=x_start[i]+step*((bh1)*coefficient[13*i+0]+(bh6)*coefficient[13*i+5]+(bh7)*coefficient[13*i+6]+(bh8)*coefficient[13*i+7]+(bh9)*coefficient[13*i+8]+(bh10)*coefficient[13*i+9]+(bh11)*coefficient[13*i+10]+(bh12)*coefficient[13*i+11]);

    x_end_h[i]=x_start[i]+step*((b1)*coefficient[13*i+0]+(b6)*coefficient[13*i+5]+(b7)*coefficient[13*i+6]+(b8)*coefficient[13*i+7]+(b9)*coefficient[13*i+8]+(b10)*coefficient[13*i+9]+(b11)*coefficient[13*i+10]+(b12)*coefficient[13*i+11]+(b13)*coefficient[13*i+12]);
    }

    err=0.;
    for(int i=0;i<dimension;i++)
    if(fabs(x_end[i]-x_end_h[i])>err)
    err=fabs(x_end[i]-x_end_h[i]);

return err;
}

double RKDOPR_8_7_B(
    int dimension,
    double t,
    const double* x_start,
    double step,
    double* x_end,
    double B0,
    double* B,
    double ( *(*F) )(double ,const double*),
    double* tmp_memory
){
    double* coefficient=tmp_memory;
    double* v_arg=coefficient+13*dimension;
    double* x_end_h=v_arg+dimension;
    double* coefficientB=x_end+dimension;
    double B_h;
    double err;

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i];
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+0]=(*(F[i]))(t,v_arg);
    }coefficientB[0]=(*(F[dimension]))(t,v_arg);

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a21)*coefficient[13*i+0]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+1]=(*(F[i]))(t+(c2)*step,v_arg);
    }coefficientB[1]=(*(F[dimension]))(t+(c2)*step,v_arg);

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a31)*coefficient[13*i+0]+(a32)*coefficient[13*i+1]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+2]=(*(F[i]))(t+(c3)*step,v_arg);
    }coefficientB[2]=(*(F[dimension]))(t+(c3)*step,v_arg);

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a41)*coefficient[13*i+0]+(a43)*coefficient[13*i+2]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+3]=(*(F[i]))(t+(c4)*step,v_arg);
    }coefficientB[3]=(*(F[dimension]))(t+(c4)*step,v_arg);

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a51)*coefficient[13*i+0]+(a53)*coefficient[13*i+2]+(a54)*coefficient[13*i+3]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+4]=(*(F[i]))(t+(c5)*step,v_arg);
    }coefficientB[4]=(*(F[dimension]))(t+(c5)*step,v_arg);

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a61)*coefficient[13*i+0]+(a64)*coefficient[13*i+3]+(a65)*coefficient[13*i+4]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+5]=(*(F[i]))(t+(c6)*step,v_arg);
    }coefficientB[5]=(*(F[dimension]))(t+(c6)*step,v_arg);

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a71)*coefficient[13*i+0]+(a74)*coefficient[13*i+3]+(a75)*coefficient[13*i+4]+(a76)*coefficient[13*i+5]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+6]=(*(F[i]))(t+(c7)*step,v_arg);
    }coefficientB[6]=(*(F[dimension]))(t+(c7)*step,v_arg);


    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a81)*coefficient[13*i+0]+(a84)*coefficient[13*i+3]+(a85)*coefficient[13*i+4]+(a86)*coefficient[13*i+5]+(a87)*coefficient[13*i+6]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+7]=(*(F[i]))(t+(c8)*step,v_arg);
    }coefficientB[7]=(*(F[dimension]))(t+(c8)*step,v_arg);

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a91)*coefficient[13*i+0]+(a94)*coefficient[13*i+3]+(a95)*coefficient[13*i+4]+(a96)*coefficient[13*i+5]+(a97)*coefficient[13*i+6]+(a98)*coefficient[13*i+7]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+8]=(*(F[i]))(t+(c9)*step,v_arg);
    }coefficientB[8]=(*(F[dimension]))(t+(c9)*step,v_arg);

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a101)*coefficient[13*i+0]+(a104)*coefficient[13*i+3]+(a105)*coefficient[13*i+4]+(a106)*coefficient[13*i+5]+(a107)*coefficient[13*i+6]+(a108)*coefficient[13*i+7]+(a109)*coefficient[13*i+8]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+9]=(*(F[i]))(t+(c10)*step,v_arg);
    }coefficientB[9]=(*(F[dimension]))(t+(c10)*step,v_arg);


    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a111)*coefficient[13*i+0]+(a114)*coefficient[13*i+3]+(a115)*coefficient[13*i+4]+(a116)*coefficient[13*i+5]+(a117)*coefficient[13*i+6]+(a118)*coefficient[13*i+7]+(a119)*coefficient[13*i+8]+(a1110)*coefficient[13*i+9]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+10]=(*(F[i]))(t+(c11)*step,v_arg);
    }coefficientB[10]=(*(F[dimension]))(t+(c11)*step,v_arg);

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a121)*coefficient[13*i+0]+(a124)*coefficient[13*i+3]+(a125)*coefficient[13*i+4]+(a126)*coefficient[13*i+5]+(a127)*coefficient[13*i+6]+(a128)*coefficient[13*i+7]+(a129)*coefficient[13*i+8]+(a1210)*coefficient[13*i+9]+(a1211)*coefficient[13*i+10]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+11]=(*(F[i]))(t+(1.)*step,v_arg);
    }coefficientB[11]=(*(F[dimension]))(t+(1.)*step,v_arg);

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a131)*coefficient[13*i+0]+(a134)*coefficient[13*i+3]+(a135)*coefficient[13*i+4]+(a136)*coefficient[13*i+5]+(a137)*coefficient[13*i+6]+(a138)*coefficient[13*i+7]+(a139)*coefficient[13*i+8]+(a1310)*coefficient[13*i+9]+(a1311)*coefficient[13*i+10]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+12]=(*(F[i]))(t+(1.)*step,v_arg);
    }coefficientB[12]=(*(F[dimension]))(t+(1.)*step,v_arg);

    for(int i=0;i<dimension;i++)
    {
    x_end[i]=x_start[i]+step*((bh1)*coefficient[13*i+0]+(bh6)*coefficient[13*i+5]+(bh7)*coefficient[13*i+6]+(bh8)*coefficient[13*i+7]+(bh9)*coefficient[13*i+8]+(bh10)*coefficient[13*i+9]+(bh11)*coefficient[13*i+10]+(bh12)*coefficient[13*i+11]);

    x_end_h[i]=x_start[i]+step*((b1)*coefficient[13*i+0]+(b6)*coefficient[13*i+5]+(b7)*coefficient[13*i+6]+(b8)*coefficient[13*i+7]+(b9)*coefficient[13*i+8]+(b10)*coefficient[13*i+9]+(b11)*coefficient[13*i+10]+(b12)*coefficient[13*i+11]+(b13)*coefficient[13*i+12]);
    }
    *B=B0+step*((bh1)*coefficientB[0]+(bh6)*coefficientB[5]+(bh7)*coefficientB[6]+(bh8)*coefficientB[7]+(bh9)*coefficientB[8]+(bh10)*coefficientB[9]+(bh11)*coefficientB[10]+(bh12)*coefficientB[11]);

    err=0.;
    for(int i=0;i<dimension;i++){
        if(fabs(x_end[i]-x_end_h[i])>err)
        err=fabs(x_end[i]-x_end_h[i]);
    }


return err;
}

double RKDOPR_8_7_with_vector_errors(
    int dimension,
    double t,
    const double* x_start,
    double step,
    double* x_end,
    double ( *(*F) )(double ,const double*),
    double* vector_errors,//len(vector_errors)=dimension
    double* tmp_memory
){
    //lens(tmp_memory)>13*dimension + dimension + dimension
    double* coefficient=tmp_memory;
    double* v_arg=coefficient+13*dimension;
    double* x_end_h=v_arg+dimension;

    double err;

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i];
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+0]=(*(F[i]))(t,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a21)*coefficient[13*i+0]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+1]=(*(F[i]))(t+(c2)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a31)*coefficient[13*i+0]+(a32)*coefficient[13*i+1]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+2]=(*(F[i]))(t+(c3)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a41)*coefficient[13*i+0]+(a43)*coefficient[13*i+2]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+3]=(*(F[i]))(t+(c4)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a51)*coefficient[13*i+0]+(a53)*coefficient[13*i+2]+(a54)*coefficient[13*i+3]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+4]=(*(F[i]))(t+(c5)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a61)*coefficient[13*i+0]+(a64)*coefficient[13*i+3]+(a65)*coefficient[13*i+4]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+5]=(*(F[i]))(t+(c6)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a71)*coefficient[13*i+0]+(a74)*coefficient[13*i+3]+(a75)*coefficient[13*i+4]+(a76)*coefficient[13*i+5]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+6]=(*(F[i]))(t+(c7)*step,v_arg);
    }


    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a81)*coefficient[13*i+0]+(a84)*coefficient[13*i+3]+(a85)*coefficient[13*i+4]+(a86)*coefficient[13*i+5]+(a87)*coefficient[13*i+6]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+7]=(*(F[i]))(t+(c8)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a91)*coefficient[13*i+0]+(a94)*coefficient[13*i+3]+(a95)*coefficient[13*i+4]+(a96)*coefficient[13*i+5]+(a97)*coefficient[13*i+6]+(a98)*coefficient[13*i+7]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+8]=(*(F[i]))(t+(c9)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a101)*coefficient[13*i+0]+(a104)*coefficient[13*i+3]+(a105)*coefficient[13*i+4]+(a106)*coefficient[13*i+5]+(a107)*coefficient[13*i+6]+(a108)*coefficient[13*i+7]+(a109)*coefficient[13*i+8]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+9]=(*(F[i]))(t+(c10)*step,v_arg);
    }


    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a111)*coefficient[13*i+0]+(a114)*coefficient[13*i+3]+(a115)*coefficient[13*i+4]+(a116)*coefficient[13*i+5]+(a117)*coefficient[13*i+6]+(a118)*coefficient[13*i+7]+(a119)*coefficient[13*i+8]+(a1110)*coefficient[13*i+9]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+10]=(*(F[i]))(t+(c11)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a121)*coefficient[13*i+0]+(a124)*coefficient[13*i+3]+(a125)*coefficient[13*i+4]+(a126)*coefficient[13*i+5]+(a127)*coefficient[13*i+6]+(a128)*coefficient[13*i+7]+(a129)*coefficient[13*i+8]+(a1210)*coefficient[13*i+9]+(a1211)*coefficient[13*i+10]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+11]=(*(F[i]))(t+(1.)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((a131)*coefficient[13*i+0]+(a134)*coefficient[13*i+3]+(a135)*coefficient[13*i+4]+(a136)*coefficient[13*i+5]+(a137)*coefficient[13*i+6]+(a138)*coefficient[13*i+7]+(a139)*coefficient[13*i+8]+(a1310)*coefficient[13*i+9]+(a1311)*coefficient[13*i+10]);
    for(int i=0;i<dimension;i++){ 
        coefficient[13*i+12]=(*(F[i]))(t+(1.)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)
    {
    x_end[i]=x_start[i]+step*((bh1)*coefficient[13*i+0]+(bh6)*coefficient[13*i+5]+(bh7)*coefficient[13*i+6]+(bh8)*coefficient[13*i+7]+(bh9)*coefficient[13*i+8]+(bh10)*coefficient[13*i+9]+(bh11)*coefficient[13*i+10]+(bh12)*coefficient[13*i+11]);

    x_end_h[i]=x_start[i]+step*((b1)*coefficient[13*i+0]+(b6)*coefficient[13*i+5]+(b7)*coefficient[13*i+6]+(b8)*coefficient[13*i+7]+(b9)*coefficient[13*i+8]+(b10)*coefficient[13*i+9]+(b11)*coefficient[13*i+10]+(b12)*coefficient[13*i+11]+(b13)*coefficient[13*i+12]);
    }

    err=0.;
    for(int i=0;i<dimension;i++){
        if(fabs(x_end[i]-x_end_h[i])>err){
            err=fabs(x_end[i]-x_end_h[i]);
        }
    vector_errors[i]=fabs(x_end[i]-x_end_h[i]);
    }

return err;
}



static double   ic2=1./5.,
                ic3=3./10.,
                ic4=4./5.,
                ic5=8./9.,
                ic6=1.,
                ic7=1.,
                ia21=1./5.,
                ia31=3./40.,
                ia32=9./40., 
                ia41=44./45.,
                ia42=-56./15.,
                ia43=32./9.,
                ia51=19372./6561.,
                ia52=-25360./2187.,
                ia53=64448./6561.,
                ia54=-212./729.,
                ia61=9017./3168.,
                ia62=-355./33.,
                ia63=46732./5247.,
                ia64=49./176.,
                ia65=-5103./18656.,
                ia71=35./384.,
                ia72=0.,
                ia73=500./1113.,
                ia74=125./192.,
                ia75=-2187./6784.,
                ia76=11./84.,
                ib1=35./384.,
                ib2=0.,
                ib3=500./1113.,
                ib4=125./192.,
                ib5=-2187./6784.,
                ib6=11./84.,
                ib7=0.,
                ibh1=5179./57600.,
                ibh2=0.,
                ibh3=7571./16695.,
                ibh4=393./640.,
                ibh5=-92097./339200.,
                ibh6=187./2100.,
                ibh7=1./40.;

double RKDOPR_5_4(
    int dimension,
    double t,
    const double* x_start,
    double step,
    double* x_end,
    double ( *(*F) )(double ,const double*),
    double* tmp_memory
){
    double* coefficient=tmp_memory;
    double* v_arg=coefficient+7*dimension;
    double* x_end_h=v_arg+dimension;

    double err;

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i];
    for(int i=0;i<dimension;i++){ 
        coefficient[7*i+0]=(*(F[i]))(t,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((ia21)*coefficient[7*i+0]);
    for(int i=0;i<dimension;i++){ 
        coefficient[7*i+1]=(*(F[i]))(t+(ic2)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((ia31)*coefficient[7*i+0]+(ia32)*coefficient[7*i+1]);
    for(int i=0;i<dimension;i++){ 
        coefficient[7*i+2]=(*(F[i]))(t+(ic3)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((ia41)*coefficient[7*i+0]+(ia42)*coefficient[7*i+1]+(ia43)*coefficient[7*i+2]);
    for(int i=0;i<dimension;i++){ 
        coefficient[7*i+3]=(*(F[i]))(t+(ic4)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((ia51)*coefficient[7*i+0]+(ia52)*coefficient[7*i+1]+(ia53)*coefficient[7*i+2]+(ia54)*coefficient[7*i+3]);
    for(int i=0;i<dimension;i++){ 
        coefficient[7*i+4]=(*(F[i]))(t+(ic5)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((ia61)*coefficient[7*i+0]+(ia62)*coefficient[7*i+1]+(ia63)*coefficient[7*i+2]+(ia64)*coefficient[7*i+3]+(ia65)*coefficient[7*i+4]);
    for(int i=0;i<dimension;i++){ 
        coefficient[7*i+5]=(*(F[i]))(t+(ic6)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)v_arg[i]=x_start[i]+step*((ia71)*coefficient[7*i+0]+(ia72)*coefficient[7*i+1]+(ia73)*coefficient[7*i+2]+(ia74)*coefficient[7*i+3]+(ia75)*coefficient[7*i+4]+(ia76)*coefficient[7*i+5]);
    for(int i=0;i<dimension;i++){ 
        coefficient[7*i+6]=(*(F[i]))(t+(ic7)*step,v_arg);
    }

    for(int i=0;i<dimension;i++)
    {
    x_end_h[i]=x_start[i]+step*((ib1)*coefficient[7*i+0]+(ib2)*coefficient[7*i+1]+(ib3)*coefficient[7*i+2]+(ib4)*coefficient[7*i+3]+(ib5)*coefficient[7*i+4]+(ib6)*coefficient[7*i+5]+(ib7)*coefficient[7*i+6]);

    x_end[i]=x_start[i]+step*((ibh1)*coefficient[7*i+0]+(ibh2)*coefficient[7*i+1]+(ibh3)*coefficient[7*i+2]+(ibh4)*coefficient[7*i+3]+(ibh5)*coefficient[7*i+4]+(ibh6)*coefficient[7*i+5]+(ibh7)*coefficient[7*i+6]);
    }

    err=0.;
    for(int i=0;i<dimension;i++)
    if(fabs(x_end[i]-x_end_h[i])>err)
    err=fabs(x_end[i]-x_end_h[i]);

return err;
}


// double* tback=sistem_parametr;
// double* xback=sistem_parametr+1;
// double* tnow=sistem_parametr+dimension+1;
// double* xnow=sistem_parametr+2+dimension;
// double* tnext=sistem_parametr+2+2*dimension;
// double* xnext=sistem_parametr+3+2*dimension;
// double T=sistem_parametr[3+3*dimension];
// double step=sistem_parametr[3+3*dimension+1];
// double* vector_err=sistem_parametr+5+3*dimension;
// double global_err=sistem_parametr[5+4*dimension];





double IntersectionSearch_with_increases_up(
int dimension,
double* sistem_parametr,
double ( *(*F) )(double ,const double*),
double* tmp_memory
){

double* tnow=sistem_parametr+dimension+1;
double* xnow=sistem_parametr+2+dimension;
double* tnext=sistem_parametr+2+2*dimension;
double* xnext=sistem_parametr+3+2*dimension;

double* xnow_tmp=tmp_memory;
vector1_equality_vector2(dimension,xnow_tmp,xnow);
double tnow_tmp=*tnow;

double* x_tmp=xnow_tmp+dimension;
double* tmp_memory1=x_tmp+dimension;
double step;

do{
    step=(*tnext-tnow_tmp)/2.;
    RKDOPR_8_7(dimension,tnow_tmp,xnow_tmp,step,x_tmp,F,tmp_memory1);

        if( fabs(   x_tmp[dimension-1]-1.  )<=tol_intersection){
            vector1_equality_vector2(dimension,xnext,x_tmp);
            *tnext=tnow_tmp+step;
            return 0.;
        }

        if(x_tmp[dimension-1]<1.){
            vector1_equality_vector2(dimension,xnow_tmp,x_tmp);
            tnow_tmp=tnow_tmp+step;
        }else if(x_tmp[dimension-1]>1.){
            vector1_equality_vector2(dimension,xnext,x_tmp);
            *tnext=tnow_tmp+step;
        }

    }while(1);

return -1.;
}

double IntersectionSearch_with_increases_bot(
int dimension,
double* sistem_parametr,
double ( *(*F) )(double ,const double*),
double* tmp_memory
){

double* tnow=sistem_parametr+dimension+1;
double* xnow=sistem_parametr+2+dimension;
double* tnext=sistem_parametr+2+2*dimension;
double* xnext=sistem_parametr+3+2*dimension;

double* xnow_tmp=tmp_memory;
vector1_equality_vector2(dimension,xnow_tmp,xnow);
double tnow_tmp=*tnow;

double* x_tmp=xnow_tmp+dimension;
double* tmp_memory1=x_tmp+dimension;
double step;

do{
    step=(*tnext-tnow_tmp)/2.;
    RKDOPR_8_7(dimension,tnow_tmp,xnow_tmp,step,x_tmp,F,tmp_memory1);

        if( fabs(   x_tmp[dimension-1]+1.  )<=tol_intersection){
            vector1_equality_vector2(dimension,xnext,x_tmp);
            *tnext=tnow_tmp+step;
            return 0.;
        }

        if(x_tmp[dimension-1]<-1.){
            vector1_equality_vector2(dimension,xnow_tmp,x_tmp);
            tnow_tmp=tnow_tmp+step;
        }else if(x_tmp[dimension-1]>-1.){
            vector1_equality_vector2(dimension,xnext,x_tmp);
            *tnext=tnow_tmp+step;
        }

    }while(1);

return -1.;
}


double IntersectionSearch_with_decreases_up(
int dimension,
double* sistem_parametr,
double ( *(*F) )(double ,const double*),
double* tmp_memory
){

double* tnow=sistem_parametr+dimension+1;
double* xnow=sistem_parametr+2+dimension;
double* tnext=sistem_parametr+2+2*dimension;
double* xnext=sistem_parametr+3+2*dimension;

double* xnow_tmp=tmp_memory;
vector1_equality_vector2(dimension,xnow_tmp,xnow);
double tnow_tmp=*tnow;

double* x_tmp=xnow_tmp+dimension;
double* tmp_memory1=x_tmp+dimension;
double step;

do{
    step=(*tnext-tnow_tmp)/2.;
    RKDOPR_8_7(dimension,tnow_tmp,xnow_tmp,step,x_tmp,F,tmp_memory1);

        if( fabs(   x_tmp[dimension-1]-1.  )<=tol_intersection){
            vector1_equality_vector2(dimension,xnext,x_tmp);
            *tnext=tnow_tmp+step;
            return 0.;
        }

        if(x_tmp[dimension-1]<1.){
            vector1_equality_vector2(dimension,xnext,x_tmp);
            *tnext=tnow_tmp+step; 
        }else if(x_tmp[dimension-1]>1.){
            vector1_equality_vector2(dimension,xnow_tmp,x_tmp);
            tnow_tmp=tnow_tmp+step;
        }

    }while(1);

return -1.;
}

double IntersectionSearch_with_decreases_bot(
int dimension,
double* sistem_parametr,
double ( *(*F) )(double ,const double*),
double* tmp_memory
){
double* tnow=sistem_parametr+dimension+1;
double* xnow=sistem_parametr+2+dimension;
double* tnext=sistem_parametr+2+2*dimension;
double* xnext=sistem_parametr+3+2*dimension;

double* xnow_tmp=tmp_memory;
vector1_equality_vector2(dimension,xnow_tmp,xnow);
double tnow_tmp=*tnow;

double* x_tmp=xnow_tmp+dimension;
double* tmp_memory1=x_tmp+dimension;
double step;

do{
    step=(*tnext-tnow_tmp)/2.;
    RKDOPR_8_7(dimension,tnow_tmp,xnow_tmp,step,x_tmp,F,tmp_memory1);

        if( fabs(   x_tmp[dimension-1]+1.  )<=tol_intersection){
            vector1_equality_vector2(dimension,xnext,x_tmp);
            *tnext=tnow_tmp+step;
            return 0.;
        }

        if(x_tmp[dimension-1]<-1.){
            vector1_equality_vector2(dimension,xnext,x_tmp);
            *tnext=tnow_tmp+step; 
        }else if(x_tmp[dimension-1]>-1.){
            vector1_equality_vector2(dimension,xnow_tmp,x_tmp);
            tnow_tmp=tnow_tmp+step;
        }

    }while(1);

return -1.;
}


void shift(
int dimension,
double* sistem_parametr
){
double* tback=sistem_parametr;
double* xback=sistem_parametr+1;
double* tnow=sistem_parametr+dimension+1;
double* xnow=sistem_parametr+2+dimension;
double* tnext=sistem_parametr+2+2*dimension;
double* xnext=sistem_parametr+3+2*dimension;

    *tback=*tnow;
    *tnow=*tnext;
    for(int i=0;i<dimension;i++){
        xback[i]=xnow[i];
        xnow[i]=xnext[i];
    }
}

FILE* file1=fopen("t1.txt","w");
FILE* file2=fopen("t2.txt","w");
FILE* file3=fopen("t3.txt","w");
FILE* file4=fopen("t4.txt","w");
FILE* file5=fopen("t5.txt","w");


void shift_print(
int dimension,
double* sistem_parametr
){
double* tback=sistem_parametr;
double* xback=sistem_parametr+1;
double* tnow=sistem_parametr+dimension+1;
double* xnow=sistem_parametr+2+dimension;
double* tnext=sistem_parametr+2+2*dimension;
double* xnext=sistem_parametr+3+2*dimension;

double u;
if(xnow[3]>1.)u=1.;
else if(xnow[3]<-1.)u=-1.;
else u=xnow[3];

fprintf(file1,"%lf   ",*tnow);for(int i=0;i<dimension;i++){fprintf(file1,"%le   ",xnow[i]);}fprintf(file1,"%le\n",u);

    *tback=*tnow;
    *tnow=*tnext;
    for(int i=0;i<dimension;i++){
        xback[i]=xnow[i];
        xnow[i]=xnext[i];
    }

return;
}

double F_x_main(double t,const double* x){
    return x[1];
}

double F_y_main(double t, const double* x){
    if(x[3]>1.){return 1.;}
    else if(x[3]<-1.){return -1.;}
    else{ return x[3];}
}

double F_px_main(double t, const double* x){
    return -x[0];
}

double F_py_main(double t, const double* x){
    return -x[2]-x[1];
}

double F_B(double t, const double* x){
    if(x[3]>1.){return 1.-x[1]*x[1]-x[0]*x[0];}
    else if(x[3]<-1.){return 1.-x[1]*x[1]-x[0]*x[0];}
    else{ return x[3]*x[3]-x[1]*x[1]-x[0]*x[0];}
}

double H(double t, const double* x){
    if(x[3]>1.){return x[2]*x[1]+x[3]-1/2.*(1.-x[1]*x[1]-x[0]*x[0]);}
    else if(x[3]<-1.){return x[2]*x[1]-x[3]-1/2.*(1.-x[1]*x[1]-x[0]*x[0]);}
    else{return x[2]*x[1]+x[3]*x[3]-1/2.*(x[3]*x[3]-x[1]*x[1]-x[0]*x[0]);}
}

double RKDOPR_8_7_interval_with_control_printf(
const int dimension,
double* sistem_parametr,
double* xstart,
double* tmp_memory
){

double* tback=sistem_parametr;
double* xback=sistem_parametr+1;
double* tnow=sistem_parametr+4+1;
double* xnow=sistem_parametr+2+4;
double* tnext=sistem_parametr+2+2*4;
double* xnext=sistem_parametr+3+2*4;

vector1_equality_vector2(4,xback,xstart);
vector1_equality_vector2(4,xnow,xstart);
vector1_equality_vector2(4,xnext,xstart);
*tback=0.;
*tnow=0.;
*tnext=0.;

double ( *(F[4]) )(double const,double const*);
F[0]=F_x_main;
F[1]=F_y_main;
F[2]=F_px_main;
F[3]=F_py_main;

double step=0.1;
double global_err=0.;
double local_err;

if(
fabs(xnow[3])<pow(10.,-9) &&
fabs(xnow[2])<pow(10.,-9)
){
return -1;
}

while(*tnow<T){
    if(*tnow+step>T)step=T-*tnow;

    local_err=RKDOPR_8_7(4,*tnow,xnow,step,xnext,F,tmp_memory);

    if( local_err<tol_local_err_RK){
        *tnext=*tnow+step;

        if( (xnow[3]>1.+tol_intersection)&&
            (xnext[3]<1.-tol_intersection)   ){
                IntersectionSearch_with_decreases_up(4,sistem_parametr,F,tmp_memory);
        }else if(fabs(xnow[3]-1.)<=tol_intersection && xnext[3]<0.){
                step=step/2.;continue;
        }else if(fabs(xnow[3])<1.-tol_intersection){
            if(xnext[3]>1.+tol_intersection)
            {
                IntersectionSearch_with_increases_up(4,sistem_parametr,F,tmp_memory);
            }else if(xnext[3]<-1.-tol_intersection){
                IntersectionSearch_with_decreases_bot(4,sistem_parametr,F,tmp_memory);
            }
        }else if( (fabs(xnow[3]+1.)<=tol_intersection)&&
                  (xnext[3]>0.)    ){
                step=step/2.;continue;
        }else if( (xnow[3]<-1.-tol_intersection)&&
                  (xnext[3]>-1.+tol_intersection)   ){
                IntersectionSearch_with_increases_bot(4,sistem_parametr,F,tmp_memory);
        }

    fprintf(file1,"%lf   ",*tnow);for(int i=0;i<dimension;i++){fprintf(file1,"%le   ",xnow[i]);}fprintf(file1,"%le",(*(F[1]))(*tnow,xnow));
    if((fabs(xnow[3]+1.)<=tol_intersection) || fabs(xnow[3]-1.)<=tol_intersection)fprintf(file1,"  1\n");
    else fprintf(file1,"  0\n");

    *tnow=*tnext;
    for(int i=0;i<4;i++){xnow[i]=xnext[i];}

    }
step=step*Min( facmax , Max(  facmin  ,  fac*pow(tol_local_err_RK/local_err,1./(p+1))  ) );
}
*tnow=*tnext;
for(int i=0;i<4;i++){xnow[i]=xnext[i];}
fprintf(file1,"%lf   ",*tnow);for(int i=0;i<dimension;i++){fprintf(file1,"%le   ",xnow[i]);}fprintf(file1,"%le",(*(F[1]))(*tnow,xnow));

if((fabs(xnow[3]+1.)<=tol_intersection) || fabs(xnow[3]-1.)<=tol_intersection)fprintf(file1,"  1\n");
else fprintf(file1,"  0\n");

fprintf(file1,"END\n");

return 1.;
}

double RKDOPR_8_7_interval_with_control(
const int dimension,
double* sistem_parametr,
double* xstart,
double* tmp_memory
){

double* tback=sistem_parametr;
double* xback=sistem_parametr+1;
double* tnow=sistem_parametr+4+1;
double* xnow=sistem_parametr+2+4;
double* tnext=sistem_parametr+2+2*4;
double* xnext=sistem_parametr+3+2*4;

vector1_equality_vector2(4,xback,xstart);
vector1_equality_vector2(4,xnow,xstart);
vector1_equality_vector2(4,xnext,xstart);
*tback=0.;
*tnow=0.;
*tnext=0.;

double ( *(F[4]) )(double const,double const*);
F[0]=F_x_main;
F[1]=F_y_main;
F[2]=F_px_main;
F[3]=F_py_main;

double step=0.1;
double global_err=0.;
double local_err;

if(
fabs(xnow[3])<pow(10.,-9) &&
fabs(xnow[2])<pow(10.,-9)
){
return -1;
}

while(*tnow<T){
    if(*tnow+step>T)step=T-*tnow;

    local_err=RKDOPR_8_7(4,*tnow,xnow,step,xnext,F,tmp_memory);

    if( local_err<tol_local_err_RK){
        *tnext=*tnow+step;

        if( (xnow[3]>1.+tol_intersection)&&
            (xnext[3]<1.-tol_intersection)   ){
                IntersectionSearch_with_decreases_up(4,sistem_parametr,F,tmp_memory);
        }else if(fabs(xnow[3]-1.)<=tol_intersection && xnext[3]<0.){
                step=step/2.;continue;
        }else if(fabs(xnow[3])<1.-tol_intersection){
            if(xnext[3]>1.+tol_intersection)
            {
                IntersectionSearch_with_increases_up(4,sistem_parametr,F,tmp_memory);
            }else if(xnext[3]<-1.-tol_intersection){
                IntersectionSearch_with_decreases_bot(4,sistem_parametr,F,tmp_memory);
            }
        }else if( (fabs(xnow[3]+1.)<=tol_intersection)&&
                  (xnext[3]>0.)    ){
                step=step/2.;continue;
        }else if( (xnow[3]<-1.-tol_intersection)&&
                  (xnext[3]>-1.+tol_intersection)   ){
                IntersectionSearch_with_increases_bot(4,sistem_parametr,F,tmp_memory);
        }

    *tnow=*tnext;
    for(int i=0;i<4;i++){xnow[i]=xnext[i];}

    }
step=step*Min( facmax , Max(  facmin  ,  fac*pow(tol_local_err_RK/local_err,1./(p+1))  ) );
}
*tnow=*tnext;
for(int i=0;i<4;i++){xnow[i]=xnext[i];}

return 1.;
}

double RKDOPR_8_7_interval_with_control_printf_B(
const int dimension,
double* sistem_parametr,
double* xstart,
double* tmp_memory
){

double* tback=sistem_parametr;
double* xback=sistem_parametr+1;
double* tnow=sistem_parametr+4+1;
double* xnow=sistem_parametr+2+4;
double* tnext=sistem_parametr+2+2*4;
double* xnext=sistem_parametr+3+2*4;

vector1_equality_vector2(4,xback,xstart);
vector1_equality_vector2(4,xnow,xstart);
vector1_equality_vector2(4,xnext,xstart);
*tback=0.;
*tnow=0.;
*tnext=0.;

double ( *(F[5]) )(double const,double const*);
F[0]=F_x_main;
F[1]=F_y_main;
F[2]=F_px_main;
F[3]=F_py_main;
F[4]=F_B;

double step=0.1;
double global_err=0.;
double local_err;
double B0,B;B0=0.;

if(
fabs(xnow[3])<pow(10.,-9) &&
fabs(xnow[2])<pow(10.,-9)
){
return -1;
}

while(*tnow<T){
    if(step>1)step=1;
    if(*tnow+step>T)step=T-*tnow;

    local_err=RKDOPR_8_7_B(4,*tnow,xnow,step,xnext,B0,&B,F,tmp_memory);

    if( local_err<tol_local_err_RK){
        *tnext=*tnow+step;

        if( (xnow[3]>1.+tol_intersection)&&
            (xnext[3]<1.-tol_intersection)   ){
                IntersectionSearch_with_decreases_up(4,sistem_parametr,F,tmp_memory);
                RKDOPR_8_7_B(4,*tnow,xnow,*tnext-*tnow,xnext,B0,&B,F,tmp_memory);
        }else if(fabs(xnow[3]-1.)<=tol_intersection && xnext[3]<0.){
                step=step/2.;continue;
        }else if(fabs(xnow[3])<1.-tol_intersection){
            if(xnext[3]>1.+tol_intersection)
            {
                IntersectionSearch_with_increases_up(4,sistem_parametr,F,tmp_memory);
                RKDOPR_8_7_B(4,*tnow,xnow,*tnext-*tnow,xnext,B0,&B,F,tmp_memory);
            }else if(xnext[3]<-1.-tol_intersection){
                IntersectionSearch_with_decreases_bot(4,sistem_parametr,F,tmp_memory);
                RKDOPR_8_7_B(4,*tnow,xnow,*tnext-*tnow,xnext,B0,&B,F,tmp_memory);
            }
        }else if( (fabs(xnow[3]+1.)<=tol_intersection)&&
                  (xnext[3]>0.)    ){
                step=step/2.;continue;
        }else if( (xnow[3]<-1.-tol_intersection)&&
                  (xnext[3]>-1.+tol_intersection)   ){
                IntersectionSearch_with_increases_bot(4,sistem_parametr,F,tmp_memory);
                RKDOPR_8_7_B(4,*tnow,xnow,*tnext-*tnow,xnext,B0,&B,F,tmp_memory);
        }

    fprintf(file3,"%lf   ",*tnow);for(int i=0;i<dimension;i++){fprintf(file3,"%le   ",xnow[i]);}fprintf(file3,"%le",(*(F[1]))(*tnow,xnow));
    if((fabs(xnow[3]+1.)<=tol_intersection) || fabs(xnow[3]-1.)<=tol_intersection)fprintf(file3,"  1\n");
    else fprintf(file3,"  0\n");

    *tnow=*tnext;
    B0=B;
    for(int i=0;i<4;i++){xnow[i]=xnext[i];}
    }
step=step*Min( facmax , Max(  facmin  ,  fac*pow(tol_local_err_RK/local_err,1./(p+1))  ) );
}
*tnow=*tnext;
B0=B;
for(int i=0;i<4;i++){xnow[i]=xnext[i];}
fprintf(file3,"%lf   ",*tnow);for(int i=0;i<dimension;i++){fprintf(file3,"%le   ",xnow[i]);}fprintf(file3,"%le",(*(F[1]))(*tnow,xnow));

if((fabs(xnow[3]+1.)<=tol_intersection) || fabs(xnow[3]-1.)<=tol_intersection)fprintf(file3,"  1\n");
else fprintf(file3,"  0\n");
fprintf(file3,"%lf\n",B0);
fprintf(file3,"END\n");

return 1.;
}


double RKDOPR_8_7_interval_with_control_err(
const int dimension,
double* sistem_parametr,
double* xstart,
double* tmp_memory
){

double* tback=sistem_parametr;
double* xback=sistem_parametr+1;
double* tnow=sistem_parametr+4+1;
double* xnow=sistem_parametr+2+4;
double* tnext=sistem_parametr+2+2*4;
double* xnext=sistem_parametr+3+2*4;

vector1_equality_vector2(4,xback,xstart);
vector1_equality_vector2(4,xnow,xstart);
vector1_equality_vector2(4,xnext,xstart);
*tback=0.;
*tnow=0.;
*tnext=0.;

double ( *(F[4]) )(double const,double const*);
F[0]=F_x_main;
F[1]=F_y_main;
F[2]=F_px_main;
F[3]=F_py_main;

double step=1;
double global_err=0.;
double local_err;
double vector_err[4];
double global_vector_err[4];
global_vector_err[0]=0.;
global_vector_err[1]=0.;
global_vector_err[2]=0.;
global_vector_err[3]=0.;

if(
fabs(xnow[3])<pow(10.,-9) &&
fabs(xnow[2])<pow(10.,-9)
){
return -1;
}
double min_h=H(0.,xstart);
double max_h=H(0.,xstart);


while(*tnow<T){
    if(*tnow+step>T)step=T-*tnow;

    local_err=RKDOPR_8_7_with_vector_errors(4,*tnow,xnow,step,xnext,F,vector_err,tmp_memory);

    if( local_err<tol_local_err_RK){
        *tnext=*tnow+step;

        if( (xnow[3]>1.+tol_intersection)&&
            (xnext[3]<1.-tol_intersection)   ){
                IntersectionSearch_with_decreases_up(4,sistem_parametr,F,tmp_memory);
        }else if(fabs(xnow[3]-1.)<=tol_intersection && xnext[3]<0.){
                step=step/2.;continue;
        }else if(fabs(xnow[3])<1.-tol_intersection){
            if(xnext[3]>1.+tol_intersection)
            {
                IntersectionSearch_with_increases_up(4,sistem_parametr,F,tmp_memory);
            }else if(xnext[3]<-1.-tol_intersection){
                IntersectionSearch_with_decreases_bot(4,sistem_parametr,F,tmp_memory);
            }
        }else if( (fabs(xnow[3]+1.)<=tol_intersection)&&
                  (xnext[3]>0.)    ){
                step=step/2.;continue;
        }else if( (xnow[3]<-1.-tol_intersection)&&
                  (xnext[3]>-1.+tol_intersection)   ){
                IntersectionSearch_with_increases_bot(4,sistem_parametr,F,tmp_memory);
        }

    if(fabs(xnow[3])<=1.)global_err=global_err*exp( ((1+sqrt(5))/4)*(*tnext-*tnow) )+local_err;
    else global_err=global_err*exp( ( 1/(sqrt(2)) )*(*tnext-*tnow) )+local_err;
    
    if(H(*tnext,xnext)<min_h)min_h=H(*tnext,xnext);
    if(H(*tnext,xnext)>max_h)max_h=H(*tnext,xnext);

    *tnow=*tnext;
    for(int i=0;i<4;i++){
        xnow[i]=xnext[i];
        global_vector_err[i]=global_vector_err[i]+vector_err[i];
        }
    }
step=step*Min( facmax , Max(  facmin  ,  fac*pow(tol_local_err_RK/local_err,1./(p+1))  ) );
}
*tnow=*tnext;
for(int i=0;i<4;i++){
    xnow[i]=xnext[i];
    global_vector_err[i]=global_vector_err[i]+vector_err[i];
}

fprintf(file4,"$%lf$&$%le$&$%le$&$%le$&$%le$&$%le$&$%le$&$%le$\n",T,tol_local_err_RK,global_vector_err[0],global_vector_err[1],global_vector_err[2],global_vector_err[3],global_err,max_h-min_h);

return 1.;
}


double fun_alpha_nevaska(
    double alpha1,
    double alpha2,
    double* x_alpha,
    double ( *(*F) )(double,const double*),
    double* sistem_parametr,
    double* tmp_memory
){
double* x_start=tmp_memory;
double* tmp_memory1=x_start+4;
x_start[0]=0.;
x_start[1]=0.;
x_start[2]=alpha1;
x_start[3]=alpha2;

if(RKDOPR_8_7_interval_with_control(4,sistem_parametr,x_start,tmp_memory1)<0.){return -1;}

double* xnow=sistem_parametr+2+4;

x_alpha[0]=xnow[0];
x_alpha[1]=xnow[3];

return 1.;
}

double fun_alpha_nevaska_printf(
    double alpha1,
    double alpha2,
    double* x_alpha,
    double ( *(*F) )(double,const double*),
    double* sistem_parametr,
    double* tmp_memory
){
double* x_start=tmp_memory;
double* tmp_memory1=x_start+4;
x_start[0]=0.;
x_start[1]=0.;
x_start[2]=alpha1;
x_start[3]=alpha2;

if(RKDOPR_8_7_interval_with_control_printf(4,sistem_parametr,x_start,tmp_memory1)<0.){return -1;}

double* xnow=sistem_parametr+2+4;

x_alpha[0]=xnow[0];
x_alpha[1]=xnow[3];

return 1.;
}


int derivative_matrix_up(
    const double* x_plus_delta1,
    const double* x_minus_delta1,
    const double* x_plus_delta2,
    const double* x_minus_delta2,
    double delta,
    double* matrix
){
    for(int i=0;i<2;i++){
        matrix[i*2+0]=(x_plus_delta1[i]-x_minus_delta1[i])/(2*delta);
    }

    for(int i=0;i<2;i++){
        matrix[i*2+1]=(x_plus_delta2[i]-x_minus_delta2[i])/(2*delta);
    }

return 1.;
}


void print_matrix(const double* matrix){
    for(int i=0;i<2;i++){
        for(int j=0;j<2;j++){
            fprintf(file2,"%le  ",matrix[i*2+j]);
        }
        fprintf(file2,"\n");
    }
}


double resh_SLAU(double* h,const double* matrix,const double* b){
    h[0]=-(-matrix[3]*b[0]+matrix[1]*b[1])/(matrix[1]*matrix[2]-matrix[0]*matrix[3]);
    h[1]=-(matrix[2]*b[0]-matrix[0]*b[1])/(matrix[1]*matrix[2]-matrix[0]*matrix[3]);
if(Max(fabs(h[0]),fabs(h[1]))>pow(10.,10)){return -1.;}
return 1.;
}

void cappa_vector(double* cappa_v,const double* matrix){
    cappa_v[0]=sqrt(matrix[0]*matrix[0]+matrix[1]*matrix[1]);
    cappa_v[1]=sqrt(matrix[2]*matrix[2]+matrix[3]*matrix[3]);
}

double S_function(const double* x,const double* cappa_v){
    return sqrt(cappa_v[0]*cappa_v[0]*x[0]*x[0]+cappa_v[1]*cappa_v[1]*x[1]*x[1]);
}

double vector_alpha_sampling(
    double* new_vector_alpha,
    const double* vector_alpha,
    const double* h,
    const double* cappa_v,
    double ( *(*F) )(double,const double*),
    double* sistem_parametr,
    double* tmp_memory
)
{
    double* x_alpha=tmp_memory;
    double* tmp_memory1=x_alpha+2;
    double Sum;

    int N=0;
    if(fun_alpha_nevaska(vector_alpha[0]+h[0]*pow(2.,-N),vector_alpha[1]+h[1]*pow(2.,-N),x_alpha,F,sistem_parametr,tmp_memory1)<0.){return -1;}
    N++;

    Sum=S_function(x_alpha,cappa_v);
    for(;N<N_MAX_vector_alpha;N++)
    {
        if(fun_alpha_nevaska(vector_alpha[0]+h[0]*pow(2.,-N),vector_alpha[1]+h[1]*pow(2.,-N),x_alpha,F,sistem_parametr,tmp_memory1)<0.){break;}

        if(S_function(x_alpha,cappa_v)>=Sum)break;
        Sum=S_function(x_alpha,cappa_v);
    }
    new_vector_alpha[0]=vector_alpha[0]+h[0]*pow(2.,-N+1);
    new_vector_alpha[1]=vector_alpha[1]+h[1]*pow(2.,-N+1);

return 1.;
}



void print_vector(double* vector){
    fprintf(file2,"%le   %le\n",vector[0],vector[1]);
}


int RKDOPR_8_7_interval_with_control_main1(
    double px0,
    double py0
){

int dimension=4;

double vector_alpha[2];
vector_alpha[0]=px0;
vector_alpha[1]=py0;
double new_vector_alpha[2];
new_vector_alpha[0]=vector_alpha[0];
new_vector_alpha[1]=vector_alpha[1];

double x_alpha_0[2];
double x_plus_delta1[2];
double x_minus_delta1[2];
double x_plus_delta2[2];
double x_minus_delta2[2];
double matrix[4];

double h[2];
double cappa_v[2];
double delta=pow(10.,-8);

double ( *(F[4]) )(double const,double const*);
F[0]=F_x_main;
F[1]=F_y_main;
F[2]=F_px_main;
F[3]=F_py_main;

double tmp_memory[1000];
double sistem_parametr[1000];

do{
    vector_alpha[0]=new_vector_alpha[0];
    vector_alpha[1]=new_vector_alpha[1];
    
    fprintf(file2,"vector_alpha   ");print_vector(vector_alpha);

    if(fun_alpha_nevaska_printf(vector_alpha[0],vector_alpha[1],x_alpha_0,F,sistem_parametr,tmp_memory)<0.){printf("GLOBAL ERR in find vector nevasok!!!\n");return -1;}
    fprintf(file2,"x_alpha_0   ");print_vector(x_alpha_0);
    if(fun_alpha_nevaska(vector_alpha[0]+delta,vector_alpha[1],x_plus_delta1,F,sistem_parametr,tmp_memory)<0.){printf("GLOBAL ERR in find vector nevasok!!!\n");return -1;}
    fprintf(file2,"x_plus_delta1   ");print_vector(x_plus_delta1);
    if(fun_alpha_nevaska(vector_alpha[0]-delta,vector_alpha[1],x_minus_delta1,F,sistem_parametr,tmp_memory)<0.){printf("GLOBAL ERR in find vector nevasok!!!\n");return -1;}
    fprintf(file2,"x_minus_delta1   ");print_vector(x_minus_delta1);
    if(fun_alpha_nevaska(vector_alpha[0],vector_alpha[1]+delta,x_plus_delta2,F,sistem_parametr,tmp_memory)<0.){printf("GLOBAL ERR in find vector nevasok!!!\n");return -1;}
    fprintf(file2,"x_plus_delta2   ");print_vector(x_plus_delta2);
    if(fun_alpha_nevaska(vector_alpha[0],vector_alpha[1]-delta,x_minus_delta2,F,sistem_parametr,tmp_memory)<0.){printf("GLOBAL ERR in find vector nevasok!!!\n");return -1;}
    fprintf(file2,"x_minus_delta2   ");print_vector(x_minus_delta2);


    derivative_matrix_up(x_plus_delta1,x_minus_delta1,x_plus_delta2,x_minus_delta2,delta,matrix);

    print_matrix(matrix);

    if(resh_SLAU(h,matrix,x_alpha_0)<0.){printf("GLOBAL ERR in find h!!!\n");return -1;}
    fprintf(file2,"h   ");print_vector(h);

    cappa_vector(cappa_v,matrix);
    fprintf(file2,"cappa_v   ");print_vector(cappa_v);

    if(vector_alpha_sampling(new_vector_alpha,vector_alpha,h,cappa_v,F,sistem_parametr,tmp_memory)<0.){printf("GLOBAL ERR in find new h!!!\n");return -1;}
    fprintf(file2,"new_vector_alpha   ");print_vector(new_vector_alpha);

    fprintf(file2,"---\n");
}while(S_function(x_alpha_0,cappa_v)>epsilon);

double x_start[4];
x_start[0]=0.;
x_start[1]=0.;
x_start[2]=vector_alpha[0];
x_start[3]=vector_alpha[1];

RKDOPR_8_7_interval_with_control_printf_B(4,sistem_parametr,x_start,tmp_memory);

// tol_local_err_RK=pow(10.,-8);T=2.5;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-8);T=5;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-8);T=7.5;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-8);T=10;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-10);T=2.5;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-10);T=5;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-10);T=7.5;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-10);T=10;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-12);T=2.5;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-12);T=5;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-12);T=7.5;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-12);T=10;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);

// tol_local_err_RK=pow(10.,-8);T=5;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-8);T=10;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-8);T=15;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-8);T=20;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-10);T=5;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-10);T=10;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-10);T=15;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-10);T=20;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-12);T=5;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-12);T=10;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-12);T=15;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);
// tol_local_err_RK=pow(10.,-12);T=20;RKDOPR_8_7_interval_with_control_err(4,sistem_parametr,x_start,tmp_memory);

double flag=1.;
int N=1000;
T=10;

while(flag>0)
{

T=T-0.1;
N=0;
do{
    N++;
    vector_alpha[0]=new_vector_alpha[0];
    vector_alpha[1]=new_vector_alpha[1];
    
    if(fun_alpha_nevaska_printf(vector_alpha[0],vector_alpha[1],x_alpha_0,F,sistem_parametr,tmp_memory)<0.){flag=-1;break;}
    if(fun_alpha_nevaska(vector_alpha[0]+delta,vector_alpha[1],x_plus_delta1,F,sistem_parametr,tmp_memory)<0.){flag=-1;break;}
    if(fun_alpha_nevaska(vector_alpha[0]-delta,vector_alpha[1],x_minus_delta1,F,sistem_parametr,tmp_memory)<0.){flag=-1;break;}
    if(fun_alpha_nevaska(vector_alpha[0],vector_alpha[1]+delta,x_plus_delta2,F,sistem_parametr,tmp_memory)<0.){flag=-1;break;}
    if(fun_alpha_nevaska(vector_alpha[0],vector_alpha[1]-delta,x_minus_delta2,F,sistem_parametr,tmp_memory)<0.){flag=-1;break;}

    derivative_matrix_up(x_plus_delta1,x_minus_delta1,x_plus_delta2,x_minus_delta2,delta,matrix);

    print_matrix(matrix);

    if(resh_SLAU(h,matrix,x_alpha_0)<0.){flag=-1;break;}

    cappa_vector(cappa_v,matrix);

    if(vector_alpha_sampling(new_vector_alpha,vector_alpha,h,cappa_v,F,sistem_parametr,tmp_memory)<0.){flag=-1;break;}

}while(S_function(x_alpha_0,cappa_v)>epsilon);

if(flag>0){printf("%lf\n",T);}
}

return 0;
}



int RKDOPR_8_7_interval_with_control_main(
    double px0,
    double py0
){

int dimension=4;

double vector_alpha[2];
vector_alpha[0]=px0;
vector_alpha[1]=py0;
double new_vector_alpha[2];
new_vector_alpha[0]=vector_alpha[0];
new_vector_alpha[1]=vector_alpha[1];

double x_alpha_0[2];
double x_plus_delta1[2];
double x_minus_delta1[2];
double x_plus_delta2[2];
double x_minus_delta2[2];
double matrix[4];

double h[2];
double cappa_v[2];
double delta=pow(10.,-8);

double ( *(F[4]) )(double const,double const*);
F[0]=F_x_main;
F[1]=F_y_main;
F[2]=F_px_main;
F[3]=F_py_main;

double tmp_memory[1000];
double sistem_parametr[1000];

do{
    vector_alpha[0]=new_vector_alpha[0];
    vector_alpha[1]=new_vector_alpha[1];
    
    fprintf(file2,"vector_alpha   ");print_vector(vector_alpha);

    if(fun_alpha_nevaska_printf(vector_alpha[0],vector_alpha[1],x_alpha_0,F,sistem_parametr,tmp_memory)<0.){printf("GLOBAL ERR in find vector nevasok!!!\n");return -1;}
    fprintf(file2,"x_alpha_0   ");print_vector(x_alpha_0);
    if(fun_alpha_nevaska(vector_alpha[0]+delta,vector_alpha[1],x_plus_delta1,F,sistem_parametr,tmp_memory)<0.){printf("GLOBAL ERR in find vector nevasok!!!\n");return -1;}
    fprintf(file2,"x_plus_delta1   ");print_vector(x_plus_delta1);
    if(fun_alpha_nevaska(vector_alpha[0]-delta,vector_alpha[1],x_minus_delta1,F,sistem_parametr,tmp_memory)<0.){printf("GLOBAL ERR in find vector nevasok!!!\n");return -1;}
    fprintf(file2,"x_minus_delta1   ");print_vector(x_minus_delta1);
    if(fun_alpha_nevaska(vector_alpha[0],vector_alpha[1]+delta,x_plus_delta2,F,sistem_parametr,tmp_memory)<0.){printf("GLOBAL ERR in find vector nevasok!!!\n");return -1;}
    fprintf(file2,"x_plus_delta2   ");print_vector(x_plus_delta2);
    if(fun_alpha_nevaska(vector_alpha[0],vector_alpha[1]-delta,x_minus_delta2,F,sistem_parametr,tmp_memory)<0.){printf("GLOBAL ERR in find vector nevasok!!!\n");return -1;}
    fprintf(file2,"x_minus_delta2   ");print_vector(x_minus_delta2);


    derivative_matrix_up(x_plus_delta1,x_minus_delta1,x_plus_delta2,x_minus_delta2,delta,matrix);

    print_matrix(matrix);

    if(resh_SLAU(h,matrix,x_alpha_0)<0.){printf("GLOBAL ERR in find h!!!\n");return -1;}
    fprintf(file2,"h   ");print_vector(h);

    cappa_vector(cappa_v,matrix);
    fprintf(file2,"cappa_v   ");print_vector(cappa_v);

    if(vector_alpha_sampling(new_vector_alpha,vector_alpha,h,cappa_v,F,sistem_parametr,tmp_memory)<0.){printf("GLOBAL ERR in find new h!!!\n");return -1;}
    fprintf(file2,"new_vector_alpha   ");print_vector(new_vector_alpha);

    fprintf(file2,"---\n");
}while(S_function(x_alpha_0,cappa_v)>epsilon);

return 0;
}