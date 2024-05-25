#include<cstdio>
#include<iostream>
#include<cmath>


double V[1000000];
double Vup[1000000];
double V_right[1000000];

void write_vector_file(FILE *file,double* vector,size_t size_vector)
{
    for(size_t i=0;i<size_vector;i++)
    {
        fprintf(file,"%e  ",vector[i]);
    }fprintf(file,"\n");
return;
}

void equalize_vector(double* vector1,const double* vector2,size_t size_vector)
{
    for(size_t i=0;i<size_vector;i++)
    {
        vector1[i]=vector2[i];
    }
return;
}

double Vright_tmp[1000000];
double Vup_cr[1000000];
double Vup_nx[1000000];
double Vup_tmp[1000000];
double alpha[1000000];
double beta[1000000];

double F(double v){return -v*v/2;}
double F_d(double v){return -v;}

int gn_sweep(int N,double* V_right,double* Vup_g)
{
    double A[N];
    double B[N];

    A[0] = 0.;
    A[1] = -beta[0];

    B[0] = 0;
    B[1] = V_right[0];

    for(int i=2;i<N;i++)
    {
        A[i]=-beta[i-1]/(alpha[i-1]*A[i-1]+1.);
        B[i]=(V_right[i-1]-alpha[i-1]*B[i-1])/(alpha[i-1]*A[i-1]+1.);
    }

    Vup_g[N-1]=(V_right[N-1]-B[N-1]*alpha[N-1])/(A[N-1]*alpha[N-1]+1.);

    for(int i=N-2;i>=0;i--)
    {
        Vup_g[i]=A[i+1]*Vup_g[i+1]+B[i+1];
    }
return 0;
}


int method_Newton(int Ntau,int Mh,double tol,double w)
{
    double tau=1./Ntau;
    double h=2./Mh;
    double nu=tau/h;
    double Tol=10;

    equalize_vector(Vup_nx,V+1,Mh-1);

    while(Tol>tol)
    {   
        equalize_vector(Vup_cr,Vup_nx,Mh-1);
        
        V_right[0]=Vup_cr[0]+nu/2.*F(Vup_cr[1])-V[1]-nu/2.*F(Vup[0])-w/tau*(V[0]-2*V[1]+V[2]);
        for(int i=1;i<Mh-2;i++)
        {
            V_right[i]=Vup_cr[i]-nu/2*F(Vup_cr[i-1])-V[i+1]+nu/2*F(Vup_cr[i+1])-w/tau*(V[i]-2*V[i+1]+V[i+2]);
        }
        V_right[Mh-2]=Vup_cr[Mh-2]-nu/2*F(Vup_cr[Mh-3])-V[Mh-1]+nu/2*F(Vup[Mh])-w/tau*(V[Mh-2]-2*V[Mh-1]+V[Mh]);
        
        for(int i=0;i<Mh-1;i++)V_right[i]=-V_right[i];

        alpha[0]=0;
        for(int i=0;i<Mh-1;i++){alpha[i]=-nu/2*F_d(Vup_cr[i-1]);}

        for(int i=0;i<Mh-2;i++){beta[i]=nu/2*F_d(Vup_cr[i+1]);}
        beta[Mh-2]=0;

        gn_sweep(Mh-1,V_right,Vup_nx);

        double capha[100000];
        capha[0]=sqrt(1+pow(beta[0],2));
        for(int i=1;i<Mh-2;i++){
        capha[i]=sqrt(1+pow(beta[i],2)+pow(alpha[i],2));
        }capha[Mh-2]=sqrt(1+pow(alpha[Mh-2],2));

        Tol=0;
        double Tol_tmp=0;
        double Tol_min=10000000;
        int K=-1;
        for(int k=0;k<12;k++)
        {
            Tol_tmp=0;
            for(int i=0;i<Mh-1;i++)
            {
                Vup_tmp[i]=Vup_cr[i]+Vup_nx[i]/pow(2,k);
            }

            Vright_tmp[0]=Vup_tmp[0]+nu/2.*F(Vup_tmp[1])-V[1]-nu/2.*F(Vup[0])-w/tau*(V[0]-2*V[1]+V[2]);
            for(int i=1;i<Mh-2;i++)
            {
                V_right[i]=Vup_tmp[i]-nu/2*F(Vup_tmp[i-1])-V[i+1]+nu/2*F(Vup_tmp[i+1])-w/tau*(V[i]-2*V[i+1]+V[i+2]);
            }
            V_right[Mh-2]=Vup_tmp[Mh-2]-nu/2*F(Vup_tmp[Mh-3])-V[Mh-1]+nu/2*F(Vup[Mh])-w/tau*(V[Mh-2]-2*V[Mh-1]+V[Mh]);

            for(int i=0;i<Mh-1;i++)
            {
                Tol_tmp=Tol_tmp+pow(capha[i]*Vright_tmp[i],2);
            }

            if(Tol_tmp<Tol_min)
            {
                Tol_min=Tol_tmp;
                K=k;
            }
        }


        for(int i=0;i<Mh-1;i++)
        {
            Vup_nx[i]=Vup_cr[i]+Vup_nx[i]/pow(2,K);
        }

        Vright_tmp[0]=Vup_nx[0]+nu/2.*F(Vup_nx[1])-V[1]-nu/2.*F(Vup[0])-w/tau*(V[0]-2*V[1]+V[2]);
        for(int i=1;i<Mh-2;i++)
        {
            V_right[i]=Vup_nx[i]-nu/2*F(Vup_nx[i-1])-V[i+1]+nu/2*F(Vup_nx[i+1])-w/tau*(V[i]-2*V[i+1]+V[i+2]);
        }
        V_right[Mh-2]=Vup_nx[Mh-2]-nu/2*F(Vup_nx[Mh-3])-V[Mh-1]+nu/2*F(Vup[Mh])-w/tau*(V[Mh-2]-2*V[Mh-1]+V[Mh]);

        for(int i=0;i<Mh-1;i++)
        {
            Tol=Tol+pow(capha[i]*Vright_tmp[i],2);
        }

    printf("%f \n",Tol);
    }
    equalize_vector(Vup+1,Vup_nx,Mh-1);
return 0;
}

double nonlinear_solution(double t,double x)
{
    return -x*(1./2);
}



int calculate_nonlinear_implicit_write_file(FILE* file,FILE* file_table, int Ntau, int Mh,double tol,double w)
{
    double h=2./Mh;
    double tau=1./Ntau;
    double nu=tau/h;
        
    for(int i=0;i<Mh+1;i++)
    {
        V[i]=-(-1.+i*h);
    }        

    for(int j=1;j<Ntau+1;j++)
    {
        Vup[0]=1./(1+tau*j);
        Vup[Mh]=-1./(1+tau*j);
        
        method_Newton(Ntau,Mh,tol,w);

        equalize_vector(V,Vup,Mh+1);
    }

    double Delta_V_Ch=0.;
    double V_Ch=0.;

    double Delta_V_L1h=0.;
    double V_L1h=0.;

    for(int i=0;i<Mh+1;i++)
    {
        if( fabs( Vup[i]-nonlinear_solution(1.,-1.+i*h) )>Delta_V_Ch )
            {
                Delta_V_Ch=fabs( Vup[i]-nonlinear_solution(1.,-1.+i*h) );
            }
        if( fabs( Vup[i] )>V_Ch )
            {
                V_Ch=fabs( Vup[i] );
            }
        Delta_V_L1h=Delta_V_L1h+fabs( Vup[i]-nonlinear_solution(1.,-1.+i*h) );
        V_L1h=V_L1h+fabs( Vup[i] );
    }
    Delta_V_L1h=Delta_V_L1h*h;
    V_L1h=V_L1h*h;

    for(int i=0;i<Mh+1;i++){Vup[i]=-1+i*h;}

    write_vector_file(file,Vup,Mh+1);
    write_vector_file(file,V,Mh+1);

    fprintf(file_table,"%.3f %.3f %e %e %e %e \n",tau,h,Delta_V_Ch,Delta_V_L1h,Delta_V_Ch/V_Ch,Delta_V_L1h/V_L1h);
return 0;

}




int main(void)
{
    FILE* file1=fopen("t1.txt","w");
    FILE* file2=fopen("t1_t.txt","w");

    calculate_nonlinear_implicit_write_file(file1,file2,100,2000,pow(10,-2),0.1);

return 0;
}
