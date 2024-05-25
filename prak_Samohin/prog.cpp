#include<cstdio>
#include<iostream>
#include<cmath>


double V[1000000];
double Vup[1000000];
double V_tmp[1000000];
double Vup_tmp[1000000];


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




double linear_solution(double t,double x)
{
    if(t>0. && 2*x-t>0.)return 1.;
    else if(t>0. && 2*x-t<=0.)return 0.;
return 0.;
}




int calculate_linear_explicit_write_file(FILE* file,FILE* file_table, int Ntau, int Mh)
{
    double h=2./Mh;
    double tau=1./Ntau;
    double nu=tau/h;
    
    for(int i=0;i<Mh+1;i++)
    {
        if(-1.+i*h>0.)V[i]=1.;
        else V[i]=0.;
    }V[Mh+1]=0.;
    
    for(int j=1;j<Ntau+1;j++)
    {
        Vup[0]=0.;
        Vup[Mh]=1.;
        Vup[Mh+1]=tau*j;

        for(int i=1;i<Mh;i++)
        {
            Vup[i]=V[i]-nu/4*(V[i+1]-V[i-1])+nu*nu/8*(V[i+1]+V[i-1]-2*V[i]);

        }

        equalize_vector(V,Vup,Mh+2);
    }

    double Delta_V_Ch=0.;
    double V_Ch=0.;

    double Delta_V_L1h=0.;
    double V_L1h=0.;

    for(int i=0;i<Mh+1;i++)
    {
        if( fabs( Vup[i]-linear_solution(Ntau*tau,-1.+i*h) )>Delta_V_Ch )
            {
                Delta_V_Ch=fabs( Vup[i]-linear_solution(Ntau*tau,-1.+i*h) );
            }
        if( fabs( Vup[i] )>V_Ch )
            {
                V_Ch=fabs( Vup[i] );
            }
        Delta_V_L1h=Delta_V_L1h+fabs( Vup[i]-linear_solution(Ntau*tau,-1.+i*h) );
        V_L1h=V_L1h+fabs( Vup[i] );
    }
    Delta_V_L1h=Delta_V_L1h*h;
    V_L1h=V_L1h*h;

    for(int i=0;i<Mh+1;i++){Vup[i]=-1+i*h;}

    write_vector_file(file,Vup,Mh+1);
    write_vector_file(file,V,Mh+1);

    fprintf(file_table,"%f %f %e %e %e %e \n",tau,h,Delta_V_Ch,Delta_V_L1h,Delta_V_Ch/V_Ch,Delta_V_L1h/V_L1h);
return 0;
}




int calculate_linear_explicit_write_file_local(FILE* file_table, int Ntau, int Mh,int K)
{
    double h=2./Mh;
    double tau=1./Ntau;
    double nu=tau/h;
    
    for(int i=0;i<Mh+1;i++)
    {
        if(-1.+i*h>0.)V[i]=1.;
        else V[i]=0.;
    }V[Mh+1]=0.;
    
    for(int j=1;j<Ntau+1;j++)
    {
        Vup[0]=0.;
        Vup[Mh]=1.;
        Vup[Mh+1]=tau*j;

        for(int i=1;i<Mh;i++)
        {
            Vup[i]=V[i]-nu/4*(V[i+1]-V[i-1])+nu*nu/8*(V[i+1]+V[i-1]-2*V[i]);

        }

        equalize_vector(V,Vup,Mh+2);
    }

    int Mh_tmp=Mh;
    int Ntau_tmp=Ntau;
    int n=1;

    for(int k=0;k<K;k++)
    {
    
    Mh_tmp=Mh_tmp*2;
    Ntau_tmp=Ntau_tmp*2;
    n=n*2;

    double h_tmp=2./Mh_tmp;
    double tau_tmp=1./Ntau_tmp;
    double nu_tmp=tau_tmp/h_tmp;

    for(int i=0;i<Mh_tmp+1;i++)
    {
        if(-1.+i*h_tmp>0.)V_tmp[i]=1.;
        else V_tmp[i]=0.;
    }V_tmp[Mh_tmp+1]=0.;
    
    for(int j=1;j<Ntau_tmp+1;j++)
    {
        Vup_tmp[0]=0.;
        Vup_tmp[Mh_tmp]=1.;
        Vup_tmp[Mh_tmp+1]=tau_tmp*j;

        for(int i=1;i<Mh_tmp;i++)
        {
            Vup_tmp[i]=V_tmp[i]-nu_tmp/4*(V_tmp[i+1]-V_tmp[i-1])+nu_tmp*nu_tmp/8*(V_tmp[i+1]+V_tmp[i-1]-2*V_tmp[i]);

        }

        equalize_vector(V_tmp,Vup_tmp,Mh_tmp+2);
    }


    double Delta_V_Ch=0.;
    double V_Ch=0.;

    double Delta_V_L1h=0.;
    double V_L1h=0.;
    
    for(int i=0;i<Mh+1;i++)
    {
        if( fabs( V[i]-V_tmp[i*n] )>Delta_V_Ch )
            {
                Delta_V_Ch=fabs( V[i]-V_tmp[i*n] );
            }
        if( fabs( V[i] )>V_Ch )
            {
                V_Ch=fabs( V[i] );
            }
        Delta_V_L1h=Delta_V_L1h+fabs( V[i]-V_tmp[i*n] );
        V_L1h=V_L1h+fabs( V[i] );
    }
    Delta_V_L1h=Delta_V_L1h*h;
    V_L1h=V_L1h*h;

    fprintf(file_table,"%f %f %e %e %e %e \n",tau_tmp,h_tmp,Delta_V_Ch,Delta_V_L1h,Delta_V_Ch/V_Ch,Delta_V_L1h/V_L1h);

    }

return 0;
}



int calculate_linear_implicit_write_file(FILE* file,FILE* file_table, int Ntau, int Mh)
{
    double h=2./Mh;
    double tau=1./Ntau;
    double nu=tau/h;
    
    for(int i=0;i<Mh+1;i++)
    {
        if(-1.+i*h>0.)V[i]=1.;
        else V[i]=0.;
    }V[Mh+1]=0.;
    
    for(int j=1;j<Ntau+1;j++)
    {
        Vup[0]=0.;
        Vup[Mh]=1.;
        Vup[Mh+1]=tau*j;

        for(int i=Mh-1;i>0;i--)
        {
            Vup[i]=(4.*V[i]-nu*(V[i]-V[i-1]+Vup[i+1]))/(4-nu);
        }
        equalize_vector(V,Vup,Mh+2);
    }

    double Delta_V_Ch=0.;
    double V_Ch=0.;

    double Delta_V_L1h=0.;
    double V_L1h=0.;

    for(int i=0;i<Mh+1;i++)
    {
        if( fabs( Vup[i]-linear_solution(Ntau*tau,-1.+i*h) )>Delta_V_Ch )
            {
                Delta_V_Ch=fabs( Vup[i]-linear_solution(Ntau*tau,-1.+i*h) );
            }
        if( fabs( Vup[i] )>V_Ch )
            {
                V_Ch=fabs( Vup[i] );
            }
        Delta_V_L1h=Delta_V_L1h+fabs( Vup[i]-linear_solution(Ntau*tau,-1.+i*h) );
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



int calculate_linear_implicit_write_file_local(FILE* file_table, int Ntau, int Mh,int K)
{
    double h=2./Mh;
    double tau=1./Ntau;
    double nu=tau/h;
    
    for(int i=0;i<Mh+1;i++)
    {
        if(-1.+i*h>0.)V[i]=1.;
        else V[i]=0.;
    }V[Mh+1]=0.;
    
    for(int j=1;j<Ntau+1;j++)
    {
        Vup[0]=0.;
        Vup[Mh]=1.;
        Vup[Mh+1]=tau*j;

        for(int i=Mh-1;i>0;i--)
        {
            Vup[i]=(4.*V[i]-nu*(V[i]-V[i-1]+Vup[i+1]))/(4-nu);
        }
        equalize_vector(V,Vup,Mh+2);
    }

    int Mh_tmp=Mh;
    int Ntau_tmp=Ntau;
    int n=1;

    for(int k=0;k<K;k++)
    {
    
    Mh_tmp=Mh_tmp*2;
    Ntau_tmp=Ntau_tmp*2;
    n=n*2;

    double h_tmp=2./Mh_tmp;
    double tau_tmp=1./Ntau_tmp;
    double nu_tmp=tau_tmp/h_tmp;

    for(int i=0;i<Mh_tmp+1;i++)
    {
        if(-1.+i*h_tmp>0.)V_tmp[i]=1.;
        else V_tmp[i]=0.;
    }V_tmp[Mh_tmp+1]=0.;
    
    for(int j=1;j<Ntau_tmp+1;j++)
    {
        Vup_tmp[0]=0.;
        Vup_tmp[Mh]=1.;
        Vup_tmp[Mh+1]=tau_tmp*j;

        for(int i=Mh_tmp-1;i>0;i--)
        {
            Vup_tmp[i]=(4.*V_tmp[i]-nu_tmp*(V_tmp[i]-V_tmp[i-1]+Vup_tmp[i+1]))/(4-nu_tmp);
        }
        equalize_vector(V_tmp,Vup_tmp,Mh_tmp+2);
    }


    double Delta_V_Ch=0.;
    double V_Ch=0.;

    double Delta_V_L1h=0.;
    double V_L1h=0.;
    
    for(int i=0;i<Mh+1;i++)
    {
        if( fabs( V[i]-V_tmp[i*n] )>Delta_V_Ch )
            {
                Delta_V_Ch=fabs( V[i]-V_tmp[i*n] );
            }
        if( fabs( V[i] )>V_Ch )
            {
                V_Ch=fabs( V[i] );
            }
        Delta_V_L1h=Delta_V_L1h+fabs( V[i]-V_tmp[i*n] );
        V_L1h=V_L1h+fabs( V[i] );
    }
    Delta_V_L1h=Delta_V_L1h*h;
    V_L1h=V_L1h*h;

    fprintf(file_table,"%f %f %e %e %e %e \n",tau_tmp,h_tmp,Delta_V_Ch,Delta_V_L1h,Delta_V_Ch/V_Ch,Delta_V_L1h/V_L1h);

    }

return 0;
}

double nonlinear_solution(double t,double x)
{
    if(t>0. && x-t>0.)return 1.;
    else if(t>0. && x<=0.)return 0.;
    else return x/t;
return 0.;
}


int calculate_nonlinear_explicit_write_file(FILE* file,FILE* file_table, int Ntau, int Mh)
{
    double h=2./Mh;
    double tau=1./Ntau;
    double nu=tau/h;
    
    for(int i=0;i<Mh+1;i++)
    {
        if(-1.+i*h>0.)V[i]=1.;
        else V[i]=0.;
    }V[Mh+1]=0.;
    
    for(int j=1;j<Ntau+1;j++)
    {
        Vup[0]=0.;
        Vup[Mh]=1.;
        Vup[Mh+1]=tau*j;
        double V1,V2;

        for(int i=1;i<Mh;i++)
        {
            // Vup[i]=V[i] - nu/4* ( pow(V[i+1],2) - pow(V[i],2) )   -   nu/4. *( V[i] - V[i-1] -nu/2.* ( pow(V[i+1],2)-2*pow(V[i],2)+pow(V[i-1],2) ))*( V[i]+V[i-1]-nu/2.*(pow(V[i+1],2)-pow(V[i-1],2)) );

            V1=V[i]-nu/2*(pow(V[i+1],2)-pow(V[i],2));
            V2=V[i-1]-nu/2*(pow(V[i],2)-pow(V[i-1],2));

            Vup[i]=V[i]-nu/4*(pow(V[i+1],2)-pow(V[i],2))-nu/4*(pow(V1,2)-pow(V2,2));
        }

        equalize_vector(V,Vup,Mh+2);
    }

    double Delta_V_Ch=0.;
    double V_Ch=0.;

    double Delta_V_L1h=0.;
    double V_L1h=0.;

    for(int i=0;i<Mh+1;i++)
    {
        if( fabs( Vup[i]-nonlinear_solution(Ntau*tau,-1.+i*h) )>Delta_V_Ch )
            {
                Delta_V_Ch=fabs( Vup[i]-nonlinear_solution(Ntau*tau,-1.+i*h) );
            }
        if( fabs( Vup[i] )>V_Ch )
            {
                V_Ch=fabs( Vup[i] );
            }
        Delta_V_L1h=Delta_V_L1h+fabs( Vup[i]-nonlinear_solution(Ntau*tau,-1.+i*h) );
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




int calculate_nonlinear_explicit_write_file_local(FILE* file_table, int Ntau, int Mh,int K)
{
    double h=2./Mh;
    double tau=1./Ntau;
    double nu=tau/h;
    
    for(int i=0;i<Mh+1;i++)
    {
        if(-1.+i*h>0.)V[i]=1.;
        else V[i]=0.;
    }V[Mh+1]=0.;
    
    for(int j=1;j<Ntau+1;j++)
    {
        Vup[0]=0.;
        Vup[Mh]=1.;
        Vup[Mh+1]=tau*j;

        for(int i=1;i<Mh;i++)
        {
            Vup[i]=V[i] - nu/4* ( pow(V[i+1],2) - pow(V[i],2) )   -   nu/4. *( V[i] - V[i-1] -nu/2.* ( pow(V[i+1],2)-2*pow(V[i],2)+pow(V[i-1],2) ))*( V[i]+V[i-1]-nu/2.*(pow(V[i+1],2)-pow(V[i-1],2)) );
        }
        equalize_vector(V,Vup,Mh+2);
    }

    int Mh_tmp=Mh;
    int Ntau_tmp=Ntau;
    int n=1;

    for(int k=0;k<K;k++)
    {
    
    Mh_tmp=Mh_tmp*2;
    Ntau_tmp=Ntau_tmp*2;
    n=n*2;

    double h_tmp=2./Mh_tmp;
    double tau_tmp=1./Ntau_tmp;
    double nu_tmp=tau_tmp/h_tmp;

    for(int i=0;i<Mh_tmp+1;i++)
    {
        if(-1.+i*h_tmp>0.)V_tmp[i]=1.;
        else V_tmp[i]=0.;
    }V_tmp[Mh_tmp+1]=0.;
    
    for(int j=1;j<Ntau_tmp+1;j++)
    {
        Vup_tmp[0]=0.;
        Vup_tmp[Mh_tmp]=1.;
        Vup_tmp[Mh_tmp+1]=tau_tmp*j;

        for(int i=1;i<Mh_tmp;i++)
        {
            Vup_tmp[i]=V_tmp[i] - nu_tmp/4* ( pow(V_tmp[i+1],2) - pow(V_tmp[i],2) )   -   nu_tmp/4. *( V_tmp[i] - V_tmp[i-1] -nu_tmp/2.* ( pow(V_tmp[i+1],2)-2*pow(V_tmp[i],2)+pow(V_tmp[i-1],2) ))*( V_tmp[i]+V_tmp[i-1]-nu_tmp/2.*(pow(V_tmp[i+1],2)-pow(V_tmp[i-1],2)) );
        }

        equalize_vector(V_tmp,Vup_tmp,Mh_tmp+2);
    }


    double Delta_V_Ch=0.;
    double V_Ch=0.;

    double Delta_V_L1h=0.;
    double V_L1h=0.;
    
    for(int i=0;i<Mh+1;i++)
    {
        if( fabs( V[i]-V_tmp[i*n] )>Delta_V_Ch )
            {
                Delta_V_Ch=fabs( V[i]-V_tmp[i*n] );
            }
        if( fabs( V[i] )>V_Ch )
            {
                V_Ch=fabs( V[i] );
            }
        Delta_V_L1h=Delta_V_L1h+fabs( V[i]-V_tmp[i*n] );
        V_L1h=V_L1h+fabs( V[i] );
    }
    Delta_V_L1h=Delta_V_L1h*h;
    V_L1h=V_L1h*h;

    fprintf(file_table,"%f %f %e %e %e %e \n",tau_tmp,h_tmp,Delta_V_Ch,Delta_V_L1h,Delta_V_Ch/V_Ch,Delta_V_L1h/V_L1h);

    }

return 0;
}



int calculate_nonlinear_implicit_write_file(FILE* file,FILE* file_table, int Ntau, int Mh)
{
    double h=2./Mh;
    double tau=1./Ntau;
    double nu=tau/h;
    
    for(int i=0;i<Mh+1;i++)
    {
        if(-1.+i*h>0.)V[i]=1.;
        else V[i]=0.;
    }V[Mh+1]=0.;
    
    for(int j=1;j<Ntau+1;j++)
    {
        Vup[0]=0.;
        Vup[Mh]=1.;
        Vup[Mh+1]=tau*j;

        for(int i=Mh-1;i>0;i--)
        {
            if(4./pow(nu,2)-4./nu*V[i]+pow(V[i],2)+pow(Vup[i+1],2)-pow(V[i-1],2)>=0.)
            {
                if(2./nu-V[i] >=0.)
                {
                    Vup[i]=2./nu-sqrt(4./pow(nu,2)-4./nu*V[i]+pow(V[i],2)+pow(Vup[i+1],2)-pow(V[i-1],2));
                }
                else
                {
                    Vup[i]=2./nu+sqrt(4./pow(nu,2)-4./nu*V[i]+pow(V[i],2)+pow(Vup[i+1],2)-pow(V[i-1],2));
                }
            }
            else
            Vup[i]=2./nu;
        }
        equalize_vector(V,Vup,Mh+2);
    }

    double Delta_V_Ch=0.;
    double V_Ch=0.;

    double Delta_V_L1h=0.;
    double V_L1h=0.;

    for(int i=0;i<Mh+1;i++)
    {
        if( fabs( Vup[i]-linear_solution(Ntau*tau,-1.+i*h) )>Delta_V_Ch )
            {
                Delta_V_Ch=fabs( Vup[i]-linear_solution(Ntau*tau,-1.+i*h) );
            }
        if( fabs( Vup[i] )>V_Ch )
            {
                V_Ch=fabs( Vup[i] );
            }
        Delta_V_L1h=Delta_V_L1h+fabs( Vup[i]-linear_solution(Ntau*tau,-1.+i*h) );
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



int calculate_nonlinear_implicit_write_file_local(FILE* file_table, int Ntau, int Mh,int K)
{
    double h=2./Mh;
    double tau=1./Ntau;
    double nu=tau/h;
    
    for(int i=0;i<Mh+1;i++)
    {
        if(-1.+i*h>0.)V[i]=1.;
        else V[i]=0.;
    }V[Mh+1]=0.;
    
    for(int j=1;j<Ntau+1;j++)
    {
        Vup[0]=0.;
        Vup[Mh]=1.;
        Vup[Mh+1]=tau*j;

        for(int i=Mh-1;i>0;i--)
        {
            if(4./pow(nu,2)-4./nu*V[i]+pow(V[i],2)+pow(Vup[i+1],2)-pow(V[i-1],2)>=0.)
            {
                if(2./nu-V[i] >=0.)
                {
                    Vup[i]=2./nu-sqrt(4./pow(nu,2)-4./nu*V[i]+pow(V[i],2)+pow(Vup[i+1],2)-pow(V[i-1],2));
                }
                else
                {
                    Vup[i]=2./nu+sqrt(4./pow(nu,2)-4./nu*V[i]+pow(V[i],2)+pow(Vup[i+1],2)-pow(V[i-1],2));
                }
            }
            else
            Vup[i]=2./nu;
        }
        equalize_vector(V,Vup,Mh+2);
    }

    int Mh_tmp=Mh;
    int Ntau_tmp=Ntau;
    int n=1;

    for(int k=0;k<K;k++)
    {
    
    Mh_tmp=Mh_tmp*2;
    Ntau_tmp=Ntau_tmp*2;
    n=n*2;

    double h_tmp=2./Mh_tmp;
    double tau_tmp=1./Ntau_tmp;
    double nu_tmp=tau_tmp/h_tmp;

    for(int i=0;i<Mh_tmp+1;i++)
    {
        if(-1.+i*h_tmp>0.)V_tmp[i]=1.;
        else V_tmp[i]=0.;
    }V_tmp[Mh_tmp+1]=0.;
    
    for(int j=1;j<Ntau_tmp+1;j++)
    {
        Vup_tmp[0]=0.;
        Vup_tmp[Mh]=1.;
        Vup_tmp[Mh+1]=tau_tmp*j;

        for(int i=Mh_tmp-1;i>0;i--)
        {
            if(4./pow(nu_tmp,2)-4./nu_tmp*V_tmp[i]+pow(V_tmp[i],2)+pow(Vup_tmp[i+1],2)-pow(V_tmp[i-1],2)>=0.)
            {
                if(2./nu_tmp-V_tmp[i] >=0.)
                {
                    Vup_tmp[i]=2./nu_tmp-sqrt(4./pow(nu_tmp,2)-4./nu_tmp*V_tmp[i]+pow(V_tmp[i],2)+pow(Vup_tmp[i+1],2)-pow(V_tmp[i-1],2));
                }
                else
                {
                    Vup_tmp[i]=2./nu_tmp+sqrt(4./pow(nu_tmp,2)-4./nu_tmp*V_tmp[i]+pow(V_tmp[i],2)+pow(Vup_tmp[i+1],2)-pow(V_tmp[i-1],2));
                }
            }
            else
            Vup_tmp[i]=2./nu_tmp;
        }
        equalize_vector(V_tmp,Vup_tmp,Mh_tmp+2);
    }


    double Delta_V_Ch=0.;
    double V_Ch=0.;

    double Delta_V_L1h=0.;
    double V_L1h=0.;
    
    for(int i=0;i<Mh+1;i++)
    {
        if( fabs( V[i]-V_tmp[i*n] )>Delta_V_Ch )
            {
                Delta_V_Ch=fabs( V[i]-V_tmp[i*n] );
            }
        if( fabs( V[i] )>V_Ch )
            {
                V_Ch=fabs( V[i] );
            }
        Delta_V_L1h=Delta_V_L1h+fabs( V[i]-V_tmp[i*n] );
        V_L1h=V_L1h+fabs( V[i] );
    }
    Delta_V_L1h=Delta_V_L1h*h;
    V_L1h=V_L1h*h;

    fprintf(file_table,"%f %f %e %e %e %e \n",tau_tmp,h_tmp,Delta_V_Ch,Delta_V_L1h,Delta_V_Ch/V_Ch,Delta_V_L1h/V_L1h);

    }

return 0;
}


//-----------------------------------


double calculate_nonlinear_explicit_write_file_tetha(FILE* file,FILE* file_table, int Ntau, int Mh,double tetha)
{
    double h=2./Mh;
    double tau=1./Ntau;
    double nu=tau/h;
    
    for(int i=0;i<Mh+1;i++)
    {
        if(-1.+i*h > tetha )V[i]=1.;
        else if(-1.+i*h>0. && -1.+i*h<=tetha)V[i]=(-1.+i*h)/tetha;
        else V[i]=0.;
    }V[Mh+1]=0.;
    
    for(int j=1;j<Ntau+1;j++)
    {
        Vup[0]=0.;

        if(tau*j<1.-tetha){Vup[Mh]=1.;}
        else{Vup[Mh]=1./(tau*j+tetha);}

        Vup[Mh+1]=tau*j;
        double V1,V2;

        for(int i=1;i<Mh;i++)
        {
            V1=V[i]-nu/2*(pow(V[i+1],2)-pow(V[i],2));
            V2=V[i-1]-nu/2*(pow(V[i],2)-pow(V[i-1],2));

            Vup[i]=V[i]-nu/4*(pow(V[i+1],2)-pow(V[i],2))-nu/4*(pow(V1,2)-pow(V2,2));
        }
        equalize_vector(V,Vup,Mh+2);
    }

    double Delta_V_Ch=0.;
    double V_Ch=0.;

    double Delta_V_L1h=0.;
    double V_L1h=0.;

    for(int i=0;i<Mh+1;i++)
    {
        if( fabs( Vup[i]-nonlinear_solution(Ntau*tau,-1.+i*h) )>Delta_V_Ch )
            {
                Delta_V_Ch=fabs( Vup[i]-nonlinear_solution(Ntau*tau,-1.+i*h) );
            }
        if( fabs( Vup[i] )>V_Ch )
            {
                V_Ch=fabs( Vup[i] );
            }
        Delta_V_L1h=Delta_V_L1h+fabs( Vup[i]-nonlinear_solution(Ntau*tau,-1.+i*h) );
        V_L1h=V_L1h+fabs( Vup[i] );
    }
    Delta_V_L1h=Delta_V_L1h*h;
    V_L1h=V_L1h*h;

    for(int i=0;i<Mh+1;i++){Vup[i]=-1+i*h;}

    // write_vector_file(file,Vup,Mh+1);
    // write_vector_file(file,V,Mh+1);
    fprintf(file_table,"%.3f %.3f %e %e %e %e %e \n",tau,h,tetha,Delta_V_Ch,Delta_V_L1h,Delta_V_Ch/V_Ch,Delta_V_L1h/V_L1h);
return Delta_V_Ch;
}

double calculate_nonlinear_implicit_write_file_tetha(FILE* file,FILE* file_table, int Ntau, int Mh,double tetha)
{
    double h=2./Mh;
    double tau=1./Ntau;
    double nu=tau/h;
    
    for(int i=0;i<Mh+1;i++)
    {
        if(-1.+i*h > tetha )V[i]=1.;
        else if(-1.+i*h>0. && -1.+i*h<=tetha)V[i]=(-1.+i*h)/tetha;
        else V[i]=0.;
    }V[Mh+1]=0.;
    
    for(int j=1;j<Ntau+1;j++)
    {
        Vup[0]=0.;

        if(tau*j<1.-tetha){Vup[Mh]=1.;}
        else{Vup[Mh]=1./(tau*j+tetha);}

        for(int i=Mh-1;i>0;i--)
        {
            if(4./pow(nu,2)-4./nu*V[i]+pow(V[i],2)+pow(Vup[i+1],2)-pow(V[i-1],2)>=0.)
            {
                if(2./nu-V[i] >=0.)
                {
                    Vup[i]=2./nu-sqrt(4./pow(nu,2)-4./nu*V[i]+pow(V[i],2)+pow(Vup[i+1],2)-pow(V[i-1],2));
                }
                else
                {
                    Vup[i]=2./nu+sqrt(4./pow(nu,2)-4./nu*V[i]+pow(V[i],2)+pow(Vup[i+1],2)-pow(V[i-1],2));
                }
            }
            else
            Vup[i]=2./nu;
        }
        equalize_vector(V,Vup,Mh+2);
    }

    double Delta_V_Ch=0.;
    double V_Ch=0.;

    double Delta_V_L1h=0.;
    double V_L1h=0.;

    for(int i=0;i<Mh+1;i++)
    {
        if( fabs( Vup[i]-linear_solution(Ntau*tau,-1.+i*h) )>Delta_V_Ch )
            {
                Delta_V_Ch=fabs( Vup[i]-linear_solution(Ntau*tau,-1.+i*h) );
            }
        if( fabs( Vup[i] )>V_Ch )
            {
                V_Ch=fabs( Vup[i] );
            }
        Delta_V_L1h=Delta_V_L1h+fabs( Vup[i]-linear_solution(Ntau*tau,-1.+i*h) );
        V_L1h=V_L1h+fabs( Vup[i] );
    }
    Delta_V_L1h=Delta_V_L1h*h;
    V_L1h=V_L1h*h;

    for(int i=0;i<Mh+1;i++){Vup[i]=-1+i*h;}

    // write_vector_file(file,Vup,Mh+1);
    // write_vector_file(file,V,Mh+1);

    fprintf(file_table,"%.3f %.3f %e %e %e %e %e \n",tau,h,tetha,Delta_V_Ch,Delta_V_L1h,Delta_V_Ch/V_Ch,Delta_V_L1h/V_L1h);
return Delta_V_Ch;
}


int main(void)
{
    FILE* file1=fopen("t1.txt","w");
    FILE* file2=fopen("t2.txt","w");
    FILE* file3=fopen("t3.txt","w");
    FILE* file4=fopen("t4.txt","w");
    FILE* file5=fopen("t5.txt","w");
    FILE* file6=fopen("t6.txt","w");
    FILE* file7=fopen("t7.txt","w");
    FILE* file8=fopen("t8.txt","w");

    FILE* File1=fopen("t_1.txt","w");
    FILE* File2=fopen("t_2.txt","w");
    FILE* File3=fopen("t_3.txt","w");
    FILE* File4=fopen("t_4.txt","w");
    FILE* File5=fopen("t_5.txt","w");
    FILE* File6=fopen("t_6.txt","w");
    FILE* File7=fopen("t_7.txt","w");
    FILE* File8=fopen("t_8.txt","w");
    
    FILE* File9=fopen("t_9.txt","w");
    FILE* File10=fopen("t_10.txt","w");
    FILE* File11=fopen("t_11.txt","w");
    FILE* File12=fopen("t_12.txt","w");
    FILE* File13=fopen("t_13.txt","w");
    FILE* File14=fopen("t_14.txt","w");
    FILE* File15=fopen("t_15.txt","w");
    FILE* File16=fopen("t_16.txt","w");

    // // Расчет линейной задачи с явной схемой
    
    //     // Таблица 1
    //     calculate_linear_explicit_write_file(file1,file2,10,20);
    //     calculate_linear_explicit_write_file(file1,file2,100,20);
    //     calculate_linear_explicit_write_file(file1,file2,1000,20);
    //     calculate_linear_explicit_write_file(file1,file2,10,200);
    //     calculate_linear_explicit_write_file(file1,file2,100,200);
    //     calculate_linear_explicit_write_file(file1,file2,1000,200);
    //     calculate_linear_explicit_write_file(file1,file2,10,2000);
    //     calculate_linear_explicit_write_file(file1,file2,100,2000);
    //     calculate_linear_explicit_write_file(file1,file2,1000,2000);
    //     // Таблица 2
    //     calculate_linear_explicit_write_file_local(file3,10,20,4);
    //     // Таблица 3
    //     calculate_linear_explicit_write_file_local(file4,100,200,4);
    
    
    // // Расчет линейной задачи с неявной схемой
    
    //     // Таблица 1
    //     calculate_linear_implicit_write_file(file5,file6,10,20);
    //     calculate_linear_implicit_write_file(file5,file6,100,20);
    //     calculate_linear_implicit_write_file(file5,file6,1000,20);
    //     calculate_linear_implicit_write_file(file5,file6,10,200);
    //     calculate_linear_implicit_write_file(file5,file6,100,200);
    //     calculate_linear_implicit_write_file(file5,file6,1000,200);
    //     calculate_linear_implicit_write_file(file5,file6,10,2000);
    //     calculate_linear_implicit_write_file(file5,file6,100,2000);
    //     calculate_linear_implicit_write_file(file5,file6,1000,2000);
    //     // Таблица 2
    //     calculate_linear_implicit_write_file_local(file7,10,20,4);
    //     // Таблица 3
    //     calculate_linear_implicit_write_file_local(file8,100,200,4);


    // // Расчет нелинейной задачи с явной схемой

    //     // Таблица 1
    //     calculate_nonlinear_explicit_write_file(File1,File2,10,20);
    //     calculate_nonlinear_explicit_write_file(File1,File2,100,20);
    //     calculate_nonlinear_explicit_write_file(File1,File2,1000,20);
    //     calculate_nonlinear_explicit_write_file(File1,File2,10,200);
    //     calculate_nonlinear_explicit_write_file(File1,File2,100,200);
    //     calculate_nonlinear_explicit_write_file(File1,File2,1000,200);
    //     calculate_nonlinear_explicit_write_file(File1,File2,10,2000);
    //     calculate_nonlinear_explicit_write_file(File1,File2,100,2000);
    //     calculate_nonlinear_explicit_write_file(File1,File2,1000,2000);
    //     // Таблица 2
    //     calculate_nonlinear_explicit_write_file_local(File3,10,20,4);
    //     // Таблица 3
    //     calculate_nonlinear_explicit_write_file_local(File4,100,200,4);


    // // Расчет нелинейной задачи с неявной схемой

    //     // Таблица 1
    //     calculate_nonlinear_implicit_write_file(File5,File6,10,20);
    //     calculate_nonlinear_implicit_write_file(File5,File6,100,20);
    //     calculate_nonlinear_implicit_write_file(File5,File6,1000,20);
    //     calculate_nonlinear_implicit_write_file(File5,File6,10,200);
    //     calculate_nonlinear_implicit_write_file(File5,File6,100,200);
    //     calculate_nonlinear_implicit_write_file(File5,File6,1000,200);
    //     calculate_nonlinear_implicit_write_file(File5,File6,10,2000);
    //     calculate_nonlinear_implicit_write_file(File5,File6,100,2000);
    //     calculate_nonlinear_implicit_write_file(File5,File6,1000,2000);
    //     // Таблица 2
    //     calculate_nonlinear_implicit_write_file_local(File7,10,20,4);
    //     // Таблица 3
    //     calculate_nonlinear_implicit_write_file_local(File8,100,200,4);

    // // Выводим таблицу для нелинейной задачи с явной схемой и непрерывными Н.У.

    //     {
    //         double p=10.;
    //         int Nh=200;
    //         int Nt=100;
    //         for(int j=0;j<5;j++){
    //         calculate_nonlinear_explicit_write_file_tetha(File9,File10,Nt,Nh,1./p);
    //         Nh=Nh*2;
    //         Nt=Nt*2;
    //         p=p*2;        
    //         }
    //     }

    // // Выводим таблицу для нелинейной задачи с неявной схемой и непрерывными Н.У.

    //     {
    //         double p=5.;
    //         int Nh=200;
    //         int Nt=100;
    //         for(int j=0;j<5;j++){
    //         calculate_nonlinear_implicit_write_file_tetha(File11,File12,Nt,Nh,1./p);
    //         Nh=Nh*2;
    //         Nt=Nt*2;
    //         p=p*2;        
    //         }
    //     }

    double p=100.;
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,100,200,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,1000,200,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,10000,200,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,100,2000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,1000,2000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,10000,2000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,100,20000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,1000,20000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,10000,20000,1./p);

    // p=1000.;
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,100,200,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,1000,200,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,10000,200,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,100,2000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,1000,2000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,10000,2000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,100,20000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,1000,20000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,10000,20000,1./p);

    // p=10000.;
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,100,200,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,1000,200,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,10000,200,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,100,2000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,1000,2000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,10000,2000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,100,20000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,1000,20000,1./p);
    // calculate_nonlinear_explicit_write_file_tetha(File13,File14,10000,20000,1./p);

    p=100.;
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,100,200,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,1000,200,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,10000,200,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,100,2000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,1000,2000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,10000,2000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,100,20000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,1000,20000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,10000,20000,1./p);

    p=1000.;
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,100,200,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,1000,200,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,10000,200,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,100,2000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,1000,2000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,10000,2000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,100,20000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,1000,20000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,10000,20000,1./p);

    p=10000.;
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,100,200,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,1000,200,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,10000,200,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,100,2000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,1000,2000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,10000,2000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,100,20000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,1000,20000,1./p);
    calculate_nonlinear_implicit_write_file_tetha(File15,File16,10000,20000,1./p);

return 0;
}
