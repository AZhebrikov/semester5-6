#include<stdio.h>
#include<math.h>

int def1(FILE* file);
int def2(FILE* file);
int def3(FILE* file);
int def4(FILE* file);
int def5(FILE* file);
int def6(FILE* file);
int def7(FILE* file,double tol,double x_start,double y_start);
int def7_1(FILE* file,double tol,double x_start,double y_start,double T);
int def8(FILE* file,double tol, long int N,double x_start,double y_start);
int def9(FILE* file,double tol,double Tol,double TOL,double x_start,double y_start);
double fabs(double x);


int main(void)
{
    // FILE* f1=fopen("DEF_1.txt","w");
    // FILE* f2=fopen("DEF_2.txt","w");
    FILE* f3=fopen("DEF_3.txt","w");
    // FILE* f4=fopen("DEF_4.txt","w");
    // FILE* f5=fopen("DEF_5.txt","w");
    // FILE* f6=fopen("DEF_6.txt","w");
    // FILE* f7=fopen("DEF_7_9.txt","w");
    // FILE* f8=fopen("DEF_7_11.txt","w");
    // FILE* f9=fopen("DEF_7_13.txt","w");
    // FILE* f10=fopen("DEF_7_1.txt","w");
    // FILE* f11=fopen("DEF_8.txt","w");
    // FILE* f12=fopen("DEF_9.txt","w");

    // def1(f1);
    // def2(f2);
    def3(f3);
    // def4(f4);
    // def5(f5);
    // def6(f6);
    // def7(f7,pow(10.,-9),0.,1.);
    // def7(f8,pow(10.,-11),0.,1.);
    // def7(f9,pow(10.,-13),0.,1.);
    // def7_1(f10,pow(10.,-9),0.,0.,10.);
    // def8(f11,pow(10.,-8),4000,0.,100.);
    // def9(f12,pow(10.,-11),pow(10.,-11),pow(10.,-9),0.,1000.);

return 0;
}

int def1(FILE* file)
{
    {
        double epsilon=1.;
        while(1.+epsilon>1.)
            {
                epsilon=epsilon/2.;
            }
        fprintf(file,"%e\n",epsilon);
    }

    {
        float epsilon=1;
        while(1+epsilon>1)
            {
                epsilon=epsilon/2;
            }
        fprintf(file,"%e\n",epsilon);
    }

    {
        double A=1.;
        double Max=1.;
        double Min=1.;

        while(A+1>A)
        {
            A=A*2;
        }

        Max=A;
        A=A/2;
        Min=A;
        
        while(Max-Min>2)
        {
            if( 
                (A+((Max-Min)/2.))+1>(A+((Max-Min)/2.)) 
            ){
            A=A+(Max-Min)/2.;
            Min=A;
            }
            else{
                Max=Min+(Max-Min)/2.;
            }
        }
        while(A+1>A)
        {
            A=A+1;
        }
        fprintf(file,"%e\n",A);
    }

    {
        double A=1.;
        double Max=1.;
        double Min=1.;

        while(A+pow(10,20)>A)
        {
            A=A*2;
        }
        Max=A;
        A=A/2;
        Min=A;
        while(Max-Min>2*pow(10,20))
        {
            if(
                (A+((Max-Min)/2))+pow(10,20)>(A+((Max-Min)/2))
            ){
                A=A+(Max-Min)/2;
                Min=A;
            }
            else{
                Max=Min+(Max-Min)/2;
            }

        }
        while(A+pow(10,20)>A)
        {
            A=A+pow(10,20);
        }
        fprintf(file,"%e\n",A);
    }

return 0;
}

int def2(FILE* file)
{
    //First way:
    double I=log(7./6.);
    for(int i=1;i<32;i++)I=1./i-6*I;
    fprintf(file,"%e\n",I);

    //Second way:
    double integral_1000=0.;
    for(int i=0;i<1000;i++)
    {
        integral_1000+=1./1000*(pow(i/1000.+1./1000,31)/(6.+1./1000+i/1000.));
    }
    double integral_10000=0.;
    for(int i=0;i<10000;i++)
    {
        integral_10000+=1./10000*(pow(i/10000.+1./10000,31)/(6.+1./10000+i/10000.));
    }
    fprintf(file,"%e\n",integral_1000);
    fprintf(file,"%e\n",integral_10000);

    //Third way:
    I=0.;
    for(int i=59;i>30;i--)I=1./(6*(i+1))-I/6.;
    fprintf(file,"%e\n",I);

return 0;
}

int def3(FILE* file)
{
    //f(x)=5 sin(x)+x^2-e^x
    //f'(x)=5cos(x)+2x-e^x
    
    double R;
    double h=1./110000000000000;
    double x_0=3;
    for(int j=0;j<51;j++)
    {
        
        R=5*cos(x_0)+2*x_0-exp(x_0)-(  (5*sin(x_0+h)+pow(x_0+h,2)-exp(x_0+h))-(5*sin(x_0)+pow(x_0,2)-exp(x_0))  )/(h);
        fprintf(file,"%e %e\n",h,fabs(R));
        h=h/1.06;
    }

return 0;
}


double f_def4(double x){return exp(x)*exp(-8)+cos(x)*10+2*x/10;}
double F_def4(double x){return (exp(x)-1)*exp(-8)+sin(x)*10+pow(x,2)/10;}

double f1_def4(double x){return x;}
double F1_def4(double x){return pow(x,2)/2.;}

int def4(FILE *file)
{
    //f(t)=e^t/e^8+cox(t)*10+2t/10
    //F(t)=(e^t-1)/e^8+sin(t)*10+t^2/10

    double y,t,t_end;
    double step[5];
    step[0]=1.;for(int i=1;i<5;i++){step[i]=step[i-1]/10;}

    double step_const[5];
    step_const[0]=1.;for(int i=1;i<3;i++){step_const[i]=step_const[i-1]/10;}

    t_end=10.;

    fprintf(file,"END\n");
    for(int i=0;i<3;i++)
    {
        y=0.;
        t=0.;
        fprintf(file,"%e %e\n",t,y);
        while(t<t_end)
        {
            if(t+step[i]>t_end)
            {
                step[i]=t_end-t;
            }
            y=y+step[i]*f_def4(t);
            t=t+step[i];
            fprintf(file,"%e %e\n",t,y);
        }
        fprintf(file,"---\n%e %e\nEND\n",step_const[i],y-F_def4(10.));
    }

    fprintf(file,"ENDEND\n");

    step[0]=1.;for(int i=1;i<3;i++){step[i]=step[i-1]/10;}

    for(int i=0;i<3;i++)
    {
        y=0.;
        t=0.;
        fprintf(file,"%e %e\n",t,y);
        while(t<t_end)
        {
            if(t+step[i]>t_end){step[i]=t_end-t;}
            y=y+step[i]*f1_def4(t);
            t=t+step[i];
            fprintf(file,"%e %e\n",t,y);
        }

        fprintf(file,"---\n%e %e\nEND\n",step_const[i],y-F1_def4(10.));
    }
return 0;
}

int def5(FILE *file)
{
    double step[4];
    step[0]=0.1;for(int i=1;i<4;i++)step[i]=step[i-1]/10;
    
    double step_const[4];
    step_const[0]=0.1;for(int i=1;i<4;i++)step_const[i]=step_const[i-1]/10;

    double T[7];
    T[0]=M_PI;for(int i=1;i<7;i++)T[i]=T[i-1]*10;

    double x,x0;
    double y,y0;
    double t;

    double T_t[7][3];

    for(int j=0;j<4;j++)
    {
        fprintf(file,"END\n");
        x0=0.;
        y0=1.;
        x=0.;
        y=1.;
        t=0.;
        for(int i=0;i<7;i++)
        {step[j]=step_const[j];
            while(t<T[i])
            {
                if(t+step[j]>T[i])
                {
                    step[j]=T[i]-t;
                }
                x=x0+step[j]*y0;
                y=y0-step[j]*x0;
                t=t+step[j];
                x0=x;
                y0=y;
            }
            T_t[i][0]=T[i];
            T_t[i][1]=x-sin(t);
            T_t[i][2]=y-cos(t);
        }
        fprintf(file,"---\n");
        for(int k=0;k<7;k++)
        {
            fprintf(file,"%e %e %e\n",T_t[k][0],T_t[k][1],T_t[k][2]);
        }
    }

return 0;
}


double fx(double t,double x,double y){return y;}
double fy(double t,double x,double y){return -x;}

int def6(FILE* file)
{

    double t,x0,x,y0,y;

    double step[4];
    step[0]=0.1;for(int i=1;i<4;i++)step[i]=step[i-1]/10;

    double step_const[4];
    step_const[0]=0.1;for(int i=1;i<4;i++)step_const[i]=step_const[i-1]/10;

    double T[7];T[0]=M_PI;for(int i=1;i<7;i++)T[i]=T[i-1]*10;

    double T_t[7][3];

    double kx1,kx2,kx3,kx4,kx5,kx6,kx7;
    double ky1,ky2,ky3,ky4,ky5,ky6,ky7;

    for(int l=0;l<1;l++)
    {
        fprintf(file,"END\n");
        x0=0.;
        y0=1.;
        x=0.;
        y=1.;
        t=0.;
        
        for(int j=0;j<7;j++)
        {step[l]=step_const[l];
            while(t<T[j])
            {
                if(t+step[l]>T[j])step[l]=T[j]-t;

                kx1=fx(t,x0,y0);
                ky1=fy(t,x0,y0);

                kx2=fx(t+(1./5)*step[l],    x0+step[l]*((1./5)*kx1),    y0+step[l]*((1./5)*ky1) );
                ky2=fy(t+(1./5)*step[l],    x0+step[l]*((1./5)*kx1),    y0+step[l]*((1./5)*ky1) );

                kx3=fx(t+(3./10)*step[l],   x0+step[l]*((3./40)*kx1+(9./40)*kx2),   y0+step[l]*((3./40)*ky1+(9./40)*ky2) );
                ky3=fy(t+(3./10)*step[l],   x0+step[l]*((3./40)*kx1+(9./40)*kx2),   y0+step[l]*((3./40)*ky1+(9./40)*ky2) );

                kx4=fx(t+(4./5)*step[l],    x0+step[l]*((44./45)*kx1+(-56./15)*kx2+(32./9)*kx3),   y0+step[l]*((44./45)*ky1+(-56./15)*ky2+(32./9)*ky3) );
                ky4=fy(t+(4./5)*step[l],    x0+step[l]*((44./45)*kx1+(-56./15)*kx2+(32./9)*kx3),   y0+step[l]*((44./45)*ky1+(-56./15)*ky2+(32./9)*ky3) );

                kx5=fx(t+(8./9)*step[l],
                    x0+step[l]*((19372./6561)*kx1+(-25360./2187)*kx2+(64448./6561)*kx3+(-212./729)*kx4),
                    y0+step[l]*((19372./6561)*ky1+(-25360./2187)*ky2+(64448./6561)*ky3+(-212./729)*ky4) 
                    );
                ky5=fy(t+(8./9)*step[l],
                    x0+step[l]*((19372./6561)*kx1+(-25360./2187)*kx2+(64448./6561)*kx3+(-212./729)*kx4),
                    y0+step[l]*((19372./6561)*ky1+(-25360./2187)*ky2+(64448./6561)*ky3+(-212./729)*ky4) 
                    );

                kx6=fx(t+step[l],
                    x0+step[l]*((9017./3168)*kx1+(-355./33)*kx2+(46732./5247)*kx3+(49./176)*kx4+(-5103./18656)*kx5),
                    y0+step[l]*((9017./3168)*ky1+(-355./33)*ky2+(46732./5247)*ky3+(49./176)*ky4+(-5103./18656)*ky5) 
                    );
                ky6=fy(t+step[l],
                    x0+step[l]*((9017./3168)*kx1+(-355./33)*kx2+(46732./5247)*kx3+(49./176)*kx4+(-5103./18656)*kx5),
                    y0+step[l]*((9017./3168)*ky1+(-355./33)*ky2+(46732./5247)*ky3+(49./176)*ky4+(-5103./18656)*ky5) 
                    );

                kx7=fx(t+step[l],
                    x0+step[l]*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6),
                    y0+step[l]*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6) 
                    );
                ky7=fy(t+step[l],
                    x0+step[l]*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6),
                    y0+step[l]*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6) 
                    );

                x=x0+step[l]*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6);
                y=y0+step[l]*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6);

                x0=x;
                y0=y;
                t=t+step[l];
            }
            T_t[j][0]=T[j];
            T_t[j][1]=x-sin(t);
            T_t[j][2]=y-cos(t);
        }
        fprintf(file,"---\n");
        for(int k=0;k<7;k++)
        {
            fprintf(file,"%e %e %e %e\n",step_const[l],T_t[k][0],T_t[k][1],T_t[k][2]);
        }
    }
return 0;
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

double fabs(double x)
{
    if(x>0)return x;
    else return -x;
}


int def7(FILE* file,double tol,double x_start,double y_start)
{
    double t;
    double x0,x,x_d;
    double y0,y,y_d;

    double T[7];T[0]=M_PI;
    for(int i=1;i<7;i++){T[i]=T[i-1]*10;}

    double step=1.;

    double kx1,kx2,kx3,kx4,kx5,kx6,kx7;
    double ky1,ky2,ky3,ky4,ky5,ky6,ky7;

    double T_t[7][3];

    double facmin,facmax,fac;
    int p=5;
    fac=0.9;
    facmax=1.3;
    facmin=0.7;

    double errx,erry,err;

    x0=x_start;
    y0=y_start;
    x=x0;
    y=y0;
    t=0.;    

    for(int j=0;j<7;j++)
    {
        while(t<T[j])
        {
            if(t+step>T[j]){step=T[j]-t;}

            kx1=fx(t,x0,y0);
            ky1=fy(t,x0,y0);

            kx2=fx(t+(1./5)*step,    x0+step*((1./5)*kx1),    y0+step*((1./5)*ky1) );
            ky2=fy(t+(1./5)*step,    x0+step*((1./5)*kx1),    y0+step*((1./5)*ky1) );

            kx3=fx(t+(3./10)*step,   x0+step*((3./40)*kx1+(9./40)*kx2),   y0+step*((3./40)*ky1+(9./40)*ky2) );
            ky3=fy(t+(3./10)*step,   x0+step*((3./40)*kx1+(9./40)*kx2),   y0+step*((3./40)*ky1+(9./40)*ky2) );

            kx4=fx(t+(4./5)*step,    x0+step*((44./45)*kx1+(-56./15)*kx2+(32./9)*kx3),   y0+step*((44./45)*ky1+(-56./15)*ky2+(32./9)*ky3) );
            ky4=fy(t+(4./5)*step,    x0+step*((44./45)*kx1+(-56./15)*kx2+(32./9)*kx3),   y0+step*((44./45)*ky1+(-56./15)*ky2+(32./9)*ky3) );

            kx5=fx(t+(8./9)*step,
                x0+step*((19372./6561)*kx1+(-25360./2187)*kx2+(64448./6561)*kx3+(-212./729)*kx4),
                y0+step*((19372./6561)*ky1+(-25360./2187)*ky2+(64448./6561)*ky3+(-212./729)*ky4) 
                );
            ky5=fy(t+(8./9)*step,
                x0+step*((19372./6561)*kx1+(-25360./2187)*kx2+(64448./6561)*kx3+(-212./729)*kx4),
                y0+step*((19372./6561)*ky1+(-25360./2187)*ky2+(64448./6561)*ky3+(-212./729)*ky4) 
                );

            kx6=fx(t+step,
                x0+step*((9017./3168)*kx1+(-355./33)*kx2+(46732./5247)*kx3+(49./176)*kx4+(-5103./18656)*kx5),
                y0+step*((9017./3168)*ky1+(-355./33)*ky2+(46732./5247)*ky3+(49./176)*ky4+(-5103./18656)*ky5) 
                );
            ky6=fy(t+step,
                x0+step*((9017./3168)*kx1+(-355./33)*kx2+(46732./5247)*kx3+(49./176)*kx4+(-5103./18656)*kx5),
                y0+step*((9017./3168)*ky1+(-355./33)*ky2+(46732./5247)*ky3+(49./176)*ky4+(-5103./18656)*ky5) 
                );

            kx7=fx(t+step,
                x0+step*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6),
                y0+step*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6) 
                );
            ky7=fy(t+step,
                x0+step*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6),
                y0+step*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6) 
                );

            x=x0+step*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6);
            y=y0+step*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6);

            x_d=x0+step*((5179./57600)*kx1+(7571./16695)*kx3+(393./640)*kx4+(-92097./339200)*kx5+(187./2100)*kx6+(1./40)*kx7 );
            y_d=y0+step*((5179./57600)*ky1+(7571./16695)*ky3+(393./640)*ky4+(-92097./339200)*ky5+(187./2100)*ky6+(1./40)*ky7 );

            errx=fabs(x_d-x);
            erry=fabs(y_d-y);

            err=Max(errx,erry);

            if(err<tol)
            {
                x0=x;
                y0=y;
                t=t+step;
            }
            step=step*Min( facmax , Max(  facmin  ,  fac*powl(tol/err,1./(p+1))  ) );
        }
        T_t[j][0]=T[j];
        T_t[j][1]=x-sin(t);
        T_t[j][2]=y-cos(t);
    }
    fprintf(file,"---\n");
    for(int k=0;k<7;k++)
    {
        fprintf(file,"%le %le %le %le\n",tol,T_t[k][0],T_t[k][1],T_t[k][2]);
    }

return 0;
}

int def7_1(FILE* file,double tol,double x_start,double y_start,double T)
{
    double t;
    double x0,x,x_d;
    double y0,y,y_d;


    double step=1.;

    double kx1,kx2,kx3,kx4,kx5,kx6,kx7;
    double ky1,ky2,ky3,ky4,ky5,ky6,ky7;


    double facmin,facmax,fac;
    int p=5;
    fac=0.9;
    facmax=1.3;
    facmin=0.7;

    double errx,erry,err;

    x0=x_start;
    y0=y_start;
    x=x0;
    y=y0;
    t=0.;
    fprintf(file,"%e %e %e\n",t,x,y);
        while(t<T)
        {
            if(t+step>T){step=T-t;}

            kx1=fx(t,x0,y0);
            ky1=fy(t,x0,y0);

            kx2=fx(t+(1./5)*step,    x0+step*((1./5)*kx1),    y0+step*((1./5)*ky1) );
            ky2=fy(t+(1./5)*step,    x0+step*((1./5)*kx1),    y0+step*((1./5)*ky1) );

            kx3=fx(t+(3./10)*step,   x0+step*((3./40)*kx1+(9./40)*kx2),   y0+step*((3./40)*ky1+(9./40)*ky2) );
            ky3=fy(t+(3./10)*step,   x0+step*((3./40)*kx1+(9./40)*kx2),   y0+step*((3./40)*ky1+(9./40)*ky2) );

            kx4=fx(t+(4./5)*step,    x0+step*((44./45)*kx1+(-56./15)*kx2+(32./9)*kx3),   y0+step*((44./45)*ky1+(-56./15)*ky2+(32./9)*ky3) );
            ky4=fy(t+(4./5)*step,    x0+step*((44./45)*kx1+(-56./15)*kx2+(32./9)*kx3),   y0+step*((44./45)*ky1+(-56./15)*ky2+(32./9)*ky3) );

            kx5=fx(t+(8./9)*step,
                x0+step*((19372./6561)*kx1+(-25360./2187)*kx2+(64448./6561)*kx3+(-212./729)*kx4),
                y0+step*((19372./6561)*ky1+(-25360./2187)*ky2+(64448./6561)*ky3+(-212./729)*ky4) 
                );
            ky5=fy(t+(8./9)*step,
                x0+step*((19372./6561)*kx1+(-25360./2187)*kx2+(64448./6561)*kx3+(-212./729)*kx4),
                y0+step*((19372./6561)*ky1+(-25360./2187)*ky2+(64448./6561)*ky3+(-212./729)*ky4) 
                );

            kx6=fx(t+step,
                x0+step*((9017./3168)*kx1+(-355./33)*kx2+(46732./5247)*kx3+(49./176)*kx4+(-5103./18656)*kx5),
                y0+step*((9017./3168)*ky1+(-355./33)*ky2+(46732./5247)*ky3+(49./176)*ky4+(-5103./18656)*ky5) 
                );
            ky6=fy(t+step,
                x0+step*((9017./3168)*kx1+(-355./33)*kx2+(46732./5247)*kx3+(49./176)*kx4+(-5103./18656)*kx5),
                y0+step*((9017./3168)*ky1+(-355./33)*ky2+(46732./5247)*ky3+(49./176)*ky4+(-5103./18656)*ky5) 
                );

            kx7=fx(t+step,
                x0+step*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6),
                y0+step*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6) 
                );
            ky7=fy(t+step,
                x0+step*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6),
                y0+step*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6) 
                );

            x=x0+step*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6);
            y=y0+step*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6);

            x_d=x0+step*((5179./57600)*kx1+(7571./16695)*kx3+(393./640)*kx4+(-92097./339200)*kx5+(187./2100)*kx6+(1./40)*kx7 );
            y_d=y0+step*((5179./57600)*ky1+(7571./16695)*ky3+(393./640)*ky4+(-92097./339200)*ky5+(187./2100)*ky6+(1./40)*ky7 );

            errx=fabs(x_d-x);
            erry=fabs(y_d-y);

            err=Max(errx,erry);

            if(err<tol)
            {
                x0=x;
                y0=y;
                t=t+step;
                fprintf(file,"%e %e %e\n",t,x,y);
            }
            step=step*Min( facmax , Max(  facmin  ,  fac*powl(tol/err,1./(p+1))  ) );
        }

return 0;
}


double Fx(double t,double x,double y)
{
    return y;
}
double Fy(double t,double x,double y)
{
    return 0.1*(1-x*x)*y-x;
}



int def8(FILE* file,double tol,long int N,double x_start,double y_start)
{
    fprintf(file,"END\n");

    double t;
    double x0,x,x_d;
    double y0,y,y_d;

    long int N_tmp;

    double step=1.;

    double kx1,kx2,kx3,kx4,kx5,kx6,kx7;
    double ky1,ky2,ky3,ky4,ky5,ky6,ky7;

    double facmin,facmax,fac;
    int p=5;
    fac=0.9;
    facmax=1.3;
    facmin=0.7;

    double errx,erry,err;

    x0=x_start;  x=x0;
    y0=y_start;  y=y0;
    t=0.;   N_tmp=0;

        while(N_tmp<N)
        {
            kx1=Fx(t,x0,y0);
            ky1=Fy(t,x0,y0);

            kx2=Fx(t+(1./5)*step,    x0+step*((1./5)*kx1),    y0+step*((1./5)*ky1) );
            ky2=Fy(t+(1./5)*step,    x0+step*((1./5)*kx1),    y0+step*((1./5)*ky1) );

            kx3=Fx(t+(3./10)*step,   x0+step*((3./40)*kx1+(9./40)*kx2),   y0+step*((3./40)*ky1+(9./40)*ky2) );
            ky3=Fy(t+(3./10)*step,   x0+step*((3./40)*kx1+(9./40)*kx2),   y0+step*((3./40)*ky1+(9./40)*ky2) );

            kx4=Fx(t+(4./5)*step,    x0+step*((44./45)*kx1+(-56./15)*kx2+(32./9)*kx3),   y0+step*((44./45)*ky1+(-56./15)*ky2+(32./9)*ky3) );
            ky4=Fy(t+(4./5)*step,    x0+step*((44./45)*kx1+(-56./15)*kx2+(32./9)*kx3),   y0+step*((44./45)*ky1+(-56./15)*ky2+(32./9)*ky3) );

            kx5=Fx(t+(8./9)*step,
                x0+step*((19372./6561)*kx1+(-25360./2187)*kx2+(64448./6561)*kx3+(-212./729)*kx4),
                y0+step*((19372./6561)*ky1+(-25360./2187)*ky2+(64448./6561)*ky3+(-212./729)*ky4) 
                );
            ky5=Fy(t+(8./9)*step,
                x0+step*((19372./6561)*kx1+(-25360./2187)*kx2+(64448./6561)*kx3+(-212./729)*kx4),
                y0+step*((19372./6561)*ky1+(-25360./2187)*ky2+(64448./6561)*ky3+(-212./729)*ky4) 
                );

            kx6=Fx(t+step,
                x0+step*((9017./3168)*kx1+(-355./33)*kx2+(46732./5247)*kx3+(49./176)*kx4+(-5103./18656)*kx5),
                y0+step*((9017./3168)*ky1+(-355./33)*ky2+(46732./5247)*ky3+(49./176)*ky4+(-5103./18656)*ky5) 
                );
            ky6=Fy(t+step,
                x0+step*((9017./3168)*kx1+(-355./33)*kx2+(46732./5247)*kx3+(49./176)*kx4+(-5103./18656)*kx5),
                y0+step*((9017./3168)*ky1+(-355./33)*ky2+(46732./5247)*ky3+(49./176)*ky4+(-5103./18656)*ky5) 
                );

            kx7=Fx(t+step,
                x0+step*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6),
                y0+step*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6) 
                );
            ky7=Fy(t+step,
                x0+step*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6),
                y0+step*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6) 
                );

            x=x0+step*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6);
            y=y0+step*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6);

            x_d=x0+step*((5179./57600)*kx1+(7571./16695)*kx3+(393./640)*kx4+(-92097./339200)*kx5+(187./2100)*kx6+(1./40)*kx7 );
            y_d=y0+step*((5179./57600)*ky1+(7571./16695)*ky3+(393./640)*ky4+(-92097./339200)*ky5+(187./2100)*ky6+(1./40)*ky7 );

            errx=fabs(x_d-x);
            erry=fabs(y_d-y);

            err=Max(errx,erry);

            if(err<tol)
            {
                x0=x;
                y0=y;
                t=t+step;
                fprintf(file,"%le %le\n",x,y);
                N_tmp++;
            }

            step=step*Min( facmax , Max(  facmin  ,  fac*powl(tol/err,1./(p+1))  ) );
        }

    fprintf(file,"END\n");

return 0;
}


int def9(FILE* file,double tol,double Tol,double TOL,double x_start,double y_start)
{
    double t,T_tmp,T;
    double x0,x,x_d,X_tmp,X,X0;
    double y0,y,y_d,Y_tmp,Y,Y0;

    double step=1.;
    double Step,step_tmp;
    
    double kx1,kx2,kx3,kx4,kx5,kx6,kx7;
    double ky1,ky2,ky3,ky4,ky5,ky6,ky7;

    double facmin,facmax,fac;
    int p=5;
    fac=0.9;
    facmax=1.3;
    facmin=0.7;

    int flag;
    int flagA;
    int flagB;

    double errx,erry,err;

    x0=x_start;  x=x0;
    y0=y_start;  y=y0;
    t=0.;   
    
    flag=1;
    flagA=0;
    flagB=10;

    double x1,x2;

    x1=y0;
    x2=y0;

    do{
                    x1=x2;
                    flagB=0;

                    while(flag==1)
                    {
                        kx1=Fx(t,x0,y0);
                        ky1=Fy(t,x0,y0);

                        kx2=Fx(t+(1./5)*step,    x0+step*((1./5)*kx1),    y0+step*((1./5)*ky1) );
                        ky2=Fy(t+(1./5)*step,    x0+step*((1./5)*kx1),    y0+step*((1./5)*ky1) );

                        kx3=Fx(t+(3./10)*step,   x0+step*((3./40)*kx1+(9./40)*kx2),   y0+step*((3./40)*ky1+(9./40)*ky2) );
                        ky3=Fy(t+(3./10)*step,   x0+step*((3./40)*kx1+(9./40)*kx2),   y0+step*((3./40)*ky1+(9./40)*ky2) );

                        kx4=Fx(t+(4./5)*step,    x0+step*((44./45)*kx1+(-56./15)*kx2+(32./9)*kx3),   y0+step*((44./45)*ky1+(-56./15)*ky2+(32./9)*ky3) );
                        ky4=Fy(t+(4./5)*step,    x0+step*((44./45)*kx1+(-56./15)*kx2+(32./9)*kx3),   y0+step*((44./45)*ky1+(-56./15)*ky2+(32./9)*ky3) );

                        kx5=Fx(t+(8./9)*step,
                            x0+step*((19372./6561)*kx1+(-25360./2187)*kx2+(64448./6561)*kx3+(-212./729)*kx4),
                            y0+step*((19372./6561)*ky1+(-25360./2187)*ky2+(64448./6561)*ky3+(-212./729)*ky4) 
                            );
                        ky5=Fy(t+(8./9)*step,
                            x0+step*((19372./6561)*kx1+(-25360./2187)*kx2+(64448./6561)*kx3+(-212./729)*kx4),
                            y0+step*((19372./6561)*ky1+(-25360./2187)*ky2+(64448./6561)*ky3+(-212./729)*ky4) 
                            );

                        kx6=Fx(t+step,
                            x0+step*((9017./3168)*kx1+(-355./33)*kx2+(46732./5247)*kx3+(49./176)*kx4+(-5103./18656)*kx5),
                            y0+step*((9017./3168)*ky1+(-355./33)*ky2+(46732./5247)*ky3+(49./176)*ky4+(-5103./18656)*ky5) 
                            );
                        ky6=Fy(t+step,
                            x0+step*((9017./3168)*kx1+(-355./33)*kx2+(46732./5247)*kx3+(49./176)*kx4+(-5103./18656)*kx5),
                            y0+step*((9017./3168)*ky1+(-355./33)*ky2+(46732./5247)*ky3+(49./176)*ky4+(-5103./18656)*ky5) 
                            );

                        kx7=Fx(t+step,
                            x0+step*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6),
                            y0+step*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6) 
                            );
                        ky7=Fy(t+step,
                            x0+step*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6),
                            y0+step*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6) 
                            );

                        x=x0+step*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6);
                        y=y0+step*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6);

                        x_d=x0+step*((5179./57600)*kx1+(7571./16695)*kx3+(393./640)*kx4+(-92097./339200)*kx5+(187./2100)*kx6+(1./40)*kx7 );
                        y_d=y0+step*((5179./57600)*ky1+(7571./16695)*ky3+(393./640)*ky4+(-92097./339200)*ky5+(187./2100)*ky6+(1./40)*ky7 );

                        errx=fabs(x_d-x);
                        erry=fabs(y_d-y);

                        err=Max(errx,erry);

                        if(err<tol)
                        {
                            if(y>0. && y0>0. && x0<=0. && x>=0. && flagB>=10)
                            {
                                flag=0;
                                X=x;
                                X0=x0;
                                Step=step;
                                T=t;
                                Y0=y0;
                                Y=y;
                            }
                                if(flagA==1)fprintf(file,"%le %le %le\n",t,x0,y0);
                                x0=x;
                                y0=y;
                                t=t+step;
                                flagB++;
                        }
                        step=step*Min( facmax , Max(  facmin  ,  fac*powl(tol/err,1./(p+1))  ) );
                    }

// T , Step, Y, Y0, X, X0

                    if(X0<0 && X>0)
                    {
                                do{
                                    T_tmp=T+(fabs(X0)*Step)/(X+fabs(X0));
                                    step_tmp=T_tmp-T;

                                    kx1=Fx(T,X0,Y0);
                                    ky1=Fy(T,X0,Y0);

                                    kx2=Fx(T+(1./5)*step_tmp,    X0+step_tmp*((1./5)*kx1),    Y0+step_tmp*((1./5)*ky1) );
                                    ky2=Fy(T+(1./5)*step_tmp,    X0+step_tmp*((1./5)*kx1),    Y0+step_tmp*((1./5)*ky1) );

                                    kx3=Fx(T+(3./10)*step_tmp,   X0+step_tmp*((3./40)*kx1+(9./40)*kx2),   Y0+step_tmp*((3./40)*ky1+(9./40)*ky2) );
                                    ky3=Fy(T+(3./10)*step_tmp,   X0+step_tmp*((3./40)*kx1+(9./40)*kx2),   Y0+step_tmp*((3./40)*ky1+(9./40)*ky2) );

                                    kx4=Fx(T+(4./5)*step_tmp,    X0+step_tmp*((44./45)*kx1+(-56./15)*kx2+(32./9)*kx3),   Y0+step_tmp*((44./45)*ky1+(-56./15)*ky2+(32./9)*ky3) );
                                    ky4=Fy(T+(4./5)*step_tmp,    X0+step_tmp*((44./45)*kx1+(-56./15)*kx2+(32./9)*kx3),   Y0+step_tmp*((44./45)*ky1+(-56./15)*ky2+(32./9)*ky3) );

                                    kx5=Fx(T+(8./9)*step_tmp,
                                        X0+step_tmp*((19372./6561)*kx1+(-25360./2187)*kx2+(64448./6561)*kx3+(-212./729)*kx4),
                                        Y0+step_tmp*((19372./6561)*ky1+(-25360./2187)*ky2+(64448./6561)*ky3+(-212./729)*ky4) 
                                        );
                                    ky5=Fy(T+(8./9)*step_tmp,
                                        X0+step_tmp*((19372./6561)*kx1+(-25360./2187)*kx2+(64448./6561)*kx3+(-212./729)*kx4),
                                        Y0+step_tmp*((19372./6561)*ky1+(-25360./2187)*ky2+(64448./6561)*ky3+(-212./729)*ky4) 
                                        );

                                    kx6=Fx(T+step_tmp,
                                        X0+step_tmp*((9017./3168)*kx1+(-355./33)*kx2+(46732./5247)*kx3+(49./176)*kx4+(-5103./18656)*kx5),
                                        Y0+step_tmp*((9017./3168)*ky1+(-355./33)*ky2+(46732./5247)*ky3+(49./176)*ky4+(-5103./18656)*ky5) 
                                        );
                                    ky6=Fy(T+step_tmp,
                                        X0+step_tmp*((9017./3168)*kx1+(-355./33)*kx2+(46732./5247)*kx3+(49./176)*kx4+(-5103./18656)*kx5),
                                        Y0+step_tmp*((9017./3168)*ky1+(-355./33)*ky2+(46732./5247)*ky3+(49./176)*ky4+(-5103./18656)*ky5) 
                                        );

                                    kx7=Fx(T+step_tmp,
                                        X0+step_tmp*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6),
                                        Y0+step_tmp*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6) 
                                        );
                                    ky7=Fy(T+step_tmp,
                                        X0+step_tmp*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6),
                                        Y0+step_tmp*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6) 
                                        );

                                    X_tmp=X0+step_tmp*((35./384)*kx1+(500./1113)*kx3+(125./192)*kx4+(-2187./6784)*kx5+(11./84)*kx6);
                                    Y_tmp=Y0+step_tmp*((35./384)*ky1+(500./1113)*ky3+(125./192)*ky4+(-2187./6784)*ky5+(11./84)*ky6);

                                    if(X_tmp<0)
                                        {
                                            X0=X_tmp;
                                            Y0=Y_tmp;
                                            Step=T+Step-T_tmp;
                                            T=T_tmp;
                                        }
                                    else
                                        {
                                            X=X_tmp;
                                            Y=Y_tmp;
                                            Step=T+Step-T_tmp;
                                        }

                                }while(fabs(X_tmp)>Tol);
                                x2=Y_tmp;
                                flag=1;
                                }
                                else if (X0>=0)
                                {
                                    x2=Y0;
                                    flag=1;
                                }
                                else 
                                {
                                    x2=Y;
                                    flag=1;
                                }

    step=1.;
    if(flagA==1)
    {
    fprintf(file,"%le %le %le\n",T_tmp,X_tmp,Y_tmp);
    flagA=2;
    }

    if((fabs(x2-x1))<TOL && flagA!=2)
    {
        flagA=1;
        fprintf(file,"%e %e %e\n",T_tmp,X_tmp,Y_tmp);
    }

    }while(flagA!=2);

    fprintf(file,"---\n %le %le %le",tol,Tol,TOL);

return 0;
}
