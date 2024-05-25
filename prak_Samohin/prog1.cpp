        {
            double f;
            double p=10.;
            int Nh=200;
            int Nt=100;
            while(p<5000.)
            {
                f=calculate_nonlinear_explicit_write_file_tetha(File9,File10,Nt,Nh,1./p);
                if(f<1./10){
                    p=p*2;
                }
                else{
                Nh=Nh*2;
                Nt=Nt*2;
                }
            }
        }

        {
            double f;
            double p=10.;
            int Nh=200;
            int Nt=100;
            while(p<2000.)
            {
                f=calculate_nonlinear_implicit_write_file_tetha(File11,File12,Nt,Nh,1./p);
                if(f<1./10){
                    p=p*2;
                }
                else{
                Nh=Nh*2;
                Nt=Nt*2;
                }
            }
        }
