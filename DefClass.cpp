/*
 * =====================================================================================
 *
 *       Filename:  DefClass.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/15/2015 10:41:58 AM
 *       Revision:  none
 *       Compiler:  icpc
 *
 *         Author:  Guanglei Wang (glwang), wgleiuni@gmail.com
 *   Organization:  Arizona State University
 *
 * =====================================================================================
 */

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex>
#include <fstream>
#include "DefClass.h"

Para::Para()
{
    std::ifstream input;
    std::string str,temp;
    double x;
    int i=1,found;

    input.open("parameter.txt",std::ifstream::in);

    while (std::getline(input,str))
    {
        found=str.find("=");
        temp.assign(str,found+1,str.length()-found-1);
        x=atof(temp.c_str());
        get(i,x);
        i++;
    }

    input.close();
    tau_=0.0;
    mag[0]=new double[(int)numdt_];
    mag[1]=new double[(int)numdt_];
    mag[2]=new double[(int)numdt_];
}

void Para::get(int N,double x)
{
    switch (N) {
        case 1:
            alpha_=x; break;
        case 2:
            omega_=x; break;
        case 3:
            h_=x; break;
        case 4:
            kf_=x; break;
        case 5:
            V_=x; break;
        case 6:
            ac_=x; break;
        case 7:
            dc_=x; break;
        case 8:
            n[0]=x; break;
        case 9:
            n[1]=x; break;
        case 10:
            n[2]=x; break;
        case 11:
            dt_=x; break;
        case 12:
            numdt_=x; break;
        case 13:
            double theta,phi;

            theta=-M_PI/2.0+M_PI/Start_*(N_/(int)End_);
            phi=2*M_PI/End_*(N_%(int)End_);
            n[0]=cos(theta)*cos(phi);
            n[1]=cos(theta)*sin(phi);
            n[2]=sin(theta);
            break;
    }
}

void Para::disp()
{
    char filename[20];
    std::ofstream output;

    sprintf(filename,"parameter%d.in",N_);

    output.open(filename,std::ostream::out);
    output << "dc=" << dc_ << std::endl;
    output.close();
}

Sigma::Sigma()
{

}

double Sigma::T2(double theta)
{
    double m[3],out;
    std::complex<double> I (0.0,1.0);
    std::complex<double> kx,ky,a,a1,a2,A,B,t;

    m[0]=h_*kf_*n[0];
    m[1]=h_*kf_*n[1];
    m[2]=h_*kf_*n[2];

    ky=kf_*sin(theta)+m[0];
    kx=sqrt(kf_*kf_-m[2]*m[2]-ky*ky);

    a=exp(I*kf_*cos(theta));
    a1=exp(I*(kx+m[1]));
    a2=exp(I*(-kx+m[1]));

    A=a2*(I*exp(-I*theta)*(ky+I*kx)-kf_-m[2])-a1*(I*exp(-I*theta)*(ky-I*kx)-kf_-m[2]);
    B=a2*(I*exp(-I*theta)*(kf_-m[2])-(ky-I*kx))-a1*(I*exp(-I*theta)*(kf_-m[2])-(ky+I*kx));

    t=-4.0*cos(theta)*kx/(a*(A+I*exp(I*theta)*B));

    out=(double) std::norm(t);

    return out;
}

void Sigma::integrate()
{
    sigma1_=0.0;
    sigma2_=0.0;

    int N=10002;
    double dtheta=M_PI/N,theta,t2;

    int i=0;

    for (i=0;i<N+1;i++)
    {
        theta=-M_PI/2+i*dtheta;
        t2=T2(theta);

        if (i%3==0)
        {
            sigma1_=sigma1_+2.0*t2*sin(theta);
            sigma2_=sigma2_+2.0*t2*cos(theta);
        }
        else if (i%3==1 || i%3==2)
        {
            sigma1_=sigma1_+3.0*t2*sin(theta);
            sigma2_=sigma2_+3.0*t2*cos(theta);
        }
    }
    sigma1_=sigma1_-T2(-M_PI/2)*sin(-M_PI/2)-T2(M_PI/2)*sin(M_PI/2);
    sigma2_=sigma2_-T2(-M_PI/2)*cos(-M_PI/2)-T2(M_PI/2)*cos(M_PI/2);

    sigma1_=3.0/8.0*dtheta*sigma1_;
    sigma2_=-3.0/8.0*dtheta*sigma2_;
}

Rk4::Rk4()
{

}

void Rk4::onestep()
{
    double k[4],l[4],j[4],p[3],s[2],tau;
    int i;

    s[0]=sigma1_;
    s[1]=sigma2_;

    double b=h_*V_*(dc_+ac_*cos(omega_*tau_));
    rsigma1_=b*sigma1_;
    rsigma2_=b*sigma2_;

    for (i=1;i<5;i++)
    {
        switch (i) {
            case 1:
                tau=tau_;
                p[0]=n[0];  p[1]=n[1];  p[2]=n[2];  break;
            case 2:
                tau=tau_+dt_/2.0;
                p[0]=n[0]+dt_/2.0*k[i-1]; p[1]=n[1]+dt_/2.0*l[i-1]; p[2]=n[2]+dt_/2.0*j[i-1]; break;
            case 3:
                tau=tau_+dt_/2.0;
                p[0]=n[0]+dt_/2.0*k[i-1]; p[1]=n[1]+dt_/2.0*l[i-1]; p[2]=n[2]+dt_/2.0*j[i-1]; break;
            case 4:
                tau=tau_+dt_;
                p[0]=n[0]+dt_*k[i-1]; p[1]=n[1]+dt_*l[i-1]; p[2]=n[2]+dt_*j[i-1]; break;
        }
        k[i]=dt(1,tau,p,s);
        l[i]=dt(2,tau,p,s);
        j[i]=dt(3,tau,p,s);
    }

    n[0]=n[0]+dt_/6.0*(k[1]+2.0*k[2]+2.0*k[3]+k[4]);
    n[1]=n[1]+dt_/6.0*(l[1]+2.0*l[2]+2.0*l[3]+l[4]);
    n[2]=n[2]+dt_/6.0*(j[1]+2.0*j[2]+2.0*j[3]+j[4]);

    tau_=tau_+dt_;
}

double Rk4::dt(int N,double tau,double* p, double* s)
{
    double out;
    double a=alpha_;
    double b=h_*V_*(dc_+ac_*cos(omega_*tau));

    switch (N) {
        case 1:
//        std::cout<< p[0] << " " << p[1] << " " << p[2] << " " << s[0] << " " << s[2] << std::endl;
            out=-(a*p[1]*p[1]+a*p[2]*p[2]-b*p[2]*s[1]-a*b*p[1]*p[1]*s[0]-a*b*p[2]*p[2]*s[0]+a*b*p[0]*p[1]*s[1])
                /(a*a*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])+1.0);
//            std::cout<< a << " " << b << " " << out << std::endl;
            break;
        case 2:
            out=(p[2]+a*p[0]*p[1]-b*p[2]*s[0]+a*b*p[0]*p[0]*s[1]+a*b*p[2]*p[2]*s[1]-a*b*p[0]*p[1]*s[0])
                /(a*a*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])+1.0);
            break;
        case 3:
            out=(a*p[0]*p[2]-p[1]-b*p[0]*s[1]+b*p[1]*s[0]-a*b*p[0]*p[2]*s[0]-a*b*p[1]*p[2]*s[1])
                /(a*a*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])+1.0);
            break;
    }

    return out;
}

STT::STT(int N,int P,double Start,double End,int Len)
{
    N_=N;
    P_=P;
    Start_=Start;
    End_=End;
    Len_=Len;

    char filename1[20],filename2[20],filename3[20],filename4[20],filename5[20];

    sprintf(filename1,"nx%d.txt",N);
    sprintf(filename2,"ny%d.txt",N);
    sprintf(filename3,"nz%d.txt",N);
    sprintf(filename4,"sx%d.txt",N);
    sprintf(filename5,"sy%d.txt",N);

    outnx_.open(filename1,std::ostream::out);
    outny_.open(filename2,std::ostream::out);
    outnz_.open(filename3,std::ostream::out);
    outs1_.open(filename4,std::ostream::out);
    outs2_.open(filename5,std::ostream::out);
}

void STT::initial()
{
    if (N_==0) return;
    if (N_>0)
    {
        Rk4::Para::get(P_,Start_+(End_-Start_)/Len_*N_);
    }
}

void STT::go()
{
    int i;
    initial();
    for (i=0;i<(int)numdt_;i++)
    {
        if (i%10==0)
        {
            record(); 
        }
        Rk4::Sigma::integrate();
        Rk4::onestep();
        *(mag[0]+i)=n[0];
        *(mag[1]+i)=n[1];
        *(mag[2]+i)=n[2];
//        if (i%100==0) std::cout<< i << std::endl;
    }
}

void STT::record()
{
    outnx_ << n[0] << std::endl;
    outny_ << n[1] << std::endl;
    outnz_ << n[2] << std::endl;
    outs1_ << sigma1_ << "\t" << rsigma1_ << std::endl;
    outs2_ << sigma2_ << "\t" << rsigma2_ << std::endl;
}

void STT::out(int N)
{
    std::ofstream output;
    char filename[20];

    sprintf(filename,"ny%d.txt",N);
    output.open(filename,std::ostream::out);

    int i;
    for (i=0;i<(int)numdt_;i++)
    {
        output << *(mag[1]+i) << std::endl;
    }
    output.close();
}

void STT::disp1()
{
    Rk4::Sigma::Para::disp();
}
