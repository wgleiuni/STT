/*
 * =====================================================================================
 *
 *       Filename:  DefClass.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/15/2015 10:34:16 AM
 *       Revision:  none
 *       Compiler:  icpc
 *
 *         Author:  Guanglei Wang (glwang), wgleiuni@gmail.com
 *   Organization:  Arizona State University
 *
 * =====================================================================================
 */

#ifndef DEFCLASS_H
#define DEFCLASS_H
#include <fstream>
class Para
{
    public:
        Para();
        void get(int,double);
        void disp();
    protected:
        int N_,P_,Len_;
        double n[3];
        double alpha_,omega_,h_,kf_,V_,ac_,dc_,dt_,numdt_,tau_,Start_,End_;
        double* mag[3];
};

class Sigma : protected Para
{
    public:
        Sigma();
    protected:
        double sigma1_,sigma2_,rsigma1_,rsigma2_;
        void integrate();
    private:
        double T2(double theta);
};

class Rk4 : protected Sigma
{
    public:
        Rk4();
        void onestep();
    private:
        double dt(int,double,double*,double*);
};

class STT : protected Rk4
{
    public:
        STT(int,int,double,double,int);
        void initial();
        void go();
        void out(int);
        void disp1();
    private:
        void record();
        std::ofstream outnx_,outny_,outnz_,outs1_,outs2_;
};
#endif
