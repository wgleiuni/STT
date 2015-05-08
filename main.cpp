/*
 * =====================================================================================
 *
 *       Filename:  main.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/16/2015 08:26:32 AM
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
#include "DefClass.h"

int main ( int argc, char *argv[] )
{
    int N=0,P=0,Len=0;
    double Start=0.0,End=0.0;
    if (argc>1)
    {
        N=(int)atof(argv[1]);
        P=(int)atof(argv[2]);
        Start=(double)atof(argv[3]);
        End=(double)atof(argv[4]);
        Len=(int)atof(argv[5]);
        std::cout << N << " " << P << " " << Start << " " << End << " " << Len << std::endl;
    }
    STT one(N,P,Start,End,Len);
    one.go();
//    one.disp1();
//    one.out(N);

    return 0;
}
