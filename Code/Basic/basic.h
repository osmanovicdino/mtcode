/* 
 * File:   basic.h
 * Author: Dino
 *
 * Created on 28 November 2010, 00:14
 */


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdarg.h>
#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <sstream>
#include <dirent.h>

using namespace std;

/* This file contains many basic functions for error checking and output,
 * as well as some common mathemtical operations */

#ifndef BASIC_H
#define	BASIC_H



const long double pii = 3.1415926535897932384626433832795028841971693993751;
const long double pi = 3.1415926535897932384626433832795028841971693993751;
const long double ee=2.7182818284590452353;
const double eee=2.7182818284590452353;
const double pid = 3.1415926535897932384626433832795028841971693993751;

template <class T>
inline void output(char *s, T x ) { cout << s << x; }
template <class T>
inline void input (char *s, T &x) { cout << s; cin >> x; }


//definitions, mainly inline functions, and mathematica plotting tools
#define check cout << "fine up to here" << endl; //basic check
#define mathematica myfile << "\nll=" << //input into mathematica
#define convert(dr) myfile << "\nTable[ll[[i]]={(i-1)*" << dr << ",ll[[i]]},{i,1,Length[ll]}];" << endl;
#define listplot3 myfile << "\nr1=ListPlot[ll,PlotJoined->True,PlotRange->All]" << endl;
#define listplot(a) myfile << "\nr1=ListPlot[ll,PlotJoined->True,PlotRange->All,PlotLabel->"<<a<<"]" << endl;
#define listplot2(a,b) a << "\nr1=ListPlot["<<b<<",PlotJoined->True,PlotRange->All]" << endl;
#define listplot3d(q) q << "\nListPlot3D[ll,PlotRange->All]"<<endl;
#define G(r0,R,t,k) myfile << "\nr2=Plot[G[r," <<r0<<","<<R<<","<<t<<","<<k<<"], {r, 0, "<<R<<"}, PlotRange -> All]" << endl;
#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))
#define SIGN(a,b)((b)>= 0 ? fabs(a) : -fabs(a))
#define ABS(a)( (a>0? (a) :(-a)) )

int inline fasterfloor( const double &x ) { return x > 0 ? (int) x : (int) x - 1; }


inline void error(const char *errmsg)
{
	cerr << "\n Runtime error ... \n";
	cerr << errmsg;
	cerr << "\n...Quitting program! \n";


        double a;
        cin >> a; //pause the program
	exit(1);
}

template <class T>
void declare(ofstream &s, T a, char *ss, bool c = false ) {
    if (c) {
    s.open(ss,ios::app);
    }
    else {
        s.open(ss);
    }
    s.precision(20);
    s << "ll=" << fixed << a << endl;
    listplot3d(s);
}

template <class T>
void mathe(ofstream &s, T a) { //output a templated object T to a file s
    if ( !s.is_open() ) error("open a file first");
    s << "ll=" << fixed << a << endl;
    listplot3d(s);
}



inline void pausel() { //pause the program (this is faster than system("pause")
    cout << "\nenter an integer: " << endl;
    int pause;
    cin >> pause;
}
#define Stop pausel();


template <class T>
void SWAP( T &a, T &b ) { //swap the values of two objects a and b
	T c = a;
	a = b;
	b = c;
}

inline void sort_doublet(int ti1, int ti2, int &i1, int &i2) {
    if (ti1 < ti2)
    {
        i1 = ti1;
        i2 = ti2;
    }
    else
    {
        i1 = ti2;
        i2 = ti1;
    }
}

inline void sort_triplet(int t1, int t2, int t3, int &i1, int &i2, int &i3) {
    if(t1 < t2) {
        if(t2 < t3) {
            i1 = t1;
            i2 = t2;
            i3 = t3;
        }
        else{
            if(t1<t3) {
                i1 = t1;
                i2 = t3;
                i3 = t2;
            }
            else{
                i1 = t3;
                i2 = t1;
                i3 = t2;
            }
        }

    }
    else{
        if(t3 < t2) {
            i1 = t3;
            i2 = t2;
            i3 = t1;
        }
        else{
            if(t3 > t1) {
                i1 = t2;
                i2 = t1;
                i3 = t3;
            }
            else{
                i1 = t2;
                i2 = t3;
                i3 = t1;
            }
        }

    }
}

inline void sort_triplet_and_save_permutation(int t1, int t2, int t3, int &i1, int &i2, int &i3, unsigned char &p1, unsigned char &p2, unsigned char &p3)
{
    if (t1 < t2)
    {
        if (t2 < t3)
        {
            i1 = t1;
            i2 = t2;
            i3 = t3;

            p1 = 1;
            p2 = 2;
            p3 = 3;
        }
        else
        {
            if (t1 < t3)
            {
                i1 = t1;
                i2 = t3;
                i3 = t2;

                p1 = 1;
                p2 = 3;
                p3 = 2;
            }
            else
            {
                i1 = t3;
                i2 = t1;
                i3 = t2;

                p1 = 3;
                p2 = 1;
                p3 = 2;
            }
        }
    }
    else
    {
        if (t3 < t2)
        {
            i1 = t3;
            i2 = t2;
            i3 = t1;

            p1 = 3;
            p2 = 2;
            p3 = 1;
        }
        else
        {
            if (t3 > t1)
            {
                i1 = t2;
                i2 = t1;
                i3 = t3;

                p1 = 2;
                p2 = 1;
                p3 = 3;
            }
            else
            {
                i1 = t2;
                i2 = t3;
                i3 = t1;

                p1 = 2;
                p2 = 3;
                p3 = 1;
            }
        }
    }
}

inline void save_permutation_triple(const bool &b1, const bool &b2, const bool &b3, const unsigned char &o1,const unsigned char &o2, const unsigned char &o3, bool &a1, bool &a2, bool &a3) {
//save the bools b1,b2,b3 in the order given by the permutation o1,o2,o3
//oi must be an integer between 1 and 3. If o1,o2,o3 are not unique, garbage will also be produced 

    if(o1 == 1 ) {
        if(o2 == 2) {
            a1 = b1;
            a2 = b2;
            a3 = b3;
        }
        else {
            a1 = b1;
            a2 = b3;
            a3 = b2;
        }
    }
    else if (o1 == 2) {
        if(o2  == 1) {
            a1 = b2;
            a2 = b1;
            a3 = b3;
        }
        else{
            a1 = b2;
            a2 = b3;
            a3 = b1;
        }
    }
    else{
        if(o2 == 1) {
            a1 = b3;
            a2 = b1;
            a3 = b2;
        }
        else{
            a1 = b3;
            a2 = b2;
            a3 = b1;
        }
    }


}

inline void what_order(const unsigned char &o1, const unsigned char &o2, const unsigned char &o3, unsigned char &no1, unsigned char &no2, unsigned char &no3)
{
//what order is it necessary to arrange o1,o2,o3 to get 1 2 3
    if(o1 == 1) {
        if(o2 == 2) {
            no1 = 1;
            no2 = 2;
            no3 = 3;
        }
        else{
            no1 = 1;
            no2 = 3;
            no3 = 2;
        }
    }
    else if (o1 == 2) {
        if(o2 == 1) {
            no1 = 2;
            no2 = 1;
            no3 = 3;
        }
        else{
            no1 = 3;
            no2 = 1;
            no3 = 2;
        }
    }
    else {
        if(o2 == 1) {
            no1 = 2;
            no2 = 3;
            no3 = 1;
        }
        else{
            no1 = 3;
            no2 = 2;
            no3 = 1;
        }

    }

}

template <class T>
int sign(T y)
{ //return -1 if negative and 1 otherwise
    if ( y>0 ) return 1;
	return -1;
}

template <class T>
T absolute(T x) { //retrun the absolute value of a function
    if (x<0) return x*T(-1.0);
    else return x;
}


template <class T, class Q>
T MAX(Q x, T y) { //return the largest of the objects x and y (operator> must be defined)
if ( x>y ) return T(x);
else return y;

}

template <class T, class Q>
T MIN(Q x, T y) { //same as function above but vice versa
if ( x<y ) return T(x);
else return y;

}

inline double arctan(double y, double x) { //arctan(y/x) allowing for singularity
    if ( x == 0 ) return sign(y)*pi;
    if ( x > 0 ) return atan(y/x);
    if( x < 0 && y >= 0 ) return atan(y/x)+pi;
    if( x < 0 && y < 0 ) return atan(y/x)-pi;
            
    else return atan(y/x);
}

inline long double Heaviside(long double x) { //Heaviside step function
    if ( x < 0) return 0;
    else return 1;

}

inline long double arccos(long double x) {
    
    if (x > 1) return 0.;
    if ( x < -1) return pi;
    else return acos(x);
}

template <class T, class Q, class F>
struct combin {
        T *func1;
        Q *func2;
        combin(T *funcc1, Q *funcc2 ) : func1(funcc1), func2(funcc2) {}
        F operator()( F x1) {
            return (*func1)(x1)*(*func2)(x1);
        }  
};


template<class H, class T,class Q>
struct combine {
    T &func1;
    Q &func2;
    combine(T &funcc1, Q &funcc2) : func1(funcc1),func2(funcc2) {}
    H operator()(H x, H y) {
        return (func1)(x,y)*(func2)(x,y);
    }
    
    
    
};


    template <class T, class Q, class F>
    struct combinefuncs {
        T *func1;
        Q *func2;
        combinefuncs(T *funcc1, Q *funcc2 ) : func1(funcc1), func2(funcc2) {}
        F operator()( F x1, F x2) {
            return (*func1)(x1,x2)*(*func2)(x1,x2);
        }
    };

    template<class T, class Q>
     struct d1a {
        Q &func;
        T R;
        d1a(Q &funcc, T RR) : func(funcc),R(RR) {}               
        long double operator()(long double x0) {
            return func(x0,R-x0);
                       }
                   };
                   
    template<class T, class Q>
     struct d1a1 {
        Q &func;
        T R;
        d1a1(Q &funcc, T RR) : func(funcc),R(RR) {}               
        long double operator()(long double x0) {
            return (R-x0)*func(x0,R-x0);
                       }
                   };
                   
    template<class T, class Q>
     struct d1a2 {
        Q &func;
        T R,z;
        d1a2(Q &funcc, T RR, T zz) : func(funcc),R(RR), z(zz) {}               
        long double operator()(long double x0) {
            return (x0-R)*func(x0,z+R-x0);
                       }
                   };
    template<class T, class Q>
     struct d1a3 {
        Q &func;
        T R,z;
        d1a3(Q &funcc, T RR, T zz) : func(funcc),R(RR), z(zz) {}               
        long double operator()(long double x0) {
            return (x0)*func(x0,z-sqrt(SQR(R)-SQR(x0)))-(x0)*func(x0,z+sqrt(SQR(R)-SQR(x0)));
                       }
                   };
    template<class T, class Q>
     struct d1a4 {
        Q &func;
        T R,z;
        d1a4(Q &funcc, T RR, T zz) : func(funcc),R(RR), z(zz) {}               
        T operator()(T x0) {
            return ((x0/(sqrt(SQR(R)-SQR(x0))))*func(x0,z-sqrt(SQR(R)-SQR(x0))))+((x0/(sqrt(SQR(R)-SQR(x0))))*func(x0,z+sqrt(SQR(R)-SQR(x0))));
                       }
                   };

       template <class T, class Q>
       struct d1a5 {
           Q &func;
           T R,z;
           d1a5(Q &funcc,T RR, T zz) : func(funcc),R(RR),z(zz) {}
           T operator()(T x0) {
               return -((SQR(x0)/(sqrt(SQR(R)-SQR(x0))))*func(x0,z-sqrt(SQR(R)-SQR(x0))))-((SQR(x0)/(sqrt(SQR(R)-SQR(x0))))*func(x0,z+sqrt(SQR(R)-SQR(x0))));
           }
           
           
       };

template <class T>
T Power(T a, int b) { //number to an integer power
	if (b<0) { cout << "b is: " << b << endl; error("b must be >= 0");  return 0; }
    if (b==0) return 1;
    if (a==0) return 0;
    if (b%2==0) {
        return Power(a*a, b/2);
    } else if (b%2==1) {
        return a*Power(a*a,b/2);
    }
    return 0;
}

int return_csv_in_current_dir(string match, vector<string> &files) {
DIR *dir;
struct dirent *diread;
//vector<string> files;

if ((dir = opendir("./")) != nullptr)
{
    while ((diread = readdir(dir)) != nullptr)
    {
        std::string fname = diread->d_name;
        if (fname.find(match) != std::string::npos && fname.find(".csv") != std::string::npos)
            files.push_back(fname);

        //files.push_back(diread->d_name);
    }
    closedir(dir);
}
else
{
    perror("opendir");
    return EXIT_FAILURE;
}

return 0;

}

int return_csv_in_dir(string directory, string match, vector<string> &files)
{
    DIR *dir;
    struct dirent *diread;
    //vector<string> files;

    if ((dir = opendir(directory.c_str())) != nullptr)
    {
        while ((diread = readdir(dir)) != nullptr)
        {
            std::string fname = diread->d_name;
            if (fname.find(match) != std::string::npos && fname.find(".csv") != std::string::npos)
                files.push_back(fname);

            //files.push_back(diread->d_name);
        }
        closedir(dir);
    }
    else
    {
        perror("opendir");
        return EXIT_FAILURE;
    }

    return 0;
}

// void handler(int sig)
// {
//     void *array[10];
//     size_t size;

//     // get void*'s for all entries on the stack
//     size = backtrace(array, 10);

//     // print out all the frames to stderr
//     fprintf(stderr, "Error: signal %d:\n", sig);
//     backtrace_symbols_fd(array, size, STDERR_FILENO);
//     exit(1);
// }

#endif	/* BASIC_H */




