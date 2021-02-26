#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "vector1.h"
#include "matrix2.h"

template <class T>
struct Bilinear2
{
    int n, m;
    vector1<T> x1, x2;
    matrix<T> y;
    bool flag;
    Bilinear2() : y(matrix<T>(2, 2))
    {
        n = 2;
        m = 2;
        vector1<T> p1(2);
        vector1<T> p2(2);
        x1 = p1;
        x2 = p2;

        //y = p3;
        flag = true;
    }

    Bilinear2(vector1<T> &x11, vector1<T> &x22, matrix<T> &yy, bool fllag = false) : n(x11.getsize()), m(x22.getsize()), x1(x11), x2(x22), y(yy)
    {
        flag = fllag;
        if (x11.getsize() != yy.getnrows())
            error("set of grid points do not agree with data dimension");
        if (x22.getsize() != yy.getncols())
            error("set of grid points do not agree with data dimension");
    }

    Bilinear2(const Bilinear2<T> &a) : n(a.n), m(a.m), x1(a.x1), x2(a.x2), y(a, y), flag(a.flag)
    {
        if (x1.getsize() != y.getnrows())
            error("set of grid points do not agree with data dimension");
        if (x2.getsize() != y.getncols())
            error("set of grid points do not agree with data dimension");
        cout << "this function called" << endl;
    }

    Bilinear2 &operator=(const Bilinear2 &a)
    {
        cout << "function called" << endl;
        n = a.n;
        m = a.m;

        x1 = a.x1;

        x2 = a.x2;
        y = a.y;
        flag = a.flag;
        return *this;
    }

    void PlotBilinear(int order = 2)
    {
        matrix<T> temp(order * n, order * m);
        vector1<T> xpos(order * n);
        vector1<T> ypos(order * m);

        xpos[0] = x1[0];
        xpos[order * n - 1] = x1[n - 1];

        ypos[0] = x2[0];
        ypos[order * m - 1] = x2[m - 1];
        double xgap = (xpos[order * n - 1] - xpos[0]) / (order * n - 1);
        double ygap = (ypos[order * m - 1] - ypos[0]) / (order * m - 1);

        for (int i = 1; i < order * n - 1; i++)
        {
            xpos[i] = xpos[0] + i * xgap;
        }
        for (int i = 1; i < order * n - 1; i++)
        {
            ypos[i] = ypos[0] + i * ygap;
        }
        for (int i = 0; i < order * n; i++)
        {
            for (int j = 0; j < order * m; j++)
            {
                temp(i, j) = (*this)(xpos[i], ypos[j]);
            }
        }

        cout << temp << endl;
    }

    void PlotBilinear(ofstream &myfile, int order = 2)
    {
        matrix<T> temp(order * n, order * m);
        vector1<T> xpos(order * n);
        vector1<T> ypos(order * m);

        xpos[0] = x1[0];
        xpos[order * n - 1] = x1[n - 1];

        ypos[0] = x2[0];
        ypos[order * m - 1] = x2[m - 1];
        double xgap = (xpos[order * n - 1] - xpos[0]) / (order * n - 1);
        double ygap = (ypos[order * m - 1] - ypos[0]) / (order * m - 1);

        for (int i = 1; i < order * n - 1; i++)
        {
            xpos[i] = xpos[0] + i * xgap;
        }
        for (int i = 1; i < order * n - 1; i++)
        {
            ypos[i] = ypos[0] + i * ygap;
        }
        for (int i = 0; i < order * n; i++)
        {
            for (int j = 0; j < order * m; j++)
            {
                temp(i, j) = (*this)(xpos[i], ypos[j]);
            }
        }

        myfile << "tempx=" << fixed << temp << endl;
        myfile << "ListPlot3D[tempx,PlotRange->All]" << endl;
    }

    T operator()(T xa, T xb)
    {
        if (xa < x1[0] || xa > x1[n - 1])
            return 0.0;
        if (xb < x2[0] || xb > x2[m - 1])
        {
            if (!flag)
                return 0.0;
            else
            {
                int j, indx1, indx2;
                if (fabs(xa - x1[0]) > fabs(xa - x1[n - 1]))
                {

                    for (int i = n - 1; i >= 0; i--)
                    {
                        T temp = xa - x1[i];
                        if (temp > 0)
                        {
                            j = i;
                            break;
                        }
                    }
                    //if(fabs(xa-x1[j]) > fabs(xa-x1[j+1])) {
                    indx1 = j;
                    indx2 = j + 1;
                    //}
                    //            else {
                    //                indx1 = j-1;
                    //                indx2 = j;
                    //            }
                }
                else
                {
                    for (int i = 0; i <= n - 1; i++)
                    {
                        T temp = xa - x1[i];
                        if (temp < 0)
                        {
                            j = i;
                            break;
                        }
                    }
                    indx1 = j - 1;
                    indx2 = j;
                }

                T f1, f2;
                if (xb < x2[0])
                {
                    f1 = y.mat[indx1 * m];
                    f2 = y.mat[indx2 * m];
                }
                else
                {
                    f1 = y.mat[indx1 * m + m - 1];
                    f2 = y.mat[indx2 * m + m - 1];
                }

                return f1 + (xa - x1[indx1]) * ((f2 - f1) / (x1[indx2] - x1[indx1]));
            }
        }

        else
        {

            int j, indx1, indx2;
            int k, indx3, indx4;
            if (fabs(xa - x1[0]) > fabs(xa - x1[n - 1]))
            {
                for (int i = n - 1; i >= 0; i--)
                {
                    T temp = xa - x1[i];
                    if (temp > 0)
                    {
                        j = i;
                        break;
                    }
                }
                //if(fabs(xa-x1[j]) > fabs(xa-x1[j+1])) {
                indx1 = j;
                indx2 = j + 1;
                //}
                //            else {
                //                indx1 = j-1;
                //                indx2 = j;
                //            }
            }
            else
            {
                for (int i = 0; i <= n - 1; i++)
                {
                    T temp = xa - x1[i];
                    if (temp < 0)
                    {
                        j = i;
                        break;
                    }
                }
                indx1 = j - 1;
                indx2 = j;
            }

            if (fabs(xb - x2[0]) > fabs(xb - x2[m - 1]))
            {
                for (int i = m - 1; i >= 0; i--)
                {
                    T temp = xb - x2[i];
                    if (temp > 0)
                    {
                        k = i;
                        break;
                    }
                }
                //if(fabs(xa-x1[j]) > fabs(xa-x1[j+1])) {
                indx3 = k;
                indx4 = k + 1;
                //}
                //            else {
                //                indx1 = j-1;
                //                indx2 = j;
                //            }
            }
            else
            {
                for (int i = 0; i <= m - 1; i++)
                {
                    T temp = xb - x2[i];
                    if (temp < 0)
                    {
                        k = i;
                        break;
                    }
                }
                indx3 = k - 1;
                indx4 = k;
            }

            T xi = x1[indx1];
            T xj = x1[indx2];
            T yi = x2[indx3];
            T yj = x2[indx4];
            T f0 = y.mat[indx1 * m + indx3];
            T f1 = y.mat[indx1 * m + indx4];
            T f2 = y.mat[indx2 * m + indx3];
            T f3 = y.mat[indx2 * m + indx4];
            T denom = ((xj - xi) * (yj - yi));

            return (-(f1 * (xa - xj) * (xb - yi)) + f0 * (xa - xj) * (xb - yj) +
                    (xa - xi) * (f3 * (xb - yi) + f2 * (-xb + yj))) /
                   ((xi - xj) * (yi - yj));
            //        cout << x1[indx1] << "," << xa << "," << x1[indx2] << endl;
            //        cout << x2[indx3] << "," << xb << "," << x2[indx4] << endl;
        }
    }
};


struct Bilinear3
{

    int ncol;
    double n, m;
    double xmax,ymax;
    double gamma_white,gamma_black;
    matrix<double> data;
    Bilinear3() : data(matrix<double>(2, 2))
    {
        n = 2;
        m = 2;
        xmax = 2.;
        ymax = 2.;
        gamma_white = 1.;
        gamma_black = 1.;
        ncol = int(m);

    }

    Bilinear3(double xmaxx, double ymaxx, const char *wherebmp) : xmax(xmaxx), ymax(ymaxx), data(matrix<double>(2, 2))
    {
        BMP a(wherebmp);


        matrix<double> b(a.bmp_info_header.width, a.bmp_info_header.height);

        int iter = 0;
        for (int i = 0; i < a.bmp_info_header.height; i++)
        {
            for (int j = 0; j < a.bmp_info_header.width; j++)
            {
                b(i, j) = a.data[iter];
                iter++;
            }
        }

        data =  b;

        n = (double)a.bmp_info_header.height;
        m = (double)a.bmp_info_header.width;
        ncol = a.bmp_info_header.width;

        gamma_white = 1.;
        gamma_black = 1.;
    }

    Bilinear3(const Bilinear3 &a) : xmax(a.xmax), ymax(a.ymax), data(a.data) , n(a.n), m(a.m), gamma_black(a.gamma_black), gamma_white(a.gamma_white)
    {
        ncol = a.ncol;
    }

    Bilinear3 &operator=(const Bilinear3 &a)
    {
        cout << "function called" << endl;
        n = a.n;
        m = a.m;

        xmax = a.xmax;

        ymax = a.ymax;
        data = a.data;
        gamma_white = a.gamma_white;
        gamma_black = a.gamma_black;
        ncol = a.ncol;
        return *this;
    }

    double operator()(double x, double y) {
        int n1 = floor(n*x/xmax);
        int m1 = floor(m*y/ymax);
        double val = data.mat[n1*ncol+m1];
        if(val >0.5) {
            return gamma_black;
        }
        else {
            return gamma_white;
        }



    }

};

#endif /* INTERPOLATION_H */
