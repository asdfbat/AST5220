#include <math.h>
#include <iostream>
#include <vector>
#include "Utils.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

using Vector = std::vector<double>;
using Vector2D = std::vector<Vector>;


int main(){
    /*
    TESTING SPLINE STUFF
    */

    auto myfunc = [&](double x){return sin(x); };
    auto myfunc_ddx = [&](double x){return cos(x); };

    int n = 10;
    Vector x = Utils::linspace(0.0, 10.0, n);
    Vector y(n);

    for(int i=0; i<n; i++){
        y[i] = myfunc(x[i]);
    }

    // plt::plot(x, y, "o");
    // plt::show();

    Spline f_spline(x, y, "myspline");

    Vector x_dense = Utils::linspace(0.0, 10.0, 10*n);
    Vector y_dense_true(10*n);
    Vector y_dense_spline(10*n);

    for(int i=0; i<10*n; i++){
        y_dense_true[i] = myfunc(x_dense[i]);
        y_dense_spline[i] = f_spline(x_dense[i]);
    }

    // plt::plot(x_dense, y_dense_spline);
    // plt::plot(x_dense, y_dense_true);
    // plt::show();

    Vector y_ddx_spline(10*n);
    Vector y_ddx_true(10*n);

    for(int i=0; i<10*n; i++){
        y_ddx_spline[i] = f_spline.deriv_x(x_dense[i]);
        y_ddx_true[i] = myfunc_ddx(x_dense[i]);
    }

    // plt::plot(x_dense, y_ddx_spline);
    // plt::plot(x_dense, y_ddx_true);
    // plt::show();


    /*
    TESTING ODESOLVER STUFF
    */

    ODEFunction somefunc = [&] (double x, const double *y, double *dydx){
        dydx[0] = y[1];
        dydx[1] = -y[0];
        return GSL_SUCCESS;
    };

    Vector y_inc = {1.5, -0.5};

    ODESolver ode;
    ode.solve(somefunc, x_dense, y_inc);

    Vector2D all_data = ode.get_data();
    Vector solution(10*n);
    for(int i=0; i<10*n; i++){
        solution[i] = all_data[i][1];
    }

    std::cout << all_data[0].size() << std::endl;

    plt::plot(x_dense, solution);
    plt::show();


    return 0;
}