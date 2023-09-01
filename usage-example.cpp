#include <iostream>
#include "bessel-library.hpp"

using namespace std;

int main()
{

    int n;
    complex<double> z;
    cout << scientific;
    cout.precision(15);

    // Jn(z), n=0, z=3.2-i1.7
    n = 0;
    z = complex<double> (3.2, -1.7);
    cout << bessel::cyl_j(n, z) << endl;

    // Jn(z), n=0, z=-15.1-i83.9
    n = 0;
    z = complex<double> (-15.1, -83.9);
    cout << bessel::cyl_j(n, z) << endl;

    // Jn(z), n=-1, z=-2.8+i7.7
    n = -1;
    z = complex<double> (-2.8, 7.7);
    cout << bessel::cyl_j(n, z) << endl;

    // Jn(z), n=3, z=4+i29
    n = 3;
    z = complex<double> (4., 29.);
    cout << bessel::cyl_j(n, z) << endl;

    // Jn(z), n=-14, z=42-i8
    n = -14;
    z = complex<double> (42., -8.);
    cout << bessel::cyl_j(n, z) << endl;

    // Jn(z), n=100, z=71+i30
    n = 100;
    z = complex<double> (71., 30.);
    cout << bessel::cyl_j(n, z) << endl;

    // Jn'(z), n=0, z=-11+i18
    n = 0;
    z = complex<double> (-11., 18.);
    cout << bessel::cyl_j_diff(n, z) << endl;

    // Jn'(z), n=-9, z=-6.6-i3.6
    n = -9;
    z = complex<double> (-6.6, -3.6);
    cout << bessel::cyl_j_diff(n, z) << endl;

    // Jn'(z), n=47, z=2.71+i3.14
    n = 47;
    z = complex<double> (2.71, 3.14);
    cout << bessel::cyl_j_diff(n, z) << endl;

}