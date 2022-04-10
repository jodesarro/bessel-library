#include <iostream>
#include "bessel-library.hpp"

using namespace std;

int main()
{

    complex<double> z;
    cout << scientific;
    cout.precision(16);

    // z=1+i0
    z = complex<double> (1, 0);
    cout << bessel::cyl_j0( z ) << endl;

    // z=5+i0
    z = complex<double> (5, 0);
    cout << bessel::cyl_j0( z ) << endl;

     // z=10+i0
    z = complex<double> (10, 0);
    cout << bessel::cyl_j0( z ) << endl;

    // z=25+i0
    z = complex<double> (25, 0);
    cout << bessel::cyl_j0( z ) << endl;

    // z=50+i0
    z = complex<double> (50, 0);
    cout << bessel::cyl_j0( z ) << endl;

    // z=100+i0
    z = complex<double> (100, 0);
    cout << bessel::cyl_j0( z ) << endl;

    // z=4+i2
    z = complex<double> (4, 2);
    cout << bessel::cyl_j0( z ) << endl;

    // z=20+i10
    z = complex<double> (20, 10);
    cout << bessel::cyl_j0( z ) << endl;

}
