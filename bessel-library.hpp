#pragma once

#define C_1_SQRTPI 0.5641895835477562869480794515607725L

#include <complex>
#include <cmath>

namespace bessel
{

template<typename T>
T cyl_j( const int _n, const T _x )
{
    int n = _n;
    T n_sign = T(1);
    if ( _n < 0 )
    {
        n = -_n;
        n_sign = T(-1);
    }

    T x = _x;
    T x_sign = T(1);
    if ( _x < T(0) )
    {
        x = -_x;
        x_sign = T(-1);
    }

    // Implementation from C++17 standard library.
    return std::pow(x_sign, n)*std::pow(n_sign, n)*std::cyl_bessel_j(n, x);
}

template<typename T>
T cyl_j_diff( const int _n, const T _x )
{
    if ( _n == 0 )
    {
        // Implementation from C++17 standard library.
        if ( _x >= 0 )
        {
            // Since J0'(x)=-J1(x).
            return -std::cyl_bessel_j(1, _x);
        }
        else
        {
            // Since J0'(-x)=J1(x).
            return std::cyl_bessel_j(1, -_x);
        }
    }
    else
    {
        // Implementation from C++17 standard library.
        return cyl_j(_n-1,_x)-T(_n)*cyl_j(_n,_x)/_x;
    }
}

template<typename T>
static std::complex<T> __ascending_series_cyl_j0( const std::complex<T> _z )
{
    // Ascending Series from G. N. Watson 'A treatise on the
    //  theory of Bessel functions', 2ed, Cambridge, 1996,
    //  Chapter II, equation (3); or from Equation 9.1.12 of
    //  M. Abramowitz, I. A. Stegun 'Handbook of Mathematical
    //  Functions'.
    const T epsilon = std::numeric_limits<T>::epsilon();
    std::complex<T> j0 = T(1);
    std::complex<T> sm = T(1);
    for ( int m = 1; std::abs(sm/j0) >= epsilon; m++ )
    {
        sm *= - _z*_z * T(0.25) / ( T(m)*T(m) );
        j0 += sm;
    }
    return j0;
}

template<typename T>
static std::complex<T> __ascending_series_cyl_j1( const std::complex<T> _z )
{
    // Ascending Series from G. N. Watson 'A treatise on the
    //  theory of Bessel functions', 2ed, Cambridge, 1996,
    //  Chapter II, equation (3); or from Equation 9.1.12 of
    //  M. Abramowitz, I. A. Stegun 'Handbook of Mathematical
    //  Functions'.
    const T epsilon = std::numeric_limits<T>::epsilon();
    std::complex<T> j1 = T(1);
    std::complex<T> sm = T(1);
    for ( int m = 1; std::abs(sm/j1) >= epsilon; m++ )
    {
        sm *= - _z*_z * T(0.25) / ( T(m)*T(m+1) );
        j1 += sm;
    }
    j1 *= T(0.5)*_z;
    return j1;
}

template<typename T>
static std::complex<T> __semiconvergent_series_cyl_j0( 
                            const std::complex<T> _z, const int _m_max )
{
    // Stokes Semiconvergent Series from A. Gray, G. B. Mathews 'A
    //  treatise on Bessel functions and their applications to
    //  physics, 1895.
    std::complex<T> Pm = T(1);
    std::complex<T> Qm = T(0.125)/_z;
    std::complex<T> P = Pm;
    std::complex<T> Q = Qm;
    for ( int m=1; m<=_m_max; m++ )
    {
        T pim = T(4*m-3)*T(4*m-3)*T(4*m-1)*T(4*m-1) / ( T(2*m-1)*T(128*m) );
        Pm = -Pm*pim/(_z*_z);

        T xim = T(4*m-1)*T(4*m-1)*T(4*m+1)*T(4*m+1) / ( T(2*m+1)*T(128*m) );
        Qm = -Qm*xim/(_z*_z);

        P += Pm;
        Q += Qm;
    }
    return T(C_1_SQRTPI)*(std::cos(_z)*(P-Q) + std::sin(_z)*(P+Q))/std::sqrt(_z);
}

template<typename T>
static std::complex<T> __semiconvergent_series_cyl_j1(
                            const std::complex<T> _z, const int _m_max )
{
    // Stokes Semiconvergent Series from A. Gray, G. B. Mathews 'A
    //  treatise on Bessel functions and their applications to
    //  physics, 1895.
    std::complex<T> Pm = T(1);
    std::complex<T> Qm = T(-0.375)/_z;
    std::complex<T> P = Pm;
    std::complex<T> Q = Qm;
    for ( int m=1; m<=_m_max; m++ )
    {
        T pim = (T(4)-T(4*m-3)*T(4*m-3))*(T(4)-T(4*m-1)*T(4*m-1))/(T(2*m-1)*T(128*m));
        Pm = -Pm*pim/(_z*_z);

        T xim = (T(4*m-1)*T(4*m-1)-T(4))*(T(4*m+1)*T(4*m+1)-T(4))/(T(2*m+1)*T(128*m));
        Qm = -Qm*xim/(_z*_z);

        P += Pm;
        Q += Qm;
    }
    return T(C_1_SQRTPI)*(-std::cos(_z)*(P+Q) + std::sin(_z)*(P-Q))/std::sqrt(_z);
}

template<typename T>
static std::complex<T> __cyl_j0( const std::complex<T> _z )
{
    std::complex<T> z = _z;

    if ( std::real(_z) < T(0) )
    {
        z = -_z;
    }

    if ( std::abs(z) == T(0) )
    {
        return std::complex<T> (T(1), T(0));
    }
    else if ( std::abs(z) <= T(14.201) )
    {
        return __ascending_series_cyl_j0(z);
    }
    else if ( std::abs(z) < T(14.35) )
    {
        return __semiconvergent_series_cyl_j0(z, 13);
    }
    else if ( std::abs(z) < T(14.5) )
    {
        return __semiconvergent_series_cyl_j0(z, 14);
    }
    else if ( std::abs(z) < T(15.5) )
    {
        return __semiconvergent_series_cyl_j0(z, 15);
    }
    else if ( std::abs(z) < T(16) )
    {
        return __semiconvergent_series_cyl_j0(z, 16);
    }
    else if ( std::abs(z) < T(16.5) )
    {
        return __semiconvergent_series_cyl_j0(z, 16.5);
    }
    else if ( std::abs(z) < T(17) )
    {
        return __semiconvergent_series_cyl_j0(z, 17);
    }
    else if ( std::abs(z) < T(17.5) )
    {
        return __semiconvergent_series_cyl_j0(z, 17.5);
    }
    else if ( std::abs(z) < T(35) )
    {
        return __semiconvergent_series_cyl_j0(z, 12);
    }
    else if ( std::abs(z) < T(50) )
    {
        return __semiconvergent_series_cyl_j0(z, 10);
    }
    else
    {
        return __semiconvergent_series_cyl_j0(z, 8);
    }
}

template<typename T>
static std::complex<T> __cyl_j1( const std::complex<T> _z )
{
    std::complex<T> z = _z;

    T z_sign = T(1);
    if ( std::real(_z) < T(0) )
    {
        z = -_z;
        z_sign = T(-1);
    }

    if ( std::abs(z) == T(0) )
    {
        return std::complex<T> (T(0), T(0));
    }
    else if ( std::abs(z) <= T(14.201) )
    {
        return z_sign*__ascending_series_cyl_j1(z);
    }
    else if ( std::abs(z) < T(14.35) )
    {
        return z_sign*__semiconvergent_series_cyl_j1(z, 13);
    }
    else if ( std::abs(z) < T(14.5) )
    {
        return z_sign*__semiconvergent_series_cyl_j1(z, 14);
    }
    else if ( std::abs(z) < T(15.5) )
    {
        return z_sign*__semiconvergent_series_cyl_j1(z, 15);
    }
    else if ( std::abs(z) < T(16) )
    {
        return z_sign*__semiconvergent_series_cyl_j1(z, 16);
    }
    else if ( std::abs(z) < T(16.5) )
    {
        return z_sign*__semiconvergent_series_cyl_j1(z, 16.5);
    }
    else if ( std::abs(z) < T(17) )
    {
        return z_sign*__semiconvergent_series_cyl_j1(z, 17);
    }
    else if ( std::abs(z) < T(17.5) )
    {
        return z_sign*__semiconvergent_series_cyl_j1(z, 17.5);
    }
    else if ( std::abs(z) < T(35) )
    {
        return z_sign*__semiconvergent_series_cyl_j1(z, 12);
    }
    else if ( std::abs(z) < T(50) )
    {
        return z_sign*__semiconvergent_series_cyl_j1(z, 10);
    }
    else
    {
        return z_sign*__semiconvergent_series_cyl_j1(z, 8);
    }
}

template<typename T>
static inline T __minus_log10_cyl_j_at_infinity( const T _n, const T _abs_z )
{
    // Auxiliary function to calculate -log( Jn(x->INF) ).
    return T(0.5)*std::log10(T(6.28)*_n) - _n*std::log10(T(1.36)*_abs_z/_n);
}

template<typename T>
static int __ini_for_br_1( const T _abs_z, const T _mg  )
{
    // Starting point for backward recurrence
    //  for when |Jn(x)|~10e-mg
    //  using the secant method.
    T n0 = std::round( T(1.1)*_abs_z + T(1) );
    T f0 = __minus_log10_cyl_j_at_infinity(n0, _abs_z) - _mg;
    T n1 = n0 + T(5);
    T f1 = __minus_log10_cyl_j_at_infinity(n1, _abs_z) - _mg;
    T nn = n1 - (n1 - n0)/(T(1) - f0/f1);
    T f = __minus_log10_cyl_j_at_infinity(nn, _abs_z) - _mg;
    while ( std::abs(nn - n1) >= T(1) )
    {
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
        nn = n1 - (n1 - n0)/(T(1) - f0/f1);
        f = __minus_log10_cyl_j_at_infinity(nn, _abs_z) - _mg;
    }
    return int(nn);   
}

template<typename T>
static int __ini_for_br_2( const T _n, const T _abs_z, const T _sd  )
{
    // Starting point for backward recurrence
    //  for when Jn(x) has sd significant digits
    //  using the secant method.
    T n0, n1, nn, obj, f0, f1, f;
    T hmp = T(0.5)*_sd;
    T ejn = __minus_log10_cyl_j_at_infinity(_n, _abs_z);
    if ( ejn <= hmp )
    {
        obj = _sd;
        n0 = std::round( T(1.1)*_abs_z );
    }
    else
    {
        obj = hmp + ejn;
        n0 = _n;
    }
    f0 = __minus_log10_cyl_j_at_infinity(n0, _abs_z) - obj;
    n1 = n0 + T(5);
    f1 = __minus_log10_cyl_j_at_infinity(n1, _abs_z) - obj;
    nn = n1 - (n1-n0)/(T(1)-f0/f1);
    f = __minus_log10_cyl_j_at_infinity(nn, _abs_z)-obj;
    while ( std::abs(nn - n1) >= T(1) )
    {
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
        nn = n1 - (n1-n0)/(T(1)-f0/f1);
        f = __minus_log10_cyl_j_at_infinity(nn, _abs_z)-obj;
    }
    return int( nn + T(10) );
}

template<typename T>
static std::complex<T> __forward_recurrence_cyl_j( const int _n, const std::complex<T> _z )
{
    std::complex<T> b0 = __cyl_j0( _z );
    std::complex<T> b1 = __cyl_j1( _z );
    std::complex<T> bn;
    for ( int k=2; k<=_n; k++ )
    {
        bn = T(2)*T(k-1)*b1/_z-b0;
        b0 = b1;
        b1 = bn;
    }
    return bn;
}

template<typename T>
static std::complex<T> __backward_recurrence_cyl_j( const int _n, const std::complex<T> _z )
{
    std::complex<T> bn;

    int nm = _n;
    int m = __ini_for_br_1( std::abs(_z), T(std::numeric_limits<T>::max_exponent10) );
    if ( m < _n )
    {
        nm = m;
    }
    else
    {
        m = __ini_for_br_2( T(_n), std::abs(_z), T(std::numeric_limits<T>::max_digits10) );
    }
  
    std::complex<T> cf2 = std::complex<T>(T(0), T(0));
    std::complex<T> cf1 = std::complex<T>(1.0e-100, T(0));
    std::complex<T> cf;

    for ( int k=m; k>=0; k-- )
    {
        cf = T(2)*T(k+1)*cf1/_z-cf2;
        if ( k == nm )
        {
            bn = cf;
        }
        cf2 = cf1;
        cf1 = cf;
    }
    
    std::complex<T> cs;
    std::complex<T> j0 = __cyl_j0(_z);
    std::complex<T> j1 = __cyl_j1(_z);
    if ( std::abs(j0) > std::abs(j1) )
    {
        cs = j0/cf;
    }
    else
    {
        cs = j1/cf2;
    }
    
    if ( !std::isnormal(abs(cs)) )
    {
        bn = __forward_recurrence_cyl_j(_n,_z);
    }
    else
    {
        bn *= cs;
    }

    return bn;
}

template<typename T>
std::complex<T> cyl_j( const int _n, const std::complex<T> _z, const bool _warnings=true )
{
    if ( abs(imag(_z)) > T(21) && _warnings )
    {
        std::cerr << "Warning: '|imag(_z)| > 21' in calculating 'cyl_j(" << _n << ", " << _z << ")'; accuracy may be lost."  << std::endl;
    }

    if ( _n == 0 )
    {
        return __cyl_j0( _z );
    }
    else if ( _n == 1 )
    {
        return __cyl_j1( _z );
    }
    else if ( _n == -1 )
    {
        return -__cyl_j1( _z );
    }
    else if ( std::abs(std::abs(_z) - T(0)) <= std::numeric_limits<T>::epsilon() )
    {
        return std::complex<T> ( T(0), T(0) );
    }
    else
    {
        int n = _n;
        T n_sign = T(1);
        if ( _n < 0 )
        {
            n = -_n;
            n_sign = T(-1);
        }

        std::complex<T> z = _z;
        T z_sign = T(1);
        if ( std::real(_z) < T(0) )
        {
            z = -_z;
            z_sign = T(-1);
        }
        
        if ( n < std::abs(z)*T(0.25) )
        {
            return std::pow(z_sign, n)*std::pow(n_sign, n)*__forward_recurrence_cyl_j(n, z);
        }
        else
        {
            return std::pow(z_sign, n)*std::pow(n_sign, n)*__backward_recurrence_cyl_j(n, z);
        }
    }
}

template<typename T>
std::complex<T> cyl_j_diff( const int _n, const std::complex<T> _z, const bool _warnings=true )
{
    if ( _n == 0 )
    {
        // Since dJ0(z)=-J1(z)
        return -cyl_j(1, _z, _warnings );
    }
    else
    {
        return cyl_j(_n-1,_z,_warnings)-T(_n)*cyl_j(_n,_z,_warnings)/_z;
    }
}

}
