#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <iomanip>
#include <variant>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <complex>
#include <map>


Matrix2d s0g()
{
    Matrix2d m;
    m << pow(Nc, 2), 0,
        0, CF * Nc / 2.;
    return m;
};

typedef Matrix<double, 9, 9> Matrix9d;
typedef Matrix<std::complex<double>, 9, 9> Matrix9cd;
Matrix9d s4g()
{
    Matrix9d m;
    double C1 = pow(Nc, 4) - 3 * pow(Nc, 2) + 3;
    double C2 = 3 - pow(Nc, 2);
    double C3 = 3 + pow(Nc, 2);
    double V = pow(Nc, 2) - 1;

    m << C1, C2, C2, C2, C2, C3, Nc * V, -Nc, Nc * V,
        C2, C1, C2, C2, C3, C2, Nc * V, Nc * V, -Nc,
        C2, C2, C1, C3, C2, C2, -Nc, Nc * V, Nc * V,
        C2, C2, C3, C1, C2, C2, -Nc, Nc * V, Nc * V,
        C2, C3, C2, C2, C1, C2, Nc * V, Nc * V, -Nc,
        C3, C2, C2, C2, C2, C1, Nc * V, -Nc, Nc * V,
        Nc * V, Nc * V, -Nc, -Nc, Nc * V, Nc * V, pow(Nc, 2) * V, pow(Nc, 2), pow(Nc, 2),
        -Nc, Nc * V, Nc * V, Nc * V, Nc * V, -Nc, pow(Nc, 2), pow(Nc, 2) * V, pow(Nc, 2),
        Nc * V, -Nc, Nc * V, Nc * V, -Nc, Nc * V, pow(Nc, 2), pow(Nc, 2), pow(Nc, 2) * V;
    m = V / pow(Nc, 2) * m;
    return m;
};

Matrix3d s2g()
{
    Matrix3d m;
    m << Nc, 0, 0,
        0, Nc / 2., 0,
        0, 0, (pow(Nc, 2) - 4) / (2. * Nc);
    m = (pow(Nc, 2) - 1.) * m;
    return m;
};


complex<double> IM{0, 1.};

MatrixXcd MH_m_qqb(double r)
{
    Matrix2cd m;

    m << 0, (CF * log(r / (1 - r))) / Nc,
        2. * log(r / (1 - r)), Nc * (IM * PI + log(r)) - (2. * log(r / (1 - r))) / Nc;
    return m * Gamma_cusp[0];
};

MatrixXcd MH_m_qq(double r)
{
    Matrix2cd m;

    m << 2. * CF * (IM * PI + log(r)) - 2. * CF * log(r / (1 - r)), (CF * (IM * PI + log(r))) / Nc,
        2. * (IM * PI + log(r)), ((-3 + pow(Nc, 2)) * (IM * PI + log(r))) / Nc + log(r / (1 - r)) / Nc;
    return m * Gamma_cusp[0];
};

MatrixXcd MH_m_ggqqb(double r)
{
    Matrix3cd m;
    m << 0, log(r / (1 - r)), 0,
        2. * log(r / (1 - r)), Nc * (IM * PI + log(r)) - (Nc * log(r / (1 - r))) / 2., ((-4 + pow(Nc, 2)) * log(r / (1 - r))) / (2. * Nc),
        0, (Nc * log(r / (1 - r))) / 2., Nc * (IM * PI + log(r)) - (Nc * log(r / (1 - r))) / 2.;
    return m * Gamma_cusp[0];
};

MatrixXcd MH_m_qgqg(double r)
{
    Matrix3cd m;
    m << (CF + Nc) * (IM * PI + log(r)), (-IM) * PI - log(r) + log(r / (1 - r)), 0,
        -2. * (IM * PI + log(r)) + 2. * log(r / (1 - r)), (-0.5 * 1 / Nc + Nc) * (IM * PI + log(r)) - (Nc * log(r / (1 - r))) / 2., ((4 - pow(Nc, 2)) * (IM * PI + log(r))) / (2. * Nc) + ((-4 + pow(Nc, 2)) * log(r / (1 - r))) / (2. * Nc),
        0, -0.5 * (Nc * (IM * PI + log(r))) + (Nc * log(r / (1 - r))) / 2., (-0.5 * 1 / Nc + Nc) * (IM * PI + log(r)) - (Nc * log(r / (1 - r))) / 2.;
    return m * Gamma_cusp[0];
};

MatrixXcd MH_m_qggq(double r)
{
    Matrix3cd m;
    m << (CF + Nc) * (IM * PI + log(1 - r)), (-IM) * PI - log(1 - r) + log((1 - r) / r), 0,
        -2. * (IM * PI + log(1 - r)) + 2. * log((1 - r) / r), (-0.5 * 1 / Nc + Nc) * (IM * PI + log(1 - r)) - (Nc * log((1 - r) / r)) / 2., ((4 - pow(Nc, 2)) * (IM * PI + log(1 - r))) / (2. * Nc) + ((-4 + pow(Nc, 2)) * log((1 - r) / r)) / (2. * Nc),
        0, -0.5 * (Nc * (IM * PI + log(1 - r))) + (Nc * log((1 - r) / r)) / 2., (-0.5 * 1 / Nc + Nc) * (IM * PI + log(1 - r)) - (Nc * log((1 - r) / r)) / 2.;
    return m * Gamma_cusp[0];
};

MatrixXcd MH_m_gggg(double r)
{
    Matrix9cd m;
    m << Nc * (IM * PI + log(r)), 0, 0, 0, 0, 0, log(r / (1 - r)), 0, (-IM) * PI - log(r) + log(r / (1 - r)),
        0, Nc * (IM * PI + log(r)) - Nc * log(r / (1 - r)), 0, 0, 0, 0, -log(r / (1 - r)), (-IM) * PI - log(r), 0,
        0, 0, 2. * Nc * (IM * PI + log(r)) - Nc * log(r / (1 - r)), 0, 0, 0, 0, IM * PI + log(r), IM * PI + log(r) - log(r / (1 - r)),
        0, 0, 0, 2. * Nc * (IM * PI + log(r)) - Nc * log(r / (1 - r)), 0, 0, 0, IM * PI + log(r), IM * PI + log(r) - log(r / (1 - r)),
        0, 0, 0, 0, Nc * (IM * PI + log(r)) - Nc * log(r / (1 - r)), 0, -log(r / (1 - r)), (-IM) * PI - log(r), 0,
        0, 0, 0, 0, 0, Nc * (IM * PI + log(r)), log(r / (1 - r)), 0, (-IM) * PI - log(r) + log(r / (1 - r)),
        (-IM) * PI - log(r) + log(r / (1 - r)), (-IM) * PI - log(r), 0, 0, (-IM) * PI - log(r), (-IM) * PI - log(r) + log(r / (1 - r)), 0, 0, 0,
        0, -log(r / (1 - r)), IM * PI + log(r) - log(r / (1 - r)), IM * PI + log(r) - log(r / (1 - r)), -log(r / (1 - r)), 0, 0, 2. * Nc * (IM * PI + log(r)) - 2. * Nc * log(r / (1 - r)), 0,
        log(r / (1 - r)), 0, IM * PI + log(r), IM * PI + log(r), 0, log(r / (1 - r)), 0, 0, 2. * Nc * (IM * PI + log(r));
    return m * Gamma_cusp[0];
};

MatrixXd f0_qqbQQb(double r)
{
    Matrix2d m;

    m << 0, 0,
        0, 2 - 4. * r + 4. * pow(r, 2);
    return m;
};

MatrixXd f0_qQbqQb(double r)
{
    Matrix2d m;

    m << 0.5 + 1 / (2. * pow(Nc, 4)) - pow(Nc, -2) + pow(r, -2) + 1 / (pow(Nc, 4) * pow(r, 2)) - 2 / (pow(Nc, 2) * pow(r, 2)) - 1 / r - 1 / (pow(Nc, 4) * r) + 2 / (pow(Nc, 2) * r), pow(Nc, -3) - 1 / Nc + 2 / (pow(Nc, 3) * pow(r, 2)) - 2 / (Nc * pow(r, 2)) - 2 / (pow(Nc, 3) * r) + 2 / (Nc * r),
        pow(Nc, -3) - 1 / Nc + 2 / (pow(Nc, 3) * pow(r, 2)) - 2 / (Nc * pow(r, 2)) - 2 / (pow(Nc, 3) * r) + 2 / (Nc * r), 2 / pow(Nc, 2) + 4 / (pow(Nc, 2) * pow(r, 2)) - 4 / (pow(Nc, 2) * r);
    return m;
};

MatrixXd f0_qQqQ(double r)
{
    Matrix2d m;

    m << 0.5 + 1 / (2. * pow(Nc, 4)) - pow(Nc, -2) + pow(r, -2) + 1 / (pow(Nc, 4) * pow(r, 2)) - 2 / (pow(Nc, 2) * pow(r, 2)) - 1 / r - 1 / (pow(Nc, 4) * r) + 2 / (pow(Nc, 2) * r), pow(Nc, -3) - 1 / Nc + 2 / (pow(Nc, 3) * pow(r, 2)) - 2 / (Nc * pow(r, 2)) - 2 / (pow(Nc, 3) * r) + 2 / (Nc * r),
        pow(Nc, -3) - 1 / Nc + 2 / (pow(Nc, 3) * pow(r, 2)) - 2 / (Nc * pow(r, 2)) - 2 / (pow(Nc, 3) * r) + 2 / (Nc * r), 2 / pow(Nc, 2) + 4 / (pow(Nc, 2) * pow(r, 2)) - 4 / (pow(Nc, 2) * r);
    return m;
};

MatrixXd f0_qQQq(double r)
{
    Matrix2d m;

    m << 0, 0,
        0, 2 + 4 / pow(1 - r, 2) - 4 / (1 - r);
    return m;
};

MatrixXd f0_qqbqqb(double r)
{
    Matrix2d m;

    m << 0.5 + 1 / (2. * pow(Nc, 4)) - pow(Nc, -2) + pow(r, -2) + 1 / (pow(Nc, 4) * pow(r, 2)) - 2 / (pow(Nc, 2) * pow(r, 2)) - 1 / r - 1 / (pow(Nc, 4) * r) + 2 / (pow(Nc, 2) * r), 2 + pow(Nc, -3) - 2 / pow(Nc, 2) - 1 / Nc + 2 / (pow(Nc, 3) * pow(r, 2)) - 2 / (Nc * pow(r, 2)) - 1 / r - 2 / (pow(Nc, 3) * r) + 1 / (pow(Nc, 2) * r) + 2 / (Nc * r) - r + r / pow(Nc, 2),
        2 + pow(Nc, -3) - 2 / pow(Nc, 2) - 1 / Nc + 2 / (pow(Nc, 3) * pow(r, 2)) - 2 / (Nc * pow(r, 2)) - 1 / r - 2 / (pow(Nc, 3) * r) + 1 / (pow(Nc, 2) * r) + 2 / (Nc * r) - r + r / pow(Nc, 2), 2 + 2 / pow(Nc, 2) - 8 / Nc + 4 / (pow(Nc, 2) * pow(r, 2)) - 4 / (pow(Nc, 2) * r) + 4 / (Nc * r) - 4 * r + (4 * r) / Nc + 4 * pow(r, 2);
    return m;
};

MatrixXd f0_qqqq(double r)
{
    Matrix2d m;

    m << 0.5 + (0.5 + pow(r, -2) - 1 / r) / pow(Nc, 4) + (-1 - 2 / pow(r, 2) + 2 / r) / pow(Nc, 2) + pow(r, -2) - 1 / r, (1 + 2 / pow(r, 2) - 2 / r) / pow(Nc, 3) + (-(1 / (1 - r)) - 1 / r) / pow(Nc, 2) + (-1 - 2 / pow(r, 2) + 2 / r) / Nc + 1 / (1 - r) + 1 / r,
        (1 + 2 / pow(r, 2) - 2 / r) / pow(Nc, 3) + (-(1 / (1 - r)) - 1 / r) / pow(Nc, 2) + (-1 - 2 / pow(r, 2) + 2 / r) / Nc + 1 / (1 - r) + 1 / r, 2 + (-4 / (1 - r) - 4 / r) / Nc + (2 + 4 / pow(r, 2) - 4 / r) / pow(Nc, 2) + 4 / pow(1 - r, 2) - 4 / (1 - r);
    return m;
};

MatrixXd f0_ggqqb(double r)
{
    Matrix3d m;
    m << (-1 + 1 / (2. * (1 - r)) + 1 / (2. * r)) / pow(Nc, 2), (-1 - 1 / (2. * (1 - r)) + 1 / (2. * r) + 2 * r) / Nc, (-1 + 1 / (2. * (1 - r)) + 1 / (2. * r)) / Nc,
        (-1 - 1 / (2. * (1 - r)) + 1 / (2. * r) + 2 * r) / Nc, -3 + 1 / (2. * (1 - r)) + 1 / (2. * r) + 4 * r - 4 * pow(r, 2), -1 - 1 / (2. * (1 - r)) + 1 / (2. * r) + 2 * r,
        (-1 + 1 / (2. * (1 - r)) + 1 / (2. * r)) / Nc, -1 - 1 / (2. * (1 - r)) + 1 / (2. * r) + 2 * r, -1 + 1 / (2. * (1 - r)) + 1 / (2. * r);
    return m;
};

MatrixXd f0_qgqg(double r)
{
    Matrix3d m;
    m << (0.5 + 1 / (2. * (1 - r)) - r / 2.) / pow(Nc, 2), (1.5 - 1 / (2. * (1 - r)) - 2 / r - r / 2.) / Nc, (0.5 + 1 / (2. * (1 - r)) - r / 2.) / Nc,
        (1.5 - 1 / (2. * (1 - r)) - 2 / r - r / 2.) / Nc, 2.5 + 1 / (2. * (1 - r)) + 4 / pow(r, 2) - 4 / r - r / 2., 1.5 - 1 / (2. * (1 - r)) - 2 / r - r / 2.,
        (0.5 + 1 / (2. * (1 - r)) - r / 2.) / Nc, 1.5 - 1 / (2. * (1 - r)) - 2 / r - r / 2., 0.5 + 1 / (2. * (1 - r)) - r / 2.;
    return m;
};

MatrixXd f0_qggq(double r)
{
    Matrix3d m;
    m << (1 / (2. * r) + r / 2.) / pow(Nc, 2), (1 - 2 / (1 - r) - 1 / (2. * r) + r / 2.) / Nc, (1 / (2. * r) + r / 2.) / Nc,
        (1 - 2 / (1 - r) - 1 / (2. * r) + r / 2.) / Nc, 2 + 4 / pow(1 - r, 2) - 4 / (1 - r) + 1 / (2. * r) + r / 2., 1 - 2 / (1 - r) - 1 / (2. * r) + r / 2.,
        (1 / (2. * r) + r / 2.) / Nc, 1 - 2 / (1 - r) - 1 / (2. * r) + r / 2., 1 / (2. * r) + r / 2.;
    return m;
};

MatrixXd f0_gggg(double r)
{
    Matrix9d m;

    m << 3 + pow(r, -2) - 2 / r - 2 * r + pow(r, 2), -2 + 1 / (1 - r) + 1 / r + r - pow(r, 2), -1 - 1 / (1 - r) - pow(r, -2) + 1 / r + r, -1 - 1 / (1 - r) - pow(r, -2) + 1 / r + r, -2 + 1 / (1 - r) + 1 / r + r - pow(r, 2), 3 + pow(r, -2) - 2 / r - 2 * r + pow(r, 2), 0, 0, 0,
        -2 + 1 / (1 - r) + 1 / r + r - pow(r, 2), 2 + pow(1 - r, -2) - 2 / (1 - r) + pow(r, 2), -pow(1 - r, -2) + 1 / (1 - r) - 1 / r - r, -pow(1 - r, -2) + 1 / (1 - r) - 1 / r - r, 2 + pow(1 - r, -2) - 2 / (1 - r) + pow(r, 2), -2 + 1 / (1 - r) + 1 / r + r - pow(r, 2), 0, 0, 0,
        -1 - 1 / (1 - r) - pow(r, -2) + 1 / r + r, -pow(1 - r, -2) + 1 / (1 - r) - 1 / r - r, 1 + pow(1 - r, -2) + pow(r, -2), 1 + pow(1 - r, -2) + pow(r, -2), -pow(1 - r, -2) + 1 / (1 - r) - 1 / r - r, -1 - 1 / (1 - r) - pow(r, -2) + 1 / r + r, 0, 0, 0,
        -1 - 1 / (1 - r) - pow(r, -2) + 1 / r + r, -pow(1 - r, -2) + 1 / (1 - r) - 1 / r - r, 1 + pow(1 - r, -2) + pow(r, -2), 1 + pow(1 - r, -2) + pow(r, -2), -pow(1 - r, -2) + 1 / (1 - r) - 1 / r - r, -1 - 1 / (1 - r) - pow(r, -2) + 1 / r + r, 0, 0, 0,
        -2 + 1 / (1 - r) + 1 / r + r - pow(r, 2), 2 + pow(1 - r, -2) - 2 / (1 - r) + pow(r, 2), -pow(1 - r, -2) + 1 / (1 - r) - 1 / r - r, -pow(1 - r, -2) + 1 / (1 - r) - 1 / r - r, 2 + pow(1 - r, -2) - 2 / (1 - r) + pow(r, 2), -2 + 1 / (1 - r) + 1 / r + r - pow(r, 2), 0, 0, 0,
        3 + pow(r, -2) - 2 / r - 2 * r + pow(r, 2), -2 + 1 / (1 - r) + 1 / r + r - pow(r, 2), -1 - 1 / (1 - r) - pow(r, -2) + 1 / r + r, -1 - 1 / (1 - r) - pow(r, -2) + 1 / r + r, -2 + 1 / (1 - r) + 1 / r + r - pow(r, 2), 3 + pow(r, -2) - 2 / r - 2 * r + pow(r, 2), 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0;
    return m;
};


