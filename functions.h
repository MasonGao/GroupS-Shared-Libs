
double Beta[4] = {(11 * CA) / 3. - (4 * nf * TF) / 3., (34 * pow(CA, 2)) / 3. - (20 * CA * nf * TF) / 3. - 4 * CF *nf *TF, (2857 * pow(CA, 3)) / 54. + ((-1415 * pow(CA, 2)) / 27. - (205 * CA * CF) / 9. + 2 * pow(CF, 2)) * nf *TF + ((158 * CA) / 27. + (44 * CF) / 9.) * pow(nf, 2) * pow(TF, 2), 149753 / 6. + (1093 * pow(nf, 3)) / 729. + 3564 * riemann_zeta(3) + pow(nf, 2) * (50065 / 162. + (6472 * riemann_zeta(3)) / 81.) - nf * (1078361 / 162. + (6508 * riemann_zeta(3)) / 27.)};

double Gamma_cusp[3] = {4, 4 * (CA * (67 / 9. - pow(PI, 2) / 3.) - (20 * nf * TF) / 9.), 4 * ((-16 * pow(nf, 2) * pow(TF, 2)) / 27. + CA * nf * TF * (-418 / 27. + (40 * pow(PI, 2)) / 27. - (56 * riemann_zeta(3)) / 3.) + pow(CA, 2) * (245 / 6. - (134 * pow(PI, 2)) / 27. + (11 * pow(PI, 4)) / 45. + (22 * riemann_zeta(3)) / 3.) + CF * nf * TF * (-55 / 3. + 16 * riemann_zeta(3)))};

double gamma_q[2] = {-3 * CF, CF *nf * (130 / 27. + (2 * pow(PI, 2)) / 3.) * TF + pow(CF, 2) * (-1.5 + 2 * pow(PI, 2) - 24 * riemann_zeta(3)) + CA *CF * (-961 / 54. - (11 * pow(PI, 2)) / 6. + 26 * riemann_zeta(3))};

double gamma_g[2] = {(-11 * CA) / 3. + (4 * nf * TF) / 3., 4 * CF *nf *TF + CA *nf * (256 / 27. - (2 * pow(PI, 2)) / 9.) * TF + pow(CA, 2) * (-692 / 27. + (11 * pow(PI, 2)) / 18. + 2 * riemann_zeta(3))};


double S_f(int o, double ga[], double mu, double nu)
{
    double r = r_as(nu, mu);
    if (o == 0)
        return (ga[0] / 4. / pow(Beta[0], 2) * (4 * PI / a_s(nu) * (1 - 1 / r - log(r))));
    if (o == 1)
        return (ga[0] / 4. / pow(Beta[0], 2) * (4 * PI / a_s(nu) * (1 - 1 / r - log(r))) + ga[0] / 4 / pow(Beta[0], 2) * ((ga[1] / ga[0] - Beta[1] / Beta[0]) * (1 - r + log(r)) + Beta[1] / 2 / Beta[0] * pow(log(r), 2)));
};

double a_f(int o, double ga[], double f, double i)
{
    if (o == -1)
        return 0.;
    if (o == 0)
        return (ga[0] / 2. / Beta[0] * (log(r_as(f, i))));
    if (o == 1)
        return (ga[0] / 2. / Beta[0] * (log(r_as(f, i))) + ga[0] / 2. / Beta[0] * (ga[1] / ga[0] - Beta[1] / Beta[0]) * (a_s(i) - a_s(f)) / 4 / PI);
};


double x1x2s(double yc, double yd, double pt)
{
    return ((exp(-yc) + exp(-yd)) * (exp(yc) + exp(yd)) * pow(pt, 2));
};

double UE(int o, int a, int b, int c, int d, double s, double mu, double nu)
{
    int ng = 0;
    (a == 0) ? ng++ : ng = ng;
    (b == 0) ? ng++ : ng = ng;
    (c == 0) ? ng++ : ng = ng;
    (d == 0) ? ng++ : ng = ng;

    double res_UE = pow((nu * nu / s), ((-1) * ((4 - ng) * CF + ng * CA) * a_f(o, Gamma_cusp, mu, nu))) * exp(2 * (((4 - ng) * CF + ng * CA) * S_f(o, Gamma_cusp, mu, nu) + (4 - ng) * a_f(o - 1, gamma_q, mu, nu) + ng * a_f(o - 1, gamma_g, mu, nu)));

    return res_UE;
};

double U_J(int o, int a, int b, int c, int d, double s, double mu, double nu)
{
    int ng = 0;
    (c == 0) ? ng++ : ng = ng;
    (d == 0) ? ng++ : ng = ng;

    double res_UJ = pow((nu * nu / s), (((2 - ng) * CF + ng * CA) * a_f(o, Gamma_cusp, mu, nu))) * exp((-2) * (((2 - ng) * CF + ng * CA) * S_f(o, Gamma_cusp, mu, nu) + (2 - ng) * a_f(o - 1, gamma_q, mu, nu) + ng * a_f(o - 1, gamma_g, mu, nu)));
    return res_UJ;
};

MatrixXcd u_f(MatrixXcd m, double mu, double nu)
{
    ComplexEigenSolver<MatrixXcd> ces(m);
    MatrixXcd v_m = ces.eigenvectors();
    double r_n = r_as(mu, nu);
    MatrixXcd g_D = v_m.inverse() * m * v_m / 2. / Beta[0];

    return (v_m * ((log(r_n) * g_D).exp()) * v_m.inverse());
};

double b_star(double b)
{
    return (b / sqrt(1 + b * b / bmax / bmax));
};

double mub_f(double b) //usually use b_star above as the b-variable
{
    return (2 / b / exp(EULERGAMMA));
};

double U_NG(int process, double mubs, double muj)
{
    double a = 0.85 * CA;
    double b = 0.86 * CA;
    double c = 1.33;
    double u = log(a_s(mubs) / a_s(muj)) / Beta[0];
    double Ccd = CA;
    if (process == 0)
    {
        Ccd = CA;
    }
    else
    {
        Ccd = CF;
    }
    return (exp((-1) * Ccd * CA * PI * PI / 3 * u * u * (1 + pow(a * u, 2)) / (1 + pow(b * u, c))));
};

double SNP(int process, int pa, double mu, double b)
{
    double Cab = CA;
    if (process == 0)
    {
        Cab = CA;
    }
    else
    {
        Cab = CF;
    }
    double g1pA = 0.106 * b * b;
    if (pa == 0) // proton
    {
        g1pA = 0.106 * b * b;
    }
    else if (pa == 1) // Pb
    {
        g1pA = (0.106 + aN * (pow(208., 1. / 3.) - 1.)) * b * b;
    }

    return (g1pA + Cab / CF * 0.84 / 2 * log(mu / sqrt(2.4)) * log(sqrt(1 + b * b / bmax / bmax)));
};
