// MV dipole model with simple matching for x > x0
// and simple evolution in x to compensate rcBK

#ifndef DIPOLE_MODELS_H
#define DIPOLE_MODELS_H

#include <cmath>

class MVDipole {
public:
    static constexpr double x0 = 0.01;
    static constexpr double Lam = 0.241;
    static constexpr double Qs02_ = 0.104;
    static constexpr double Qs02_e = 0.060;
    static constexpr double ec = 18.9;
    static constexpr double Qs02_g = 0.169;
    static constexpr double gamma = 1.119;

    static double Qs2(double x, double Qs02) {
        return Qs02 * std::pow(x0 / x, 0.3);
    }

    // MV model
    static double MV(double r, double x) {
        if (x >= 1.0) {
            return 0.0;
        } else if (x >= x0) {
            double result = std::exp(-r * r * Qs02_ / 4.0 * std::log(1.0 / (Lam * r) + std::exp(1.0)));
            result *= std::pow((1.0 - x) / (1.0 - x0), 4.0) * std::pow(x0 / x, 0.15);
            return result;
        } else {
            double result = std::exp(-r * r * Qs2(x, Qs02_) / 4.0 * std::log(1.0 / (Lam * r) + std::exp(1.0)));
            return result;
        }
    }


    // MV_gamma model
    static double MVg(double r, double x) {
        if (x >= 1.0) {
            return 0.0;
        } else if (x >= x0) {
            double result = std::exp(-std::pow(r * r * Qs02_g, gamma) / 4.0 * std::log(1.0 / (Lam * r) + std::exp(1.0)));
            result *= std::pow((1.0 - x) / (1.0 - x0), 4.0) * std::pow(x0 / x, 0.15);
            return result;
        } else {
            double result = std::exp(-std::pow(r * r * Qs2(x, Qs02_g), gamma) / 4.0 * std::log(1.0 / (Lam * r) + std::exp(1.0)));
            return result;
        }
    }

    // MV_e model
    static double MVe(double r, double x) {
        if (x >= 1.0) {
            return 0.0;
        } else if (x >= x0) {
            double result = std::exp(-r * r * Qs02_e / 4.0 * std::log(1.0 / (Lam * r) + ec * std::exp(1.0)));
            result *= std::pow((1.0 - x) / (1.0 - x0), 4.0) * std::pow(x0 / x, 0.15);
            return result;
        } else {
            double result = std::exp(-r * r * Qs2(x, Qs02_e) / 4.0 * std::log(1.0 / (Lam * r) + ec * std::exp(1.0)));
            return result;
        }
    }


    static double rlaplaceN(double r, double x) {
        if (x >= 1.0 || r == 0.0) {
            return 0.0;
        } else if (x < x0) {
            double Q = std::sqrt(Qs2(x, Qs02_));
            double result = -(Q * Q * std::pow(1.0 / (Lam * r) + std::exp(1.0), -1.0 / 4.0 * Q * Q * r * r)
                * (-4.0 * (4.0 * std::pow(std::exp(1.0), 2.0) * Lam * Lam * r * r
                    + std::exp(1.0) * Lam * r * (Q * Q * r * r + 8.0) + Q * Q * r * r + 4.0) * std::log(1.0 / (Lam * r) + std::exp(1.0))
                    + 4.0 * Q * Q * r * r * std::pow(std::exp(1.0) * Lam * r + 1.0, 2.0) * std::pow(std::log(1.0 / (Lam * r) + std::exp(1.0)), 2.0)
                    + 12.0 * std::exp(1.0) * Lam * r + Q * Q * r * r + 16.0)) / (16.0 * std::pow(std::exp(1.0) * Lam * r + 1.0, 2.0));
            return result;
        } else {
            double Q = std::sqrt(Qs02_);
            double result = -(Q * Q * std::pow(1.0 / (Lam * r) + std::exp(1.0), -1.0 / 4.0 * Q * Q * r * r)
                * (-4.0 * (4.0 * std::pow(std::exp(1.0), 2.0) * Lam * Lam * r * r
                    + std::exp(1.0) * Lam * r * (Q * Q * r * r + 8.0) + Q * Q * r * r + 4.0) * std::log(1.0 / (Lam * r) + std::exp(1.0))
                    + 4.0 * Q * Q * r * r * std::pow(std::exp(1.0) * Lam * r + 1.0, 2.0) * std::pow(std::log(1.0 / (Lam * r) + std::exp(1.0)), 2.0)
                    + 12.0 * std::exp(1.0) * Lam * r + Q * Q * r * r + 16.0)) / (16.0 * std::pow(std::exp(1.0) * Lam * r + 1.0, 2.0));
            result *= std::pow((1.0 - x) / (1.0 - x0), 4.0) * std::pow(x0 / x, 0.15);
            return result;
        }
    }

    static double rlaplaceNe(double r, double x) {
        if (x >= 1.0 || r == 0.0) {
            return 0.0;
        } else if (x < x0) {
            double Q = std::sqrt(Qs2(x, Qs02_e));
            double result = -(Q * Q * std::pow(1.0 / (Lam * r) + ec * std::exp(1.0), -1.0 / 4.0 * Q * Q * r * r)
                * (-4.0 * (4.0 * std::pow(ec * std::exp(1.0), 2.0) * Lam * Lam * r * r
                    + ec * std::exp(1.0) * Lam * r * (Q * Q * r * r + 8.0) + Q * Q * r * r + 4.0) * std::log(1.0 / (Lam * r) + ec * std::exp(1.0))
                    + 4.0 * Q * Q * r * r * std::pow(ec * std::exp(1.0) * Lam * r + 1.0, 2.0) * std::pow(std::log(1.0 / (Lam * r) + ec * std::exp(1.0)), 2.0)
                    + 12.0 * ec * std::exp(1.0) * Lam * r + Q * Q * r * r + 16.0)) / (16.0 * std::pow(ec * std::exp(1.0) * Lam * r + 1.0, 2.0));
            return result;
        } else {
            double Q = std::sqrt(Qs02_e);
            double result = -(Q * Q * std::pow(1.0 / (Lam * r) + ec * std::exp(1.0), -1.0 / 4.0 * Q * Q * r * r)
                * (-4.0 * (4.0 * std::pow(ec * std::exp(1.0), 2.0) * Lam * Lam * r * r
                    + ec * std::exp(1.0) * Lam * r * (Q * Q * r * r + 8.0) + Q * Q * r * r + 4.0) * std::log(1.0 / (Lam * r) + ec * std::exp(1.0))
                    + 4.0 * Q * Q * r * r * std::pow(ec * std::exp(1.0) * Lam * r + 1.0, 2.0) * std::pow(std::log(1.0 / (Lam * r) + ec * std::exp(1.0)), 2.0)
                    + 12.0 * ec * std::exp(1.0) * Lam * r + Q * Q * r * r + 16.0)) / (16.0 * std::pow(ec * std::exp(1.0) * Lam * r + 1.0, 2.0));
            result *= std::pow((1.0 - x) / (1.0 - x0), 4.0) * std::pow(x0 / x, 0.15);
            return result;
        }
    }

    static double rlaplaceNg(double r, double x) {
        if (x >= 1.0 || r == 0.0) {
            return 0.0;
        } else if (x < x0) {
            double Q = std::sqrt(Qs2(x, Qs02_g));
            double result = -((std::pow(Q * r, 2.0 * gamma)) * std::pow(1.0 / (Lam * r) + std::exp(1.0), -1.0 / 4.0 * std::pow(Q * r, 2.0 * gamma))
                * (4.0 * gamma * (std::exp(1.0) * Lam * r + 1.0) * std::log(1.0 / (Lam * r) + std::exp(1.0))
                    * ((std::pow(Q * r, 2.0 * gamma) * ((std::exp(1.0) * gamma * Lam * r + gamma)
                        * std::log(1.0 / (Lam * r) + std::exp(1.0)) - 1.0) - 2.0 * (1.0 + 2.0 * gamma)
                        * (std::exp(1.0) * Lam * r + 1.0)) + std::pow(Q * r, 2.0 * gamma) + 16.0 * gamma
                        * (std::exp(1.0) * Lam * r + 1.0) + 4.0))) / (16.0 * r * r * std::pow(std::exp(1.0) * Lam * r + 1.0, 2.0));
            return result;
        } else {
            double Q = std::sqrt(Qs02_g);
            double result = -((std::pow(Q * r, 2.0 * gamma)) * std::pow(1.0 / (Lam * r) + std::exp(1.0), -1.0 / 4.0 * std::pow(Q * r, 2.0 * gamma))
                * (4.0 * gamma * (std::exp(1.0) * Lam * r + 1.0) * std::log(1.0 / (Lam * r) + std::exp(1.0))
                    * ((std::pow(Q * r, 2.0 * gamma) * ((std::exp(1.0) * gamma * Lam * r + gamma)
                        * std::log(1.0 / (Lam * r) + std::exp(1.0)) - 1.0) - 2.0 * (1.0 + 2.0 * gamma)
                        * (std::exp(1.0) * Lam * r + 1.0)) + std::pow(Q * r, 2.0 * gamma) + 16.0 * gamma
                        * (std::exp(1.0) * Lam * r + 1.0) + 4.0))) / (16.0 * r * r * std::pow(std::exp(1.0) * Lam * r + 1.0, 2.0));
            result *= std::pow((1.0 - x) / (1.0 - x0), 4.0) * std::pow(x0 / x, 0.15);
            return result;
        }
    }
};

#endif // DIPOLE_MODELS_H