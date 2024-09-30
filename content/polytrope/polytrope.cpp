#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>

struct State {
    double theta;
    double dtheta_dxi;
};

class Polytrope {

public:

    std::vector<double> xi;
    std::vector<State> sol;

    double n;

    Polytrope(const double _n=1.5) {
        assert(_n >= 0);
        n = _n;
    }

    State rhs(const double x, const State& s) {

        State f{};

        f.theta = s.dtheta_dxi;

        if (x == 0.0) {
            f.dtheta_dxi = (2.0/3.0) - std::pow(s.theta, n);
        } else {
            f.dtheta_dxi = -2.0 * s.dtheta_dxi / x - std::pow(s.theta, n);
        }

        return f;
    }

    int npts() {
        return xi.size();
    }

    void integrate(const double h0=1.e-2, const double tol=1.e-12) {

        // initial conditions

        State y{1.0, 0.0};

        // storage for the stages

        State tmp{};

        auto h = h0;

        double _xi = 0.0;

        while (h > tol) {

            // 4th order RK integration

            auto k1 = rhs(_xi, y);

            tmp.theta = y.theta + 0.5 * h * k1.theta;
            tmp.dtheta_dxi = y.dtheta_dxi + 0.5 * h * k1.dtheta_dxi;

            auto k2 = rhs(_xi + 0.5*h, tmp);

            tmp.theta = y.theta + 0.5 * h * k2.theta;
            tmp.dtheta_dxi = y.dtheta_dxi + 0.5 * h * k2.dtheta_dxi;

            auto k3 = rhs(_xi + 0.5*h, tmp);

            tmp.theta = y.theta + h * k3.theta;
            tmp.dtheta_dxi = y.dtheta_dxi + h * k3.dtheta_dxi;

            auto k4 = rhs(_xi + h, tmp);

            y.theta += (1./6.) * h * (k1.theta + 2.0*k2.theta + 2.0*k3.theta + k4.theta);
            y.dtheta_dxi += (1./6.) * h * (k1.dtheta_dxi + 2.0*k2.dtheta_dxi + 2.0*k3.dtheta_dxi + k4.dtheta_dxi);


            _xi += h;

            // set the new stepsize--our systems is always convex
            // (theta'' < 0), so the intersection of theta' with the
            // x-axis will always be a conservative estimate of the
            // radius of the star.  Make sure that the stepsize does
            // not take us past that.

            double R_est = _xi - y.theta/y.dtheta_dxi;

            if (_xi + h > R_est) {
                h = -y.theta/y.dtheta_dxi;
            }

            // store the solution:

            xi.push_back(_xi);
            sol.push_back(y);

        }
    }

    double get_xi1() {
        if (!xi.empty()) {
            return xi.back();
        } else {
            return -1;
        }
    }

    double get_minus_xisq_dtheta_dxi() {
        if (!xi.empty() && !sol.empty()) {
            auto xi_last = xi.back();
            auto sol_last = sol.back();

            return -xi_last * xi_last * sol_last.dtheta_dxi;
        } else {
            return -1;
        }
    }
};


int main() {

    for (int imodel = 0; imodel <= 9; ++imodel) {

        double n = static_cast<double>(imodel) / 2.0;

        Polytrope p(n);
        p.integrate();

        if (imodel == 0) {
            std::cout << "#" << std::setw(4) << "n" << std::setw(22) << "xi_1" << std::setw(22) << "-xi**2 dtheta/dxi |_xi1" << std::endl;
        }
        auto xi1 = p.get_xi1();
        auto minus_xisq_dtheta_dxi = p.get_minus_xisq_dtheta_dxi();

        std::cout << std::setw(4) << n << std::setw(22) << xi1 << std::setw(22) << minus_xisq_dtheta_dxi << std::endl;
    }
}
