#include <cmath>

#include "../../../SpaceHub/src/spaceHub.hpp"
#include "../../../SpaceHub/src/taskflow/taskflow.hpp"
using namespace hub;
using namespace orbit;
using namespace callback;
using f = force::Interactions<force::NewtonianGrav>;
using Solver = methods::DefaultMethod<f, particles::SizeParticles>;
using Particle = Solver::Particle;
using namespace random;

size_t N = 20;
std::vector<double> f1s;
std::vector<double> f2s;
double M1 = 1;
double M2 = 1;
double mj = 1e-3;
void full(size_t i, size_t j, size_t tot_num) {
    double fv1 = f1s[i];
    double fv2 = f2s[j];

    double a2 = 10;
    double a1 = a2 * fv2;

    double v_orb = sqrt(consts::G * (M1 + mj) / a2);

    double v_inf = v_orb * fv1;

    double a_impulse = consts::G * (M1 + M2) / (v_inf * v_inf);

    double rp_max = 3 * a2;
    //double b_max = 3 * std::max(a_impulse, a2);
    double b_max  = rp_max * sqrt(1 + 2*consts::G*M1/(rp_max*v_inf*v_inf));

    double r_start = 10 * a2;

    std::ofstream file{"full-" + std::to_string(i) + "-" + std::to_string(j) + ".txt"};
    // file << std::scientific << std::setprecision(6);

    for (size_t k = 0; k < tot_num; ++k) {
        Particle star{M1, unit::Rs}, j1{mj, unit::Rj}, j2{mj, unit::Rj};

        Particle intruder{M2, unit::Rs};

        double nu1 = random::Uniform(0, 2 * consts::pi);

        double nu2 = random::Uniform(0, 2 * consts::pi);

        double I = random::Uniform(0, 2 * consts::pi);

        double Omega = random::Uniform(0, 2 * consts::pi);

        auto j1_orbit = Elliptic{star.mass, j1.mass, a1, 0.0, I, Omega, 0.0, nu1};

        auto j2_orbit = Elliptic{star.mass, j2.mass, a2, 0.0, I, Omega, 0.0, nu2};

        orbit::move_particles(j1_orbit, j1);

        orbit::move_particles(j2_orbit, j2);

        move_to_COM_frame(star, j1, j2);

        double b = sqrt(random::Uniform(unit::Rs*unit::Rs, b_max*b_max));

        double w = random::Uniform(0, 2 * consts::pi);
        
        auto orb =
            orbit::Hyperbolic(M_tot(star, j1, j2), intruder.mass, v_inf, b, w, 0.0, 0.0, r_start, orbit::Hyper::in);

        orbit::move_particles(orb, intruder);

        orbit::move_to_COM_frame(star, j1, j2, intruder);

        Solver sim{0, star, j1, j2, intruder};

        double t_end = 4 * orbit::time_to_periapsis(group(star, j1, j2), intruder);

        Solver::RunArgs args;

        args.add_stop_condition(t_end);

        int status = 0;

        auto collision = [&](auto& ptc, auto h) {
            size_t num = ptc.number();
            for (size_t i = 0; i < num; ++i) {
                for (size_t j = i + 1; j < num; ++j) {
                    if (norm(ptc.pos(i) - ptc.pos(j)) <= ptc.radius(i) + ptc.radius(j)) {
                        status = -1;
                        return true;
                    }
                }
            }
            return false;
        };
        args.atol = 1e-9;
        args.rtol = 1e-9;
        args.add_stop_condition(collision);

        args.add_stop_point_operation([&](auto& ptc, auto h) {
            if (status == -1) {
                file << status << ' ' << b_max << ' ' << b << ' ' << 0 << ' ' << 0 << ' ' << 0 << '\n';
                return;
            }
            auto [a10, e10] = calc_a_e(ptc.mass(0) + ptc.mass(1), ptc.pos(1) - ptc.pos(0), ptc.vel(1) - ptc.vel(0));
            auto [a20, e20] = calc_a_e(ptc.mass(0) + ptc.mass(2), ptc.pos(2) - ptc.pos(0), ptc.vel(2) - ptc.vel(0));
            auto [a13, e13] = calc_a_e(ptc.mass(3) + ptc.mass(1), ptc.pos(1) - ptc.pos(3), ptc.vel(1) - ptc.vel(3));
            auto [a23, e23] = calc_a_e(ptc.mass(3) + ptc.mass(2), ptc.pos(2) - ptc.pos(3), ptc.vel(2) - ptc.vel(3));
            auto [a12, e12] = calc_a_e(ptc.mass(1) + ptc.mass(2), ptc.pos(2) - ptc.pos(1), ptc.vel(2) - ptc.vel(1));
            auto vcm = norm((ptc.mass(1)*ptc.vel(1)+ptc.mass(2)*ptc.vel(2))/(ptc.mass(1)+ptc.mass(2))-ptc.vel(0));

            if (a10 > 0 && a20 > 0) {
                status = 0;
            } else if (a10 < 0 && a20 > 0 && a13 < 0 && a23 < 0) {
                status = 1;  // eject j1
            } else if (a10 > 0 && a20 < 0 && a13 < 0 && a23 < 0) {
                status = 2;  // eject j2
            } else if (a10 < 0 && a20 > 0 && a13 > 0 && a23 < 0) {
                status = 3;  // capture j1
            } else if (a10 > 0 && a20 < 0 && a13 < 0 && a23 > 0) {
                status = 4;  // capture j2
            } else if (a10 < 0 && a20 < 0 && a13 < 0 && a23 > 0) {
                status = 5;  // eject both, capture j2
            } else if (a10 < 0 && a20 < 0 && a13 > 0 && a23 < 0) {
                status = 6;  // eject both, capture j1
            } else if (a10 < 0 && a20 < 0 && a13 > 0 && a23 > 0) {
                status = 7;  // eject both capture both
            } else if (a10 < 0 && a20 < 0 && a13 < 0 && a23 < 0 && a12 < 0) {
                status = 8;  // eject both capture none, no binary
            } else if (a10 < 0 && a20 < 0 && a13 < 0 && a23 < 0 && a12 > 0) {
                status = 9;  // eject both capture none, binary!
            }

	    if(status != 9){
		vcm = 0;
	    }
            file << status << ' ' << b_max << ' ' << b << ' ' << a12 << ' ' << e12 << ' ' << vcm/v_orb << '\n';
        });

        sim.run(args);
    }
    std::cout << "Done " << i << ' ' << j << '\n';
}

int main() {
    double f1min = -1;
    double f1max = 1;
    double df1 = (f1max - f1min) / (N - 1);
    double f2min = 0.2;
    double f2max = 0.8;//sqrt(1 - 2 * pow(2 * mj / 3 / M1, 1.0 / 3));
    double df2 = (f2max - f2min) / (N - 1);

    for (size_t i = 0; i < N; ++i) {
        f1s.push_back(pow(10, f1min + i * df1));
        f2s.push_back(f2min + i * df2);
    }

    tf::Executor executor;  // create multi thread executor; thread number = cpu logical core number

    size_t tot_run = 10000000;
    for (int i = N - 1; i >= 0; i--) {
        for (int j = N - 1; j >= 0; j--) {
            executor.silent_async(full, i, j, tot_run);
        }
    }
    executor.wait_for_all();
    return 0;
}
