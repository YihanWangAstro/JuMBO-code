#include <cmath>

#include "SpaceHub/src/spaceHub.hpp"
#include "SpaceHub/src/taskflow/taskflow.hpp"
using namespace hub;
using namespace orbit;
using namespace callback;
using f = force::Interactions<force::NewtonianGrav>;
using Solver = methods::DefaultMethod<f, particles::SizeParticles>;
using Particle = Solver::Particle;
using namespace random;

double M1 = 1;
double M2 = 1;
double mj = 1e-3;
void fixed_I(size_t tot_num, int I_deg, int Omega_deg) {
    double fv1 = 0.1;

    double fv2 = 0.7;

    double a2 = 10;

    double a1 = a2 * fv2 * fv2;

    double v_orb = sqrt(consts::G * (M1 + mj) / a2);

    double v_inf = v_orb * fv1;

    double a_impulse = consts::G * (M1 + M2) / (v_inf * v_inf);

    double rp_max = 3 * a2;
    // double b_max = 3 * std::max(a_impulse, a2);
    double b_max = rp_max * sqrt(1 + 2 * consts::G * M1 / (rp_max * v_inf * v_inf));

    double r_start = 10 * a2;

    std::ofstream file{"fixed-" + std::to_string(I_deg) + "-" + std::to_string(Omega_deg) + ".txt"};

    for (size_t k = 0; k < tot_num; ++k) {
        Particle star{M1, unit::Rs}, j1{mj, unit::Rj}, j2{mj, unit::Rj};

        Particle intruder{M2, unit::Rs};

        double nu1 = random::Uniform(0, 2 * consts::pi);

        double nu2 = random::Uniform(0, 2 * consts::pi);

        double I = I_deg * consts::pi / 180;

        double Omega = Omega_deg * consts::pi / 180;

        auto j1_orbit = Elliptic{star.mass, j1.mass, a1, 0.0, I, Omega, 0.0, nu1};

        auto j2_orbit = Elliptic{star.mass, j2.mass, a2, 0.0, I, Omega, 0.0, nu2};

        orbit::move_particles(j1_orbit, j1);

        orbit::move_particles(j2_orbit, j2);

        move_to_COM_frame(star, j1, j2);

        double b = random::Uniform(0, b_max);

        double w = random::Uniform(0, 2 * consts::pi);

        auto orb =
            orbit::Hyperbolic(M_tot(star, j1, j2), intruder.mass, v_inf, b, w, 0.0, 0.0, r_start, orbit::Hyper::in);

        orbit::move_particles(orb, intruder);

        orbit::move_to_COM_frame(star, j1, j2, intruder);

        Solver sim{0, star, j1, j2, intruder};

        double t_end = 2 * orbit::time_to_periapsis(group(star, j1, j2), intruder);

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
                file << status << ' ' << b_max << ' ' << b << ' ' << 0 << ' ' << 0 << '\n';
                return;
            }
            auto [a10, e10] = calc_a_e(ptc.mass(0) + ptc.mass(1), ptc.pos(1) - ptc.pos(0), ptc.vel(1) - ptc.vel(0));
            auto [a20, e20] = calc_a_e(ptc.mass(0) + ptc.mass(2), ptc.pos(2) - ptc.pos(0), ptc.vel(2) - ptc.vel(0));
            auto [a13, e13] = calc_a_e(ptc.mass(3) + ptc.mass(1), ptc.pos(1) - ptc.pos(3), ptc.vel(1) - ptc.vel(3));
            auto [a23, e23] = calc_a_e(ptc.mass(3) + ptc.mass(2), ptc.pos(2) - ptc.pos(3), ptc.vel(2) - ptc.vel(3));
            auto [a12, e12] = calc_a_e(ptc.mass(1) + ptc.mass(2), ptc.pos(2) - ptc.pos(1), ptc.vel(2) - ptc.vel(1));

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
            file << status << ' ' << b_max << ' ' << b << ' ' << a12 << ' ' << e12 << '\n';
        });

        sim.run(args);
    }
}

int main() {
    tf::Executor executor;  // create multi thread executor; thread number = cpu logical core number

    size_t tot_run = 1000000;
    for (int i = 0; i <= 180; i += 15) {
        for (int j = 0; j <= 360; j += 15) {
            executor.silent_async(fixed_I, tot_run, i, j);
        }
    }
    executor.wait_for_all();
    return 0;
}
