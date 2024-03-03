#include <cmath>

#include "SpaceHub/src/spaceHub.hpp"
#include "SpaceHub/src/taskflow/taskflow.hpp"
using namespace hub;
using namespace orbit;
using namespace callback;
using f = force::Interactions<force::NewtonianGrav>;
using Solver = methods::Sym4<f, particles::SizeParticles>;
using Particle = Solver::Particle;
using namespace random;

double M1 = 1;
double M2 = 1;
double mj = 1e-3;

void traj(size_t id, double a2, double a1, double nu1, double nu2, double v_inf, double b, double w, double r_start) {
    Particle star{M1, unit::Rs}, j1{mj, unit::Rj}, j2{mj, unit::Rj};

    Particle intruder{M2, unit::Rs};

    double I = 0 * consts::pi / 180;

    double Omega = 0 * consts::pi / 180;

    auto j1_orbit = Elliptic{star.mass, j1.mass, a1, 0.0, I, Omega, 0.0, nu1};

    auto j2_orbit = Elliptic{star.mass, j2.mass, a2, 0.0, I, Omega, 0.0, nu2};

    orbit::move_particles(j1_orbit, j1);

    orbit::move_particles(j2_orbit, j2);

    move_to_COM_frame(star, j1, j2);

    auto orb = orbit::Hyperbolic(M_tot(star, j1, j2), intruder.mass, v_inf, b, w, 0.0, 0.0, r_start, orbit::Hyper::in);

    orbit::move_particles(orb, intruder);

    orbit::move_to_COM_frame(star, j1, j2, intruder);

    Solver sim{0, star, j1, j2, intruder};

    double t_end = 4 * orbit::time_to_periapsis(star, intruder);

    Solver::RunArgs args;

    args.add_stop_condition(t_end);

    args.atol = 1e-11;

    args.rtol = 1e-11;

    args.add_operation(DefaultWriter("trajectories" + std::to_string(id) + ".txt"));

    sim.run(args);
}
void fixed_I(size_t id) {
    double fv1 = 1;

    double fv2 = sqrt(0.7);

    double a2 = 1;

    double a1 = a2 * fv2 * fv2;

    double v_orb = sqrt(consts::G * (M1 + mj) / a2);

    double v_inf = v_orb * fv1;

    double a_impulse = consts::G * (M1 + M2) / (v_inf * v_inf);

    double rp_max = 3 * a2;
    // double b_max = 3 * std::max(a_impulse, a2);
    double b_max = rp_max * sqrt(1 + 2 * consts::G * M1 / (rp_max * v_inf * v_inf));

    double r_start = 20 * a2;

    bool found_one = false;
    for (size_t k = 0; !found_one; ++k) {
        Particle star{M1, unit::Rs}, j1{mj, unit::Rj}, j2{mj, unit::Rj};

        Particle intruder{M2, unit::Rs};

        double nu1 = random::Uniform(0, 2 * consts::pi);

        double nu2 = random::Uniform(0, 2 * consts::pi);

        double I = 0 * consts::pi / 180;

        double Omega = 0 * consts::pi / 180;

        auto j1_orbit = Elliptic{star.mass, j1.mass, a1, 0.0, I, Omega, 0.0, nu1};

        auto j2_orbit = Elliptic{star.mass, j2.mass, a2, 0.0, I, Omega, 0.0, nu2};

        orbit::move_particles(j1_orbit, j1);

        orbit::move_particles(j2_orbit, j2);

        move_to_COM_frame(star, j1, j2);

        double b = random::Uniform(0, b_max);

        double w = 0;
        // random::Uniform(0, 2 * consts::pi);

        auto orb =
            orbit::Hyperbolic(M_tot(star, j1, j2), intruder.mass, v_inf, b, w, 0.0, 0.0, r_start, orbit::Hyper::in);

        orbit::move_particles(orb, intruder);

        orbit::move_to_COM_frame(star, j1, j2, intruder);

        Solver sim{0, star, j1, j2, intruder};

        double t_end = 2 * orbit::time_to_periapsis(star, intruder);

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
        args.atol = 1e-11;
        args.rtol = 1e-11;
        args.add_stop_condition(collision);

        args.add_stop_point_operation([&](auto& ptc, auto h) {
            if (status == -1) {
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
                /*if (e23 < 0.5 && e13 < 0.5) {
                    traj(id, a2, a1, nu1, nu2, v_inf, b, w, r_start);
                    found_one = true;
                }*/
            } else if (a10 < 0 && a20 < 0 && a13 < 0 && a23 < 0 && a12 < 0) {
                status = 8;  // eject both capture none, no binary
            } else if (a10 < 0 && a20 < 0 && a13 < 0 && a23 < 0 && a12 > 0) {
                status = 9;  // eject both capture none, binary!
                             // if (e12 < 0.5) {
                traj(id, a2, a1, nu1, nu2, v_inf, b, w, r_start);
                found_one = true;
                //}
            }
        });

        sim.run(args);
    }
}

int main() {
    tf::Executor executor;  // create multi thread executor; thread number = cpu logical core number

    executor.silent_async(fixed_I, 0);

    executor.wait_for_all();
    return 0;
}
