#include <cmath>

#include "SpaceHub/src/spaceHub.hpp"
#include "SpaceHub/src/taskflow/taskflow.hpp"
using namespace hub;
using namespace unit;
using namespace orbit;
using namespace callback;
using f = force::Interactions<force::NewtonianGrav>;
using Solver = methods::DefaultMethod<f, particles::SizeParticles>;
using Particle = Solver::Particle;
using namespace random;

void job_run(size_t i, size_t tot_num) {
    double v_inf = 2_kms;

    double b_max = 50_AU;

    double r_start = 100_AU;

    std::ofstream file("results/out-" + std::to_string(i) + ".txt");

    for (size_t i = 0; i < tot_num; ++i) {
        Particle star{Ms, Rs}, j1{Mj, Rj}, j2{Mj, Rj};

        Particle intruder{Ms, Rs};

        double a1 = 10_AU;

        double a2 = random::Uniform(11_AU, 20_AU);

        double nu1 = random::Uniform(0, 2 * consts::pi);

        double nu2 = random::Uniform(0, 2 * consts::pi);

        auto j1_orbit = Elliptic{star.mass, j1.mass, a1, 0.0, 0.0, 0.0, 0.0, nu1};

        auto j2_orbit = Elliptic{star.mass, j2.mass, a2, 0.0, 0.0, 0.0, 0.0, nu2};

        orbit::move_particles(j1_orbit, j1);

        orbit::move_particles(j2_orbit, j2);

        move_to_COM_frame(star, j1, j2);

        double b = random::Uniform(-b_max, b_max);

        double w = 0;

        if (b < 0) {
            w = consts::pi;
        }

        auto orb =
            orbit::Hyperbolic(M_tot(star, j1, j2), intruder.mass, v_inf, b, w, 0.0, 0.0, r_start, orbit::Hyper::in);

        orbit::move_particles(orb, intruder);

        orbit::move_to_COM_frame(star, j1, j2, intruder);

        Solver sim{0, star, j1, j2, intruder};

        double t_end = 5 * orbit::time_to_periapsis(group(star, j1, j2), intruder);

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
        args.atol = 1e-12;
        args.rtol = 1e-12;
        args.add_stop_condition(collision);

        args.add_stop_point_operation([&](auto& ptc, auto h) {
            if (status == -1) {
                file << status << ' ' << a2 << ' ' << nu1 << ' ' << nu2 << ' ' << b << ' ' << 0 << ' ' << 0 << '\n';
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
                // std::cout << "b =" << b << ", a2 = " << a2 << ", nu1 = " << nu1 << ", nu2= " << nu2 << ", a_p = " <<
                // a12
                //           << ", e_p=" << e12 << " floating binary!!\n";
            }
            file << status << ' ' << a2 << ' ' << nu1 << ' ' << nu2 << ' ' << b << ' ' << a12 << ' ' << e12 << '\n';
        });

        sim.run(args);
    }
}

int main() {
    tf::Executor executor;  // create multi thread executor; thread number = cpu logical core number
    size_t thread_num = 4;
    size_t tot_run = 10000000;
    for (size_t i = 0; i < thread_num; i++) {
        executor.silent_async(job_run, i, tot_run / thread_num);
    }
    executor.wait_for_all();
    return 0;
}
