#include <cmath>
#include <cstdio>
#include <fstream>

#include "../../SpaceHub/src/taskflow/taskflow.hpp"
FILE *outfile;
void full(size_t i, size_t j) {
    double a2 = 10;
    std::ifstream file{"fixed-" + std::to_string(i) + "-" + std::to_string(j) + ".txt"};
    double a_tot = 0;
    double e_tot = 0;
    int n[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    int stat = 0;
    double bmax = 0;
    double b = 0;
    double a = 0;
    double e = 0;

    int N = 0;

    for (; file; ++N) {
        file >> stat >> bmax >> b >> a >> e;
        n[stat + 1] += 1;
        if (stat == 9) {
            a_tot += log10(a / a2);
            e_tot += e;
        }
    }
    double a_ave = a_tot / n[10];
    double e_ave = e_tot / n[10];
    file.close();

    file.open("fixed-" + std::to_string(i) + "-" + std::to_string(j) + ".txt");

    double a_std = 0;
    double e_std = 0;
    for (size_t k = 0; file; ++k) {
        file >> stat >> bmax >> b >> a >> e;
        if (stat == 9) {
            a_std += (log10(a / a2) - a_ave) * (log10(a / a2) - a_ave);
            e_std += (e - e_ave) * (e - e_ave);
        }
    }
    a_std /= (n[10] - 1);
    a_std = sqrt(a_std);
    e_std /= (n[10] - 1);
    e_std = sqrt(e_std);

    fprintf(outfile, "%zu %zu %lf %lf %lf %lf %lf %d %d %d %d %d %d %d %d %d %d %d %d\n", i, j, bmax, a_ave, a_std,
            e_ave, e_std, N, n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10]);
    std::cout << i << ' ' << j << '\n';
}

int main() {
    //int N = 20;
    tf::Executor executor;
  //  for(int rn = 2; rn<=10;++rn){
    //char outname[100];
    //sprintf(outname,"stat%d.txt",rn);
    outfile = fopen("angle-stat.txt", "w");
    for (int i = 0; i <= 180; i+=9) {
        for (int j = 0; j <= 180; j+=9) {
            executor.silent_async(full, i, j);
        }
    }
    executor.wait_for_all();
    //fclose(outfile);
    //}
    return 0;
}
