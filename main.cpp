#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <random>

struct result
{
    double h;
    double dist;
    std::vector<double> al;
    double w;
    double d;
};

std::vector<double> alpha(int r) { //вектор шумов
    std::vector<double> al(r);
    std::mt19937 gen(time(0));
    std::uniform_real_distribution<double> urd(0, 1);
    double A = urd(gen);
    al[(r - 1) / 2] = A;

    for (auto i = (r - 1) / 2 - 1; i > 0; i--)
    {
        double sum = 0;
        for (auto j = i; j < r - i - 1; j++)
            sum += al[j];
        std::uniform_real_distribution<double> urd1(0, 1 - sum);
        double B = urd1(gen);
        al[i] = 0.5 * B;
    }
    for (auto i = (r - 1) / 2 + 1; i < r - 1; i++)
        al[i] = al[r - i - 1];
    double sum = 0;
    for (auto i = 1; i < r - 1; i++)
        sum += al[i];
    al[0] = 0.5 * (1 - sum);
    al[r - 1] = al[0];
    return al;
}

std::vector<double> filter(std::vector<double> F, std::vector<double> al, int r){
    std::vector<double> F_res(101);
    int M = (r - 1) / 2;
    for (auto i = 0; i <= 100; i++){
        double pr = 1;
        for (auto j = i - M; j <= i + M; j++){
            try {
                F.at(j);
                pr *= pow(F[j], al[j + M - i]);
            }
            catch (const std::out_of_range& s) {
                pr *= 1;
            }
        }
        F_res[i] = pr;
    }
    return F_res;
}

double w(std::vector<double> F)
{
    double max = 0;
    double mod;
    for (auto i = 1; i <= 100; i++) {
        mod = abs((F[i] - F[i - 1]));
        if (max < mod) max = mod;
    }
    return max;
}

double d(std::vector<double> F, std::vector<double> F_noise)
{
    double max = 0;
    for (auto i = 1; i <= 100; i++) {
        if (abs(F[i] - F_noise[i - 1]) > max) max = abs(F[i] - F_noise[i]);
    }
    return max;
}

double dist( double w, double d)
{
    if (w > d) return w;
    else return d;
}

double J(double h, double w, double d)
{
    return h * w + (1 - h) * d;
}

double n(double P, double E,  double X_max,  double X_min)
{
    return (log(1 - P) / log(1 - E / (X_max - X_min)));
}

int main() {
    int r3 = 3, r5 = 5;
    double X_min = 0, X_max = M_PI;
    int K = 100;
    double a1 = -0.25, a2 = 0.25;
    int L = 10;
    double H = 0;
    double P = 0.95;
    double E = 0.01;
    double W, D, minW, minD;
    double Dist, minDist = 100000;
    std::vector<double> minAL;
    std::vector<result> res;
    std::vector<double> al;
    double N = n(P, E, X_max, X_min);

    std::vector<double> X(101);
    std::vector<double> F(101);
    std::vector<double> F_noise(101);
    for (auto i = 0; i < K; i++) X[i] = X_min + i * (X_max - X_min) / K;
    for (auto i = 0; i <= K; i++) F[i] = sin(X[i]) + 0.5;
    for (auto i = 0; i <= K; i++)
    {
        srand(time(NULL));
        double A__ = rand() % 50;
        double A = a1 + (A__ / 100);
        F_noise[i] = F[i] + A;
    }
    std::cout << "without noise F[i]|with noise F[i]"<<std::endl;
    for (auto i = 0; i <= K; i++) {
        std::cout << std::setw(18) << F[i] << "|" << F_noise[i]<<std::endl;
    }

    std::vector<double> F_r3;
    std::cout <<std::endl<< " h|   dist|            alpha|             w|       d" << std::endl;
    for (auto j = 0; j <= L; j++){
        for (auto i = 0; i < K / 2; i++){
            al = alpha(r3);
            F_r3 = filter(F_noise, al, r3);
            W = w(F_r3);
            D = d(F_r3, F_noise);
            Dist = dist(W, D);

            if (Dist < minDist){
                minDist = Dist;
                minW = W;
                minD = D;
                minAL = al;
            }
        }
        H = (double)j / L;
        res.push_back({ H, minDist, minAL, minW, minD });
        std::cout << std::fixed << std::setprecision(1) << H << " |" << std::setprecision(4) << minDist << " | [ ";
        for (auto i : minAL)
            std::cout << i << " ";
        std::cout << "] |" << minW << "  |" << minD << std::endl;
        minDist = 100000;
    }

    std::sort(res.begin(), res.end(), [](auto x, auto y) {return x.dist < y.dist; });
    std::cout <<std::endl<< "h* = "<< res[0].h << "  J = " << J(res[0].h, res[0].w, res[0].d) << "  w = " << res[0].w << "  d = " << res[0].d << std::endl <<std::endl;
    F_r3 = filter(F_noise, minAL, r3);
    for (auto i : F_r3) std::cout << i << "; ";
    std::cout << std::endl;
    res.clear();
    minAL.clear();
    minD = 0;
    minW = 0;
    //////////////////////////////////////////////////////////////////////////
    std::vector<double> F_r5;
    std::cout << std::endl << " h|   dist|            alpha|                          w|       d" << std::endl;
    for (auto j = 0; j <= L; j++){
        for (auto i = 0; i < K / 2; i++){
            al = alpha(r5);
            F_r5 = filter(F_noise, al, r5);
            W = w(F_r5);
            D = d(F_r5, F_noise);
            Dist = dist(W, D);

            if (Dist < minDist){
                minDist = Dist;
                minW = W;
                minD = D;
                minAL = al;
            }
        }
        H = (double)j / L;
        res.push_back({ H, minDist, minAL, minW, minD });
        std::cout << std::fixed << std::setprecision(1) << H << " |" << std::setprecision(4) << minDist << " | [ ";
        for (auto i : minAL) std::cout << i << " ";
        std::cout << "] |" << minW << "  |" << minD << std::endl;
        minDist = 100000;
    }

    std::sort(res.begin(), res.end(), [](auto x, auto y) {return x.dist < y.dist; });
    std::cout << std::endl << "h* = " << res[0].h << "  J = " << J(res[0].h, res[0].w, res[0].d) << "  w = " << res[0].w << "  d = " << res[0].d << std::endl << std::endl;
    F_r5 = filter(F_noise, minAL, r5);
    for (auto i : F_r5) std::cout << i << "; ";
    std::cout << std::endl;

    return 0;
}
