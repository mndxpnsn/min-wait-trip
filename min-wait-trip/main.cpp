//
//  main.cpp
//  min-wait-trip
//
//  Created by Derek Harrison on 11/06/2022.
//

#include <iostream>
#include <time.h>

int ops = 0;

double min(double x, double y) {
    double res = 0.0;
    
    // Cut out zeros if one argument is nonzero
    if(x > 0 && y == 0) {
        return x;
    }
    
    if(x == 0 && y > 0) {
        return y;
    }
    
    if(x < y) res = x;
    else res = y;
    
    return res;
}

void init_t(int n, double * t) {
    double min_t = 0;
    double max_t = 20.0;
    
    srand((unsigned) time(NULL));

    for(int i = 0; i < n; ++i) {
        double f1 = (double) rand() / RAND_MAX;
        double t_loc = min_t + f1 * (max_t - min_t);
        
        t[i] = t_loc;
    }
}

void init_d(int n, double * d, double D, double dmax) {
    double min_x = 0;
    double max_x = dmax;
    double min_b = 0.0;
    
    srand((unsigned) time(NULL));

    for(int i = 0; i < n; ++i) {
        double f1 = (double) rand() / RAND_MAX;
        double x = min_b + f1 * (max_x - min_x);

        d[i] = x;
        
        std::cout << "x: " << x << std::endl;
        
        min_b = x;
    }
    
    
}

double min_time_rec(double * t, double * d, int n, double D, double dmax, int s, double dist, double * dp) {
    double res = 0.0;
    
    ops++;
    
    // Store data in memo table
    if(dp[s + 1] != 0) {
        return dp[s + 1];
    }
    
    // Determine stations in range
    int ns = s + 1;
    while(ns < n && d[ns] - dist <= dmax) {
        ns++;
    }
    
    // Compute optimum route
    for(int st = s + 1; st < n && st <= ns; ++st) {
        double cost = t[st] + min_time_rec(t, d, n, D, dmax, st, d[st], dp);
        res = min(res, cost);
    }
    
    // Store data in memo table
    dp[s + 1] = res;
    
    return res;
}

double min_time_rec_no_dp(double * t, double * d, int n, double D, double dmax, int s, double dist) {
    double res = 0.0;
    
    // Determine stations in range
    int ns = s + 1;
    while(ns < n && d[ns] - dist <= dmax) {
        ns++;
    }
    
    // Compute optimum route
    for(int st = s + 1; st < n && st <= ns; ++st) {
        double cost = t[st] + min_time_rec_no_dp(t, d, n, D, dmax, st, d[st]);
        res = min(res, cost);
    }
    
    return res;
}

double * init_dp(int n) {
    
    double * dp = new double[n];
    
    for(int i = 0; i < n; ++i) {
        dp[i] = 0;
    }
    
    return dp;
}

void free_dp(double * dp) {
    
    delete [] dp;
}

double min_time(double * t, double * d, int n, double D, double dmax) {
    double res = 0.0;
    
    // Initialize memo table
    double * dp = init_dp(n + 1);
    
    // Compute minimum wait time
    res = min_time_rec(t, d, n, D, dmax, -1, 0.0, dp);
    
    // Free memo table
    free_dp(dp);
    
    return res;
}

double min_time_no_dp(double * t, double * d, int n, double D, double dmax) {
    double res = 0.0;
    
    // Compute minimum wait time
    res = min_time_rec_no_dp(t, d, n, D, dmax, -1, 0.0);
    
    return res;
}

int main(int argc, const char * argv[]) {
    
    // Number of pump stations
    int n = 29;
    
    // Drive distance
    double D = 100.0;
    
    // Max drive distance on full tank
    double dmax = 10.0;
    
    // Declare and initialize wait time and station distances
    double * t = new double[n];
    double * d = new double[n];
    
    init_t(n, t);
    init_d(n, d, D, dmax);
    
    // Compute minimum waiting time
    double min_wait = min_time(t, d, n, D, dmax);
    double min_wait_no_dp = min_time_no_dp(t, d, n, D, dmax);
    
    // Print results
    std::cout << "min wait time: " << min_wait << std::endl;
    std::cout << "min wait time no dp: " << min_wait_no_dp << std::endl;
    std::cout << "number of operations: " << ops << std::endl;
    std::cout << "n * n: " << n * n << std::endl;
    
    return 0;
}
