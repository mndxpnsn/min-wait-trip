//
//  main.cpp
//  min-wait-trip
//
//  Created by Derek Harrison on 11/06/2022.
//

#include <iostream>
#include <map>
#include <set>
#include <time.h>
#include <vector>

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

    for(int i = 0; i < n; ++i) {
        d[i] = (i + 1) * D / (n + 1);
    }
}

double min_time_rev(double * t, double * d, int n, double D, double dmax, int m, double * dp) {

    double res = 0.0;
    
    int index = n - m - 1;
    
    // Compute number of operations
    ops++;
    
    // Get data from memo table if available
    if(dp[index] != 0.0) {
        return dp[index];
    }
    
    // Boundary
    if(m < 0) {
        return 0;
    }
    
    // Final distance within reach
    if(D - d[index] <= dmax) {
        return t[index];
    }

    // Compute stations in range
    int j = index;
    while(j < n && d[j + 1] - d[index] <= dmax) {
        j++;
    }
    
    int diff = j - index;
    
    // Compute minimum wait time
    for(int del = 1; del <= diff; ++del) {
        double cost = t[index] + min_time_rev(t, d, n, D, dmax, m - del, dp);
        res = min(res, cost);
    }
    
    // Store data in memo table
    dp[index] = res;
    
    return res;
}

double min_time_rec_no_dp(double * t, double * d, int n, double D, double dmax, int m) {
    
    double res = 0.0;
    
    int index = n - m - 1;
    
    // Boundary
    if(m < 0) {
        return 0;
    }
    
    // Last station
    if(m == 0) {
        return t[index];
    }
    
    // Final distance within reach
    if(D - d[index] <= dmax) {
        return t[index];
    }

    // Compute stations in reach
    int j = index;
    while(j < n && d[j + 1] - d[index] <= dmax) {
        j++;
    }
    
    int diff = j - index;
    
    // Compute minimum wait time
    for(int del = 1; del <= diff; ++del) {
        double cost = t[index] + min_time_rec_no_dp(t, d, n, D, dmax, m - del);
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
    
    // Compute stations in range
    int j = 0;
    while(j < n && d[j] <= dmax) {
        j++;
    }
    
    int diff = j;
    
    // Compute minimum wait time
    for(int del = 1; del <= diff; ++del) {
        double cost = min_time_rev(t, d, n, D, dmax, n - del, dp);
        res = min(res, cost);
    }
    
    // Free memo table
    free_dp(dp);
    
    return res;
}

double min_time_no_dp(double * t, double * d, int n, double D, double dmax) {
    double res = 0.0;
    
    int j = 0;
    while(j < n && d[j] <= dmax) {
        j++;
    }
    
    int diff = j;
    
    // Compute minimum wait time
    for(int del = 1; del <= diff; ++del) {
        double cost = min_time_rec_no_dp(t, d, n, D, dmax, n - del);
        res = min(res, cost);
    }
    
    return res;
}

int main(int argc, const char * argv[]) {
    
    // Number of pump stations
    int n = 8;
    
    // Drive distance
    double D = 100.0;
    
    // Max drive distance on full tank
    double dmax = 30.0;
    
    // Declare and initialize wait time and station distances
    double * t = new double[n];
    double * d = new double[n];
    std::vector<int> route;
    
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
