#include "distlink-2.0/distlink.h"
#include "iostream"

int main(){
    // Define two orbits
    COrbitData o1 = COrbitData(7149.23810, 0.01951, 1.7164, 2.853, 0.7038);
    COrbitData o2 = COrbitData(7254.82582, 0.02770, 1.7123, 2.85, 0.188);

    // Change min mutal inclination value
    SMOIDResult res = MOID_fast(o1, o2, 1e-10, 1e-10);

    std::cout << "Result" << std::endl;
    std::cout << "Distance: " << res.distance << "(" << res.distance_error << ")" << std::endl;
    std::cout << "u1: " << res.u1 << "(" << res.u1_error << ")" << std::endl;
    std::cout << "u2: " << res.u2 << "(" << res.u2_error << ")" << std::endl;
}