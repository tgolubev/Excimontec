#ifndef SETUP_POISSON_BILAYER_H
#define SETUP_POISSON_BILAYER_H

#include <vector>

std::vector<double> potential(int num_elements, std::vector<double> &epsilon, double V_leftBC, double V_rightBC, double left_int_charge, double right_int_charge, int z_int);


#endif // SETUP_POISSON_BILAYER_H
