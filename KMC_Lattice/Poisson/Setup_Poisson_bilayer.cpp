
//March 30, 2018
//Timofey Golubev


#include <iostream>
#include "Set_diagonals.h"
#include "Tridiag_solver.h"

//This sets up and solves the 1D Poisson equation (using supporting functions) for a bilayer device.
//This is made for use with the IQE test, to adjust the electric potential to reflect charge accumulation at the bilayer donor/acceptor interface.

//The layers are assumed to be stacked in the z-direction
//Function requires:
//number of nodes in the device = num_elements (this is # of pts -2 b/c of endpts)
//array of relative dielectric const: epsilon
//left and right bndry voltages: V_leftBC, V_rightBC
//charge DENSITIES (in Coloumbs) at the left and right boundaries of the interface: left_int_charge, right_int_charge
//location of interface: z_int  (this is location of START of 2nd layer, i.e. 1st layer from z =  0  - 24nm, 2nd layer from 25nm-75nm then z_int = 25).

double * potential(int num_elements, double* epsilon, double V_leftBC, double V_rightBC, double left_int_charge, double right_int_charge, int z_int)
{

    //double e = 1.6e-19;  //elementary charge (C)
    //int num_elements = 100;  //number of nodes in the system
    //double *epsilon = new double[num_elements+2];
    double *a = new double[num_elements]; //main diag
    double *b = new double[num_elements-1]; //upper diag
    double *c = new double[num_elements-1]; //lower diag
    double *potential = new double[num_elements];  //array for solution, electric potential
    double *rhs = new double[num_elements];

    //double V_leftBC = 0.;
    //double V_rightBC= 0.9;
    /*
    for(int i=0;i<=num_elements+1;i++){
        epsilon[i] = 3.8;
    }
    */

    a = set_main_diag(epsilon, a,  num_elements);
    b = set_upper_diag(epsilon, b, num_elements);
    c = set_lower_diag(epsilon, c, num_elements);

    //setup rhs of Poisson eqn.
    for (int i = 1;i<= num_elements;i++){
         rhs[i] = 0;  //set all charge densities to 0, then add the interface_charge
    }
    rhs[z_int-1] = left_int_charge;
    rhs[z_int] = right_int_charge;

    //test i.e. having a dipole at accross 24-25
    //rhs[24] = 0.681; //corresponds to having 5 holes/(75nm)^2 plane at the interface, with z-mesh size of 1nm
    //rhs[25] = -0.681;

    //for bndrys
    rhs[1] = rhs[1] - epsilon[1]*V_leftBC;
    rhs[num_elements] = rhs[num_elements] - epsilon[num_elements]*V_rightBC;

    potential = TriCRSSolver(a, b, c, rhs,  num_elements);

    /*
    //output the result to terminal
    for(int i = 1;i<= num_elements; i++){
        std::cout << potential[i] << " " <<std::endl;
    }
    */

    return potential;
}
