
//March 30, 2018
//Timofey Golubev


#include <iostream>
#include "Set_diagonals.h"
#include "Tridiag_solver.h"
#include <vector>

//This sets up and solves the 1D Poisson equation (using supporting functions) for a bilayer device.
//This is made for use with the IQE test, to adjust the electric potential to reflect charge accumulation at the bilayer donor/acceptor interface.

//The layers are assumed to be stacked in the z-direction
//Function requires:
//number of nodes in the device = num_elements (this is # of pts -2 b/c of endpts)
//array of relative dielectric const: epsilon
//left and right bndry voltages: V_leftBC, V_rightBC
//charge DENSITIES (in Coloumbs) at the left and right boundaries of the interface: left_int_charge, right_int_charge
//location of interface: z_int  (this is location of START of 2nd layer, i.e. 1st layer from z =  0  - 24nm, 2nd layer from 25nm-75nm then z_int = 25).

std::vector<double> potential(int num_elements, std::vector<double> &epsilon, double V_leftBC, double V_rightBC, double left_int_charge, double right_int_charge, int z_int)
{
    double epsilon_0 = 8.854e-12; //SI units
    double dx = 1e-9;   //for now, will always be using a 1nm grid spacing

    double CV = -dx*dx/epsilon_0;

    std::vector<double> a (num_elements+1);//main diag
    std::vector<double> b(num_elements); //upper diag, size can be = num_elements b/c off-diags are 1 element less than  main diag
    std::vector<double> c(num_elements);//lower diag
    std::vector<double> potential(num_elements+1); //vector for solution, electric potential
    std::vector<double> rhs(num_elements+1);

    a = set_main_diag(epsilon, a);
    b = set_upper_diag(epsilon, b);
    c = set_lower_diag(epsilon, c);

    //setup rhs of Poisson eqn.
    for (int i = 1;i<= num_elements;i++){
         rhs[i] = 0;  //set all charge densities to 0, then add the interface_charge
    }
    rhs[z_int-1] = CV*left_int_charge;
    rhs[z_int] = CV*right_int_charge;

    std::cout  << "left int charge density" << left_int_charge <<std::endl;

    //test i.e. having a dipole at accross 24-25
    //rhs[24] = 0.681; //corresponds to having 5 holes/(75nm)^2 plane at the interface, with z-mesh size of 1nm
    //rhs[25] = -0.681;

    //for bndrys
    rhs[1] = rhs[1] - epsilon[1]*V_leftBC;
    rhs[num_elements] = rhs[num_elements] - epsilon[num_elements]*V_rightBC;

    potential = TriCRSSolver(a, b, c, rhs);

    //Trying to add these deletes, causes BREAK, doesn't even reach "E potential was recalculated"
    //delete[] a;
    //delete[] b;
    //delete[] c;
    //delete[] rhs;


    /*
    //output the result to terminal
    for(int i = 1;i<= num_elements; i++){
        std::cout << potential[i] << " " <<std::endl;
    }
    */

    return potential;
}
