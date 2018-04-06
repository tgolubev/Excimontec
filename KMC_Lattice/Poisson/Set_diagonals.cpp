#include <iostream>
#include<vector>


//this is a in tridiag_solver
std::vector<double> set_main_diag(std::vector<double> &epsilon, std::vector<double> &main_diag){
    int num_elements = main_diag.size() - 1;

    for (int i=1; i<= num_elements;i++){
        main_diag[i] = -2.*epsilon[i];
    }
    return main_diag;
}

//this is b in tridiag_solver
std::vector<double> set_upper_diag(std::vector<double> &epsilon, std::vector<double> &upper_diag){
    int num_elements = upper_diag.size() -1;  //-1 since vectors fill from 0, but I'm indexing from 1
    //NOTE: num_elements here, is actually # of elements in the off-diagonal
    
    for (int i = 1; i<=num_elements; i++){
        upper_diag[i] = epsilon[i+1];    //1st element here corresponds to 2nd row of matrix
    }
    return upper_diag;
}

//this is c in tridiag_solver
std::vector<double> set_lower_diag(std::vector<double> &epsilon, std::vector<double> &lower_diag){
    int num_elements = lower_diag.size() -1;
    
    for (int i = 1; i<=num_elements; i++){
        lower_diag[i] = epsilon[i-1];    //1st element here corresponds to the 1st row of matrix: so need to use i.e. epsilon corresponding to fullV(3)
    }
    return lower_diag;
    
}
