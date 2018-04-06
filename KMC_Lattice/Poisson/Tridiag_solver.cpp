// % !==================================================================================
// % ! CRSSolver = Compressed Row Storage (CRS] + Linear Equations Solver
// % !
// % ! This is a solver for a set of linear equations in the form of matrix
// % ! Ax = b  where A = nxn matrix, x and b is vector.
// % ! where matrix A is a "TRIDIAGONAL matrix".
// % ! 
// % !
//    This uses the technique from "Inversion of Jacobi's Tridiagonal Matrix" R. Usmani, 1993
//   The variable notation is almost the same, except, here zeta is theta in Usmani.
//
//% ! NOTE: Inverse matrix will be calculated row by row, then used to calculate the
//% ! answer x element by element. I'll not keep the whole inverse matrix. 


//TriCRSSolver will return a pointer to the answer array x. Takes pointers to the arrays as input
// a = array containing elements of main diagonal. indices: (a1.....an)
// b = array containing elements of upper diagonal. indices (b1....b_n-1)
//c = array containing elements of lower diagonal. indices (c1...c_n-1)

// % !==================================================================================


#include <iostream>
#include <numeric>  //allows to use inner product fnc.
#include <vector> //using STL vector prevents issues with pointers, information about size of the array is not lost as with C-style arrays. Also has a size() fnc
//and allows to allocate with a variable length
#include <cmath>


std::vector<double> TriCRSSolver(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &rhs){

    const int num_elements = a.size() -1; //# of elements in a is +1 of # of elements in diagonal b/c a[0] is not used but exists

    // valid indices in a vector fall into the range 0 to vectorSize - 1. So if want to index 1 to num_elements, vectors must be num_elements+1 sized
    std::vector<double> zeta(num_elements+2);
    std::vector<double> phi(num_elements+2);
    std::vector<double> A_inverse(num_elements+1);
    std::vector<double> x(num_elements+1);
    double temp;

    zeta[0] = 1.0;
    zeta[1] = a[1];
    for(int i = 2; i<=num_elements; i++){
        zeta[i] = a[i]*zeta[i-1] - b[i-1]*c[i-1]*zeta[i-2];
    }
   
    phi[num_elements+1] = 1.0;
    phi[num_elements] = a[num_elements];
    
    for(int i = num_elements-1; i>0; i--){  //decreasing iterations.., for recursion relation
        phi[i] = a[i]*phi[i+1] - b[i]*c[i]*phi[i+2];
    }

   //calculate inverse matrix 
    for(int i = 1; i<=num_elements; i++){
        for (int j = 1; j<= i-1; j++){
            temp = 1.0;
            for(int k = j; k<=i-1; k++){  //k<= i-1 b/c there's 1 less element in off-diagonals than main diagonal (which is indexed by i)
                temp = temp*c[k];   //multiplies all the c's together--> also is what sometimes causes Usmani method to blow up
            }
            A_inverse[j] = pow(-1., i+j)*temp*zeta[j-1]*phi[i+1]/zeta[num_elements];
            //std::cout <<"a inverse " << A_inverse[j] <<std::endl;
            //here seems fine
        }
        for(int j = i; j<=num_elements;j++){
            temp = 1.0;
            for(int k = i; k<=j-1;k++){
                temp = temp * b[k];
            }
            A_inverse[j] = pow(-1., i+j)*temp*zeta[i-1]*phi[j+1]/zeta[num_elements];
            //std::cout <<"a inverse " << A_inverse[j] <<std::endl;
        }
        //hardcode the dot product
        x[i] = 0;
        for(int m = 1; m<= num_elements;m++){
            x[i] += A_inverse[m]*rhs[m];
        }
    }
    //%test[i] = A_inverse*rhs
return x;
}

