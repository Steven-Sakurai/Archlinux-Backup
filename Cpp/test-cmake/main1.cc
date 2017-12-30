#include <iostream>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

int main() {
    arma_rng::set_seed_random();  
    
	// Create a 4x4 random matrix
    Mat<double> A = randu(4,4);
    cout << "A:\n" << A << "\n";  
	   	
	// Multiply A with his transpose:
   	cout << "A * A.t() =\n";
   	cout << A * A.t() << "\n";  
	
	// Access/Modify rows and columns from the array:
   	A.row(0) = A.row(1) + A.row(3);
   	A.col(3).zeros();
   	cout << "add rows 1 and 3, store result in row 0, also fill 4th column with zeros:\n";
   	cout << "A:\n" << A << "\n";     
   		
	// Create a new diagonal matrix using the main diagonal of A:
   	Mat<double>B = diagmat(A);
	cout << "B:\n" << B << "\n";   
	
	// Save matrices A and B:
    ofstream file1("A.csv");
    file1 << "X1, X2, X3, X4" << endl;
    ofstream file2("B.csv");
    file2 << "Y1, Y2, Y3, Y4" << endl;
	A.save(file1, csv_ascii);
	B.save(file2, csv_ascii);
	file1.close();
    file2.close();
    return 0;
}
