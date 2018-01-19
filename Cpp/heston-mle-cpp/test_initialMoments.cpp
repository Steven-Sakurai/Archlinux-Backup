#include "InitialMoments.hpp"

int main() {
	stdMat M(4, stdVec(21, 0));
	Parameter my_par;
	initialMoments(my_par, M, true);
	return 0;
}
