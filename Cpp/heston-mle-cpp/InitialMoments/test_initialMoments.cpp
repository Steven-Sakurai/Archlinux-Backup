#include "InitialMoments.hpp"

int main() {
	stdVec M(4*21 - 1);
	Parameter my_par;
	initialMoments(my_par, M, true);
	return 0;
}
