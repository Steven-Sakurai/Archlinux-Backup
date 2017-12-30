#include <cstdlib>
#include <string>

using namespace std;

int main() {
    string r_script = "/usr/local/bin/Rscript --vanilla ./ISLR.r";
    system(r_script.c_str());
    return 0;
}
