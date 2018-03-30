#include <boost/python.hpp>

char const* Cgreet() {
    return "hello, world!";
}

BOOST_PYTHON_MODULE(hello_pylib)
{
    using namespace boost::python;
    def("greet", Cgreet);
}


