#include <boost/python.hpp>
#include "HeurClsp.hpp"

using namespace boost::python;

BOOST_PYTHON_MODULE(heurclsp)
{
    class_<HeurClsp>("heurclsp",init<list, list, list, list, list, list, list, int, int, int, int, float, float>())
        .def("thomas", &HeurClsp::thomas)
        .def("heursolver", &HeurClsp::heursolver)
        .def("objective", &HeurClsp::objective)
    ;
}
