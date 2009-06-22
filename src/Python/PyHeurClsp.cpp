#include <boost/python.hpp>
#include "HeurClsp.hpp"

using namespace boost::python;

BOOST_PYTHON_MODULE(heurclsp)
{
    class_<HeurClsp>("heurclsp",init<list, list, list, list, list, list, list, int, int, int, int, float, float>())
        .def_readonly("price",&HeurClsp::getPrice)
        .def_readonly("hold",&HeurClsp::getHold)
        .def_readonly("prod",&HeurClsp::getProd)
        .def_readonly("setup",&HeurClsp::getSetup)
        .def_readonly("kkt",&HeurClsp::getCoef)
        .def("thomas", &HeurClsp::thomas)
        .def("heursolver", &HeurClsp::heursolver)
        .def("objective", &HeurClsp::objective)
    ;
}
