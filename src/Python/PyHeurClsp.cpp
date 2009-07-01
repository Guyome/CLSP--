#include <boost/python.hpp>
#include "HeurClsp.hpp"

using namespace boost::python;

void (HeurClsp::*thomasx1)() = &HeurClsp::thomas;
void (HeurClsp::*thomasx2)(double) = &HeurClsp::thomas;
void (HeurClsp::*thomasx3)(list) = &HeurClsp::thomas;

BOOST_PYTHON_MODULE(heurclsp)
{
    class_<HeurClsp>("heurclsp",init<list, list, list, list, list, list, list, int, int, int, int, float, float>())
        .add_property("verbose",&HeurClsp::getVerbosity,&HeurClsp::setVerbosity)
        .add_property("nbcycle",&HeurClsp::getNbCycle,&HeurClsp::setNbCycle)
        .add_property("mindiff",&HeurClsp::getStopDiff,&HeurClsp::setStopDiff)
        .add_property("smoothparam",&HeurClsp::getSmooth,&HeurClsp::setSmooth)
        .def_readonly("price",&HeurClsp::getPrice)
        .def_readonly("hold",&HeurClsp::getHold)
        .def_readonly("prod",&HeurClsp::getProd)
        .def_readonly("setup",&HeurClsp::getSetup)
        .def_readonly("kkt",&HeurClsp::getCoef)
        .def("thomas", thomasx1)
        .def("thomas", thomasx2)
        .def("thomas", thomasx3)
        .def("heursolver", &HeurClsp::heursolver)
        .def("objective", &HeurClsp::objective)
        .def("useheur", &HeurClsp::setHeur)
        .def("isWW",&HeurClsp::isWW)
        .def("isDiscret",&HeurClsp::isDiscret)
    ;
}
