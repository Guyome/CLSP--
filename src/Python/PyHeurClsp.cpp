#include <boost/python.hpp>
#include "HeurClsp.hpp"

using namespace boost::python;

void (HeurClsp::*thomasx1)() = &HeurClsp::thomas;
void (HeurClsp::*thomasx2)(double) = &HeurClsp::thomas;
void (HeurClsp::*thomasx3)(list) = &HeurClsp::thomas;

BOOST_PYTHON_MODULE(heurclsp)
{
    scope().attr("__doc__") =
        "This library implements some algorithms\n"
        " for lot-size, capacited lot-size, profit maximizing lot-size\n"
        " and profit maximizing  capacited lot-size problem for continus\n"
        " and discret prices.\n\n"
        " It uses Ipopt solver from COINOR project.\n";
    class_<HeurClsp>("pclsp",init<list, list, list, list, list, list, list, int, int, int, int, float, float>())
        .add_property("verbose",&HeurClsp::getVerbosity,&HeurClsp::setVerbosity)
        .add_property("nb_cycle",&HeurClsp::getNbCycle,&HeurClsp::setNbCycle)
        .add_property("min_diff",&HeurClsp::getStopDiff,&HeurClsp::setStopDiff)
        .add_property("smooth_param",&HeurClsp::getSmooth,&HeurClsp::setSmooth)
        .def_readonly("price",&HeurClsp::getPrice)
        .def_readonly("hold",&HeurClsp::getHold)
        .def_readonly("prod",&HeurClsp::getProd)
        .def_readonly("setup",&HeurClsp::getSetup)
        .def_readonly("kkt",&HeurClsp::getCoef)
        .def_readonly("gap",&HeurClsp::getGAP)
        .def("Objective", &HeurClsp::objective)
        .def("noQP", &HeurClsp::setHeur, "Use an heuristic in spite of QP solver in PCLSPSolvHeur()\n")
        .def("isLSP",&HeurClsp::isWW, "Returns True if you are in lot-size case\n")
        .def("isDiscret",&HeurClsp::isDiscret, "Returns True if you are using discret price")
        .def("DynProgSolv", thomasx1, "Solser based on dynamic programming\n"
            "Solve profit lot-size problem\nSee. J. Thomas(1970)\n")
        .def("DynProgSolv", thomasx2, "Solser based on dynamic programming\n"
            "Solve lot-size problem\nSee Wagner and Within(1958)\n")
        .def("DynProgSolv", thomasx3, "Solser based on dynamic programming\n"
            "Modificatied Thomas's algorithm for discrets prices'\n")
        .def("PCLSPSolvHeur", &HeurClsp::heursolver, "Heuristic for PCLSP developped by K.K. Haugen and al.")
    ;
}
