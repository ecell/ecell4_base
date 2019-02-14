#ifndef ECELL4_PYTHON_API_OBSERVER_HPP
#define ECELL4_PYTHON_API_OBSERVER_HPP

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <ecell4/core/observers.hpp>

namespace py = pybind11;

namespace ecell4
{

namespace python_api
{

    template<class Base = ecell4::Observer>
    class PyObserver: public Base
    {
    public:
        using Base::Base;

        const Real next_time() const override
        {
            PYBIND11_OVERLOAD(const Real, Base, next_time,);
        }

        void reset() override
        {
            PYBIND11_OVERLOAD(void, Base, reset,);
        }
    };

    class FixedIntervalPythonHooker
        : public FixedIntervalObserver
    {
    public:
        using base_type = FixedIntervalObserver;
        using callback_t = std::function<bool (const boost::shared_ptr<WorldInterface>&, bool)>;

        FixedIntervalPythonHooker(const Real& dt, callback_t callback)
            : base_type(dt), callback_(callback)
        {
        }

        bool fire(const Simulator* sim, const boost::shared_ptr<WorldInterface>& world) override
        {
            return callback_(world, sim->check_reaction()) && base_type::fire(sim, world);
        }

    protected:
        callback_t callback_;

    };

    void define_observers(py::module& m)
    {
        py::class_<Observer, PyObserver<>, boost::shared_ptr<Observer>>(m, "Observer")
            .def("next_time", &Observer::next_time)
            .def("reset", &Observer::reset)
            .def("num_steps", &Observer::num_steps);

        py::class_<FixedIntervalPythonHooker, PyObserver<FixedIntervalPythonHooker>,
            boost::shared_ptr<FixedIntervalPythonHooker>>(m, "FixedIntervalPythonHooker")
            .def(py::init<const Real&, FixedIntervalPythonHooker::callback_t>());

        py::class_<FixedIntervalNumberObserver, PyObserver<FixedIntervalNumberObserver>,
            boost::shared_ptr<FixedIntervalNumberObserver>>(m, "FixedIntervalNumberObserver")
            .def(py::init<const Real&>())
            .def(py::init<const Real&, const std::vector<std::string>&>())
            .def("data", &FixedIntervalNumberObserver::data)
            .def("targets", &FixedIntervalNumberObserver::targets)
            .def("save", &FixedIntervalNumberObserver::save);

        py::class_<NumberObserver, PyObserver<NumberObserver>,
            boost::shared_ptr<NumberObserver>>(m, "NumberObserver")
            .def(py::init<>())
            .def(py::init<const std::vector<std::string>&>())
            .def("data", &NumberObserver::data)
            .def("targets", &NumberObserver::targets)
            .def("save", &NumberObserver::save);

        py::class_<TimingNumberObserver, PyObserver<TimingNumberObserver>,
            boost::shared_ptr<TimingNumberObserver>>(m, "TimingNumberObserver")
            .def(py::init<const std::vector<Real>&>())
            .def(py::init<const std::vector<Real>&, const std::vector<std::string>&>())
            .def("data", &TimingNumberObserver::data)
            .def("targets", &TimingNumberObserver::targets)
            .def("save", &TimingNumberObserver::save);

        py::class_<FixedIntervalHDF5Observer, PyObserver<FixedIntervalHDF5Observer>,
            boost::shared_ptr<FixedIntervalHDF5Observer>>(m, "FixedIntervalHDF5Observer")
            .def(py::init<const Real&, const std::string&>())
            .def("prefix", &FixedIntervalHDF5Observer::prefix)
            .def("filename", (const std::string (FixedIntervalHDF5Observer::*)() const) &FixedIntervalHDF5Observer::filename)
            .def("filename", (const std::string (FixedIntervalHDF5Observer::*)(const Integer) const) &FixedIntervalHDF5Observer::filename);

        py::class_<FixedIntervalCSVObserver, PyObserver<FixedIntervalCSVObserver>,
            boost::shared_ptr<FixedIntervalCSVObserver>>(m, "FixedIntervalCSVObserver")
            .def(py::init<const Real&, const std::string&>())
            .def(py::init<const Real&, const std::string&, std::vector<std::string>&>())
            .def("log", &FixedIntervalCSVObserver::log)
            .def("filename", &FixedIntervalCSVObserver::filename)
            .def("set_header", &FixedIntervalCSVObserver::set_header)
            .def("set_formatter", &FixedIntervalCSVObserver::set_formatter);

        py::class_<CSVObserver, PyObserver<CSVObserver>, boost::shared_ptr<CSVObserver>>(m, "CSVObserver")
            .def(py::init<const std::string&>())
            .def(py::init<const std::string&, std::vector<std::string>&>())
            .def("log", &CSVObserver::log)
            .def("filename", &CSVObserver::filename)
            .def("set_header", &CSVObserver::set_header)
            .def("set_formatter", &CSVObserver::set_formatter);

        py::class_<FixedIntervalTrajectoryObserver, PyObserver<FixedIntervalTrajectoryObserver>,
            boost::shared_ptr<FixedIntervalTrajectoryObserver>>(m, "FixedIntervalTrajectoryObserver")
            .def(py::init<const Real&, const std::vector<ParticleID>&>())
            .def(py::init<const Real&, const std::vector<ParticleID>&, const bool>())
            .def(py::init<const Real&, const std::vector<ParticleID>&, const bool, const Real>())
            .def(py::init<const Real&, const bool>())
            .def(py::init<const Real&, const bool, const Real>())
            .def("data", &FixedIntervalTrajectoryObserver::data)
            .def("num_tracers", &FixedIntervalTrajectoryObserver::num_tracers)
            .def("t", &FixedIntervalTrajectoryObserver::t);

        py::class_<TimingTrajectoryObserver, PyObserver<TimingTrajectoryObserver>,
            boost::shared_ptr<TimingTrajectoryObserver>>(m, "TimingTrajectoryObserver")
            .def(py::init<const std::vector<Real>&, const std::vector<ParticleID>&>())
            .def(py::init<const std::vector<Real>&, const std::vector<ParticleID>&, const bool>())
            .def(py::init<const std::vector<Real>&, const std::vector<ParticleID>&, const bool, const Real>())
            .def(py::init<const std::vector<Real>&, const bool>())
            .def(py::init<const std::vector<Real>&, const bool, const Real>())
            .def("data", &TimingTrajectoryObserver::data)
            .def("num_tracers", &TimingTrajectoryObserver::num_tracers)
            .def("t", &TimingTrajectoryObserver::t);

        py::class_<FixedIntervalTrackingObserver, PyObserver<FixedIntervalTrackingObserver>,
            boost::shared_ptr<FixedIntervalTrackingObserver>>(m, "FixedIntervalTrackingObserver")
            .def(py::init<const Real&, const std::vector<Species>&>())
            .def(py::init<const Real&, const std::vector<Species>&, const bool&>())
            .def(py::init<const Real&, const std::vector<Species>&, const bool&, const Real>())
            .def(py::init<const Real&, const std::vector<Species>&, const bool&, const Real, const Real>())
            .def("data", &FixedIntervalTrackingObserver::data)
            .def("num_tracers", &FixedIntervalTrackingObserver::num_tracers)
            .def("t", &FixedIntervalTrackingObserver::t);

        py::class_<TimeoutObserver, PyObserver<TimeoutObserver>, boost::shared_ptr<TimeoutObserver>>(m, "TimeoutObserver")
            .def(py::init<>())
            .def(py::init<const Real>())
            .def("interval", &TimeoutObserver::interval)
            .def("duration", &TimeoutObserver::duration)
            .def("accumulation", &TimeoutObserver::accumulation);
    }

}

}

#endif /* ECELL4_PYTHON_API_OBSERVER_HPP */
