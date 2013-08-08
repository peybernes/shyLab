/// Simulation events, e.g. mainly outputs, and initialisation of
/// variables for now.
#include "event.hpp"

#include "serialize.hpp"

Event::Event(Simulation* sim_ptr, const Timetable& timetable, 
	     const EventType& type):
  m_type(type), m_timetable(timetable), m_sim_ptr(sim_ptr) {}

void Event::Save(ptree &pt) {

  ptree pt_timetable;
  m_timetable.Save(pt_timetable);

  pt.put_child("Timetable", pt_timetable);

}

EventType Event::Type() {
  return m_type;
}

const std::string Event::name() const {
  return m_name;
}

void Event::UpdateTimetableClock(RealType dt) {
  m_timetable.UpdateLocalClock(dt);
}

void Event::ResetTimetableClock() {
  m_timetable.ResetLocalClock();
}


Timetable Event::Timetable2() {
  return m_timetable;
}

void Event::Finalize() {
}


