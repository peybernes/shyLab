#include "timetable.hpp"

#include <iostream>

#include "config.hpp"
#include "serialize.hpp"

Timetable::Timetable():
  m_start(0.0), m_end(INFTY), m_step(INFTY),
  m_istart(0), m_iend(MAXINT), m_istep(MAXINT) {}

Timetable::Timetable(int istart, int istep,
	     int iend, RealType start,
	     RealType step, RealType end):
  m_istart(istart), m_istep(istep), m_iend(iend),
  m_start(start), m_step(step), m_end(end) {};

void Timetable::UpdateLocalClock(RealType dt) {

  m_local_clock.Update(dt);

}

void Timetable::ResetLocalClock() {
  m_local_clock.Reset();
}

bool Timetable::Happens(const Clock& clock) {

  bool result =
    (clock.time() >= m_start) &&
    (clock.time() <= m_end) &&
    (clock.iter() >= m_istart) &&
    (clock.iter() <= m_iend) &&
    ((clock.iter() == 0) ||
     (m_istep == 1) ||
     ((m_local_clock.iter() == m_istep) ||
      (m_local_clock.time() >= m_step)));

  return result;
}

bool Timetable::IsFinished(const Clock& clock) {

  return (clock.time() > m_end) || (clock.iter() > m_iend);

}

void Timetable::Save(ptree &pt) {
  
  pt.put("start", m_start);
  pt.put("end", m_end);
  pt.put("step", m_step);
  pt.put("istart", m_istart);
  pt.put("iend", m_iend);
  pt.put("istep", m_istep);

}

void Timetable::Load(ptree &pt) {
  
  m_start = pt.get<RealType>("start", 0.0);
  m_end = pt.get<RealType>("end", INFTY);
  m_step = pt.get<RealType>("step", INFTY);
  m_istart = pt.get<int>("istart", 0);
  m_iend = pt.get<int>("iend", MAXINT);
  m_istep = pt.get<int>("istep", MAXINT);

}

int Timetable::istart() {
  return m_istart;
}

int Timetable::istart() const {
  return m_istart;
}

int Timetable::istep() {
  return m_istep;
}

int Timetable::istep() const {
  return m_istep;
}
