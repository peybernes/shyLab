#include "clock.hpp"
#include "serialize.hpp"

Clock::Clock():
  m_time(0.0), m_iter(0) {}

Clock::Clock(RealType t, int i):
  m_time(t), m_iter(i) {}

void Clock::Update(RealType dt) {

  m_time += dt;
  m_iter += 1;

}
 
void Clock::Reset() {

    m_time = 0.0;
    m_iter = 0;

}

void Clock::Save(ptree &pt) {

  pt.put("time", m_time);
  pt.put("iteration", m_iter);

}

void Clock::Load(ptree &pt) {

  m_time = pt.get<RealType>("time", 0.0);
  m_iter = pt.get<int>("iteration", 0);

}

RealType Clock::time() {
  return m_time;
}

RealType Clock::time() const {
  return m_time;
}

int Clock::iter() {
  return m_iter;
}

int Clock::iter() const {
  return m_iter;
}
