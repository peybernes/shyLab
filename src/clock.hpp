#ifndef CLOCK_HPP
#define CLOCK_HPP

#include "serialize_fwd.hpp"

class Clock {
public:
  Clock();
  Clock(RealType t, int i);
  void Update(RealType dt);
  void Reset();
  void Save(ptree &pt);
  void Load(ptree &pt);
  RealType time();
  RealType time() const;
  int iter();
  int iter() const;
private:
  RealType m_time;
  int m_iter;
};

#endif // CLOCK_HPP
