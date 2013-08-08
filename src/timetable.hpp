/// Timetable controls when an event starts, stops and happens.
#ifndef TIMETABLE_HPP_
#define TIMETABLE_HPP_

#include "clock.hpp"
#include "serialize_fwd.hpp"

class Timetable {
public:
  Timetable();
  Timetable(int, int, int, RealType, RealType, RealType);
  void Update(RealType);
  void ResetLocalClock();
  void UpdateLocalClock(RealType);
  bool Happens(const Clock&);
  bool IsFinished(const Clock&);
  void Save(ptree&);
  void Load(ptree&);
  int istart();
  int istart() const;
private:
  RealType m_start, m_end, m_step;
  int m_istart, m_iend, m_istep;
  /// Local clock. Stores time and iteration count since last time
  /// event was executed.
  Clock m_local_clock;
};

#endif // TIMETABLE_HPP_
