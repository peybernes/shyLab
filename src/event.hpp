#ifndef EVENT_HPP
#define EVENT_HPP

#include "timetable.hpp"
#include "serialize_fwd.hpp"

// before or after solution update ?
typedef enum _EventType {PRE, POST} EventType;

class Simulation;

class Event {
public:
  Event();
  Event(Simulation* sim_ptr, const Timetable&, const EventType&);
  virtual ~Event() {};
  virtual void Execute() = 0;
  virtual void Save(ptree&);
  virtual void Finalize();
  void UpdateTimetableClock(RealType dt);
  void ResetTimetableClock();
  Timetable Timetable2();
  EventType Type();
  const std::string name() const;
protected:
  Simulation* m_sim_ptr;
  Timetable m_timetable;
  std::string m_name;
  EventType m_type;
};

#endif // EVENT_HPP
