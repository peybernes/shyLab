#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <vector>
#include <string>

#include "event.hpp"
#include "serialize_fwd.hpp"

enum IoFormatType {ASCII, BINARY};

class Output : public Event {
public:
  Output(Simulation*, std::string&, const Timetable&);
  virtual ~Output() {};
  virtual void Execute() = 0;
  virtual void Save(ptree&);
protected:
  std::string m_stream_name;
  IoFormatType m_io_format;
};

class OutputSimulation : public Output {
public:
  OutputSimulation(Simulation*, std::string&, const std::string&, const Timetable&);
  void Execute();
};

class OutputMinMax : public Output {
public:
  OutputMinMax(Simulation*, std::string&, const Timetable&);
  void Execute();
};

class OutputLocation : public Output {
public:
  OutputLocation(Simulation*, std::string&, const Timetable&,
		 const RealType, const RealType);
  void Execute();
  RealType m_x, m_y;
private:
  int m_cell_id;
};

class OutputTime : public Output {
public:
  OutputTime(Simulation*, std::string&, const Timetable&);
  void Execute();
};

class OutputNorms : public Output {
public:
  OutputNorms(Simulation*, 
	      std::string&, 
	      const Timetable&,
	      const std::string& name1, 
	      const std::string& name2);
  void Execute();
  void Save(ptree& pt);
private:
  std::string m_name1;
  std::string m_name2;
  RealType* var1;
  RealType* var2;
};

#endif // OUTPUT_HPP

