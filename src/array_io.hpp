#ifndef ARRAY_IO_HPP
#define ARRAY_IO_HPP

#include <string>
#include <vector>

void ReadTxtAsciiScalar(int n, std::istream* is_ptr, RealType* tab);

#ifdef CATALINA2
void ReadBCFile(const std::ifstream& ifs, 
		std::vector<RealType>* times_ptr, 
		std::vector <RealType>* values_ptr);
#endif // CATALINA2

void WriteVTKAsciiTab(int n, const RealType* tab,
		      const std::string& variable_name, 
		      std::ostream* os_ptr);

void WriteVTKAsciiVectorTab(int n, const RealType* tab_x,
			    const RealType* tab_y, 
			    const std::string& variable_name,
			    std::ostream* os_ptr);

void WriteVtkXmlScalarHeader(const std::string& variable_name,
			     std::ostream* os_ptr);

void WriteVtkXmlAsciiScalar(int n, const RealType* v, 
			    const std::string& variable_name,
			    std::ostream* os_ptr);

void WriteVtkXmlVectorHeader(const std::string& variable_name,
			     std::ostream* os_ptr);

void WriteVtkXmlAsciiVector(int n, const RealType* vx, 
			    const RealType* vy,
			    const RealType* vz, 
			    const std::string& variable_name,
			    std::ostream* os_ptr);

void WriteVTKBinaryTab(int n, const RealType *v,
		       const std::string& variable_name, 
		       bool swap, std::ostream* os_ptr);

#endif // ARRAY_IO_HPP_
