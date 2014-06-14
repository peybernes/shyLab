/// \file arrayIo.hpp Input/output routines for arrays.

#include "array_io.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "byte_swap.hpp"
#include "token_parser.hpp"

#define STRINGIFY(s) #s

void ReadTxtAsciiScalar(int n, std::istream *is_ptr, RealType* tab) {
  
  if (is_ptr->bad())
    std::cerr << "Could not read from file. Aborting...\n";

  RealType value = 0.0;

  // String holding the current line.
  std::string line = "";
  
  int cnt = 0;

  while (getline(*is_ptr, line)) {

    if (line.substr(0, 2) == std::string("//")) // ignore comments
      {}

    else {

      std::istringstream isstream(line);
      std::string token;
      getline(isstream, token);
	
      if (!ParseToken<RealType>(token, &value))
	std::cerr << "Bad input: \"" << token << "\", expected a floating point number\n";

      tab[cnt] = value;
    }
    
    ++cnt;
  }

  if (n != cnt) {

    std::cerr << "ERROR: n=" << n << ", cnt=" << cnt << "\n";

    assert(n == cnt);

  }

}

void WriteTxtAsciiScalar(int n, const RealType* tab, std::ostream* os_ptr) {

  if (os_ptr->bad())
    std::cerr << "Could not write to file. Aborting...\n";
  
  for (int i = 0; i < n; ++i)
    *os_ptr << tab[i] << "\n";

}

void WriteVTKAsciiTab(int n, const RealType* v,
		      const std::string& variable_name,
		      std::ostream* os_ptr) {

  *os_ptr << "SCALARS " << variable_name << " "
	  << STRINGIFY(RealType)
    // no space here, putting a space bugs with new Paraview versions
	  << n << "\n"
	  << "LOOKUP_TABLE DEFAULT\n";
 
  for (int i = 0; i < n; ++i)
    *os_ptr << v[i] << "\n";

  *os_ptr << "\n";
}

/// Write vector value in VTK ascii format.
void WriteVTKAsciiVectorTab(int n, 
			    const RealType* tab_x,
			    const RealType* tab_y,
			    const std::string& variable_name,
			    std::ostream* os_ptr) {

  *os_ptr << "VECTORS " << variable_name << " float" << n << "\n";

  for (int i = 0; i < n; ++i)
    *os_ptr << tab_x[i] << " " << tab_y[i] << " 0.0\n";

  *os_ptr << "\n\n";

}

void WriteVTKBinaryTab(std::ostream &os, int n,
		       const RealType *v,
		       const std::string& var_name, bool swap) {

  os << "SCALARS " << var_name << " "
     << STRINGIFY(s)
     << n << "\n"
     << "LOOKUP_TABLE DEFAULT\n";
 
  for (int i = 0; i < n; ++i) {

    RealType value = v[i];

    // swap bytes before writing.
    if (swap == true)
      value = ByteSwap<RealType>(value);

    os.write((char*)&value, sizeof(RealType));

  }    

  os << "\n";
}

void WriteVtkXmlScalarHeader(const std::string& variable_name,
			     std::ostream* os_ptr) {

  *os_ptr << "type=\"Float64\" "
	  << "Name=\"" << variable_name << "\" "
	  << "NumberOfComponents=\"1\" "
	  << "format=\"ascii\"";

}

void WriteVtkXmlAsciiScalar(int n, const RealType* v, const
			    std::string& variable_name,
			    std::ostream* os_ptr) {

  *os_ptr << "<DataArray ";

  WriteVtkXmlScalarHeader(variable_name, os_ptr);

  *os_ptr << ">\n";

  for (int i = 0; i < n; ++i)
    *os_ptr << v[i] << "\n";

  *os_ptr << "</DataArray>\n";

}

void WriteVtkXmlVectorHeader(const std::string& variable_name,
			     std::ostream* os_ptr) {

  *os_ptr << "type=\"Float64\" "
	  << "Name=\"" << variable_name << "\" "
	  << "NumberOfComponents=\"3\" "
	  << "format=\"ascii\"";

}

void WriteVtkXmlAsciiVector(int n, const RealType* vx, 
			    const RealType* vy, 
			    const RealType* vz, 
			    const std::string& variable_name,
			    std::ostream* os_ptr) {

  *os_ptr << "<DataArray ";

  WriteVtkXmlVectorHeader(variable_name, os_ptr);

  *os_ptr << ">\n";

  for (int i = 0; i < n; ++i) {

    const RealType z = (vz == NULL ? 0.0 : vz[i]);

    *os_ptr << vx[i] << " " << vy[i] << " " << z << "\n";

  }

  *os_ptr << "</DataArray>\n";

}



