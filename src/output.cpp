#include "output.hpp"

#include <fstream>
#include <iomanip>
#include <limits>
#ifdef HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI
#include <string>

#include "array_io.hpp"
#include "cell_variable_metadata.hpp"
#include "serialize.hpp"
#include "simulation.hpp"
#include "variable_attribute.hpp"
#include "vertice_variable_metadata.hpp"

Output::Output(Simulation *sim_ptr, std::string &stream_name, const Timetable &timetable):
  Event(sim_ptr, timetable, POST), m_stream_name(stream_name) {}

OutputSimulation::OutputSimulation(Simulation *sim_ptr,
				   std::string &stream_name, 
				   const std::string& io_format_string,
				   const Timetable & timetable):

  Output(sim_ptr, stream_name, timetable) {
  
  m_name = "OutputSimulation";

  if (io_format_string == "binary")
    m_io_format = BINARY;

  // Default value.
  else
    m_io_format = ASCII;

}

void OutputSimulation::Execute() {
  
  int nb_processes = 1;
  int process_rank = 0;

  std::cerr << "Writing simulation results...\n";

#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &nb_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
#endif // HAVE_MPI

  std::stringstream meta_iter_stream;
  std::stringstream iter_stream;
  std::stringstream suffix_iter_stream;

  meta_iter_stream  << "-" << std::setfill('0') << std::setw(4)
		    << m_sim_ptr->clock.iter()
		    << ".pvtu";
  
  iter_stream << "-PE" << std::setfill('0') << std::setw(4) 
	      << process_rank;

  iter_stream << "-" << std::setfill('0') << std::setw(4) 
	      << m_sim_ptr->clock.iter()
	      << ".vtr";
  
  const std::string iter_string = iter_stream.str();

  const std::string meta_iter_string = meta_iter_stream.str();

  // FIXME
  std::string filename =  m_stream_name;

  std::string meta_filename = m_stream_name;

  const std::string search_string = "\%i";

  std::string::size_type pos = 0;
  while ((pos = filename.find(search_string, pos)) != std::string::npos) {
    filename.replace(pos, search_string.size(), iter_string);
    ++pos;
  }
  
  pos = 0;
  while ((pos = meta_filename.find(search_string, pos)) != std::string::npos) {
    meta_filename.replace(pos, search_string.size(), meta_iter_string);
    ++pos;
  }

  if (process_rank == 0) {

    std::ofstream meta_stream(meta_filename.c_str());

    // Write header.
    meta_stream << "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type=\"PUnstructuredGrid\">\n";

    meta_stream << "<PUnstructuredGrid GhostLevel=\"2\">\n";

    meta_stream << "<PPointData>\n";
    meta_stream << "</PPointData>\n";

    meta_stream << "<PCellData>\n";

    for (VariableDatabase::const_iterator it = m_sim_ptr->variables_database.begin();
	 it != m_sim_ptr->variables_database.end(); ++it) {

      const std::string variable_name = it->first;
      const VariableEntry current_entry = it->second;
      
      if (current_entry.written() && current_entry.support() == CELL) {

	meta_stream << "<PDataArray ";
	WriteVtkXmlScalarHeader(variable_name, &meta_stream);
	meta_stream << "/>\n";

      }

    }

    for (int id_var = 0; id_var < NB_CELL_VECTOR_VALUES; ++id_var) {

      const int x_coord_id =
	cell_vector_variable_coordinates[id_var][0];
      const int y_coord_id = 
	cell_vector_variable_coordinates[id_var][1];

      if (cell_vector_variable_attributes[id_var] & WRITTEN) {

	meta_stream << "<PDataArray ";
	WriteVtkXmlVectorHeader(cell_vector_variable_names[id_var], &meta_stream);
	meta_stream << "/>\n";

      }      
    
    }
    
    meta_stream << "</PCellData>\n";

    meta_stream << "<PPoints>\n";
    
    meta_stream << "<PDataArray "
		<< "type=\"Float32\" "
		<< "NumberOfComponents=\"3\" "
		<< "format=\"ascii\"/>\n";
    
    meta_stream << "</PPoints>\n";

    for (int i = 0; i < nb_processes; ++i) {

      std::stringstream sstream;

      sstream << "-PE" << std::setfill('0') << std::setw(4) << i;

      sstream << "-" << std::setfill('0') << std::setw(4) 
	      << m_sim_ptr->clock.iter() << ".vtu";

      std::string suffix_string = sstream.str();
      
      std::string filename = "sim%i";

      std::string::size_type pos = 0;
      while ((pos = filename.find(search_string, pos)) != std::string::npos) {
	filename.replace(pos, search_string.size(), suffix_string);
	++pos;
      }

      meta_stream << "<Piece Source=\"";

      meta_stream << filename;
      
      meta_stream << "\"/>\n";

    }

    // Write footer.
    meta_stream << "</PUnstructuredGrid>\n";
    
    meta_stream << "</VTKFile>\n";

    meta_stream.close();

  }

  // // Now, dump mesh faces data.

  // std::stringstream meta_iter_stream_faces;
  // std::stringstream iter_stream_faces;
  // std::stringstream suffix_iter_stream_faces;

  // meta_iter_stream_faces  << "-" << std::setfill('0') << std::setw(4)
  // 			  << m_sim_ptr->clock.iter()
  // 			  << ".faces.pvtu";
  
  // iter_stream_faces << "-PE" << std::setfill('0') << std::setw(4) 
  // 		    << process_rank;

  // iter_stream_faces << "-" << std::setfill('0') << std::setw(4) 
  // 		    << m_sim_ptr->clock.iter()
  // 		    << ".faces.vtu";
  
  // const std::string iter_string_faces = iter_stream_faces.str();

  // const std::string meta_iter_string_faces = meta_iter_stream_faces.str();

  // // FIXME
  // std::string filename_faces =  m_stream_name + ".faces";

  // std::string meta_filename_faces = m_stream_name + ".faces";

  // const std::string search_string_faces = "\%i";

  // std::string::size_type pos_faces = 0;
  // while ((pos_faces = filename_faces.find(search_string, pos_faces)) 
  // 	 != std::string::npos) {
  //   filename_faces.replace(pos_faces, search_string_faces.size(), iter_string_faces);
  //   ++pos_faces;
  // }
  
  // pos_faces = 0;
  // while ((pos = meta_filename_faces.find(search_string_faces,
  // 					 pos_faces)) 
  // 	 != std::string::npos) {
  //   meta_filename_faces.replace(pos_faces, search_string_faces.size(), meta_iter_string_faces);
  //   ++pos_faces;
  // }

  //   if (process_rank == 0) {

  //     std::ofstream meta_stream_faces(meta_filename.c_str());

  //     // Write header.
  //     meta_stream_faces << "<?xml version=\"1.0\"?>\n"
  // 		  << "<VTKFile type=\"PUnstructuredGrid\">\n";

  //     meta_stream_faces << "<PUnstructuredGrid GhostLevel=\"2\">\n";

  //     meta_stream_faces << "<PPointData>\n";
  //     meta_stream_faces << "</PPointData>\n";

  //     meta_stream_faces << "<PCellData>\n";

  //     for (int id_var = 0; id_var < NB_CELL_VALUES; ++id_var) {

  // 	if (cell_variable_attributes[id_var] & WRITTEN) {

  // 	  meta_stream_faces << "<PDataArray ";
  // 	  WriteVtkXmlScalarHeader(cell_variable_names[id_var], &meta_stream_faces);
  // 	  meta_stream_faces << "/>\n";

  // 	}

  //     }

  //     for (int id_var = 0; id_var < NB_CELL_VECTOR_VALUES; ++id_var) {

  // 	const int x_coord_id =
  // 	  cell_vector_variable_coordinates[id_var][0];
  // 	const int y_coord_id = 
  // 	  cell_vector_variable_coordinates[id_var][1];

  // 	if (cell_vector_variable_attributes[id_var] & WRITTEN) {

  // 	  meta_stream_faces << "<PDataArray ";
  // 	  WriteVtkXmlVectorHeader(cell_vector_variable_names[id_var], &meta_stream_faces);
  // 	  meta_stream_faces << "/>\n";

  // 	}      
    
  //     }
    
  //     meta_stream_faces << "</PCellData>\n";

  //     meta_stream_faces << "<PPoints>\n";
    
  //     meta_stream_faces << "<PDataArray "
  // 		  << "type=\"Float32\" "
  // 		  << "NumberOfComponents=\"3\" "
  // 		  << "format=\"ascii\"/>\n";
    
  //     meta_stream_faces << "</PPoints>\n";

  //     for (int i = 0; i < nb_processes; ++i) {

  // 	std::stringstream sstream_faces;

  // 	sstream_faces << "-PE" << std::setfill('0') << std::setw(4) << i;

  // 	sstream_faces << "-" << std::setfill('0') << std::setw(4) 
  // 		<< m_sim_ptr->clock.iter() << ".vtu";

  // 	std::string suffix_string = sstream_faces.str();
      
  // 	std::string filename = "sim%i";

  // 	std::string::size_type pos = 0;
  // 	while ((pos = filename.find(search_string, pos)) != std::string::npos) {
  // 	  filename.replace(pos, search_string.size(), suffix_string);
  // 	  ++pos;
  // 	}

  // 	meta_stream_faces << "<Piece Source=\"";

  // 	meta_stream_faces << filename;
      
  // 	meta_stream_faces << "\"/>\n";

  //     }

  //     // Write footer.
  //     meta_stream_faces << "</PUnstructuredGrid>\n";
    
  //     meta_stream_faces << "</VTKFile>\n";

  //     meta_stream_faces.close();

  //   }
  // }

  std::ofstream stream(filename.c_str());

  if (!stream.good())
    std::cerr << "\nWARNING: could not write to file " << filename.c_str() << "\n";

  else {

    m_sim_ptr->m_grid.WriteHeaderVTKXmlAscii(&stream);
    m_sim_ptr->m_grid.WriteVTKXmlAscii(&stream);

    // m_sim_ptr->mesh.WriteHeaderVTKXmlAscii(&stream);
    // m_sim_ptr->mesh.WriteVTKXmlAscii(&stream);
    
    stream << "<CellData>\n";

    const int nb_cells = m_sim_ptr->cell_variables.nb_elements();

    assert(nb_cells == m_sim_ptr->m_grid.nb_cells());

    for (VariableDatabase::const_iterator it = m_sim_ptr->variables_database.begin();
	 it != m_sim_ptr->variables_database.end(); ++it) {

      const std::string variable_name = it->first;
      const VariableEntry current_entry = it->second;
      
      if (current_entry.written() && current_entry.support() == CELL) {
    	WriteVtkXmlAsciiScalar(nb_cells, 
    			       m_sim_ptr->cell_variables(current_entry.id()),
    			       variable_name,
    			       &stream);

      }
    }

    for (int id_var = 0; id_var < NB_CELL_VECTOR_VALUES; ++id_var) {

      const int x_coord_id =
    	cell_vector_variable_coordinates[id_var][0];
      const int y_coord_id = 
    	cell_vector_variable_coordinates[id_var][1];
      const int z_coord_id = 
    	cell_vector_variable_coordinates[id_var][2];

      RealType* z_coords = NULL;

      if (z_coord_id >= 0)
  	z_coords = m_sim_ptr->cell_variables(z_coord_id);

      if (cell_vector_variable_attributes[id_var] & WRITTEN)
    	WriteVtkXmlAsciiVector(nb_cells, 
    			       m_sim_ptr->cell_variables(x_coord_id),
    			       m_sim_ptr->cell_variables(y_coord_id),
  			       z_coords,
    			       cell_vector_variable_names[id_var], 
    			       &stream);
    
    }

    stream << "</CellData>\n";

    stream << "<PointData>\n";

    const int nb_vertices = m_sim_ptr->vertice_variables.nb_elements();

    assert(nb_vertices == m_sim_ptr->m_grid.nb_vertices());

    for (VariableDatabase::const_iterator it = m_sim_ptr->variables_database.begin();
	 it != m_sim_ptr->variables_database.end(); ++it) {

      const std::string variable_name = it->first;
      const VariableEntry current_entry = it->second;
      
      if (current_entry.written() && current_entry.support() == VERTICE) {
    	WriteVtkXmlAsciiScalar(nb_vertices, 
    			       m_sim_ptr->vertice_variables(current_entry.id()),
			       variable_name,
    			       &stream);

      }
    }

    for (int id_var = 0; id_var < NB_VERTICE_VECTOR_VALUES; ++id_var) {

      const int x_coord_id =
    	vertice_vector_variable_coordinates[id_var][0];
      const int y_coord_id = 
    	vertice_vector_variable_coordinates[id_var][1];
      const int z_coord_id = 
    	vertice_vector_variable_coordinates[id_var][2];

      RealType* z_coords = NULL;

      if (z_coord_id >= 0)
  	z_coords = m_sim_ptr->vertice_variables(z_coord_id);

      if (vertice_vector_variable_attributes[id_var] & WRITTEN)
    	WriteVtkXmlAsciiVector(nb_vertices, 
    			       m_sim_ptr->vertice_variables(x_coord_id),
    			       m_sim_ptr->vertice_variables(y_coord_id),
  			       z_coords,
    			       vertice_vector_variable_names[id_var], 
    			       &stream);
    
    }


    stream << "</PointData>\n";

    m_sim_ptr->m_grid.WriteFooterVTKXmlAscii(&stream);
    //m_sim_ptr->mesh.WriteFooterVTKXmlAscii(&stream);

    stream.close();
  }
}

void Output::Save(ptree &pt) {
  
  ptree pt_timetable;
  m_timetable.Save(pt_timetable);

  pt.put("filename", m_stream_name);
  pt.put("format", m_io_format);
  pt.put_child("Timetable", pt_timetable);
}

OutputLocation::OutputLocation(Simulation *sim_ptr, std::string &stream_name, 
			       const Timetable &timetable,
			       const RealType x, const RealType y):
  Output(sim_ptr, stream_name, timetable), m_x(x), m_y(y), m_cell_id(-1) {};

void OutputLocation::Execute() {

  // erase the file if it already exists during first simulation
  // iteration.
  if (m_timetable.istart() == m_sim_ptr->clock.iter()) {

    std::ofstream stream(m_stream_name.c_str());

    if (!stream.good())
      std::cerr << "\nWARNING: could not write to file " << m_stream_name.c_str() << "\n";

    else {
      stream << "# Shy output at location (x, y) = (" 
	     << m_x << ", " << m_y << ")\n"
	     << "# Time ";

      for (int id_var = 0; id_var < NB_CELL_VALUES; ++id_var)
	stream << cell_variable_names[id_var] << " ";

      stream << "\n";
      stream.close();
    }

    m_cell_id = m_sim_ptr->mesh.TriangleIndex(m_x, m_y, 0.0);
    assert((0 <= m_cell_id) && (m_cell_id < m_sim_ptr->mesh.nb_cells()));

  }

  std::ofstream stream(m_stream_name.c_str(), std::ofstream::app);

  if (!stream.good())
    std::cerr << "\nWARNING: could not write to file " << m_stream_name.c_str() << "\n";

  else {

    stream << m_sim_ptr->clock.time() << " ";

    for (int id_var = 0; id_var < NB_CELL_VALUES; ++id_var) {

      if (cell_variable_attributes[id_var] & WRITTEN)
	stream << m_sim_ptr->cell_variables(id_var)[m_cell_id]
	       << " ";

    }
    
    stream.close();
  }
}

OutputMinMax::OutputMinMax(Simulation* sim_ptr, std::string &stream_name, 
			   const Timetable & timetable):
  Output(sim_ptr, stream_name, timetable) {
  m_name = "OutputMinMax";
}

void OutputMinMax::Execute() {

  const RealType maxfloat = std::numeric_limits<RealType>::max();
  const RealType minfloat = std::numeric_limits<RealType>::min();

  RealType in_H_min, in_U_min, in_V_min, in_Z_min = maxfloat;
  RealType in_H_max, in_U_max, in_V_max, in_Z_max = minfloat;

  RealType out_H_min, out_U_min, out_V_min, out_Z_min = maxfloat;
  RealType out_H_max, out_U_max, out_V_max, out_Z_max = minfloat;

  RealType face_fluxes_H_min, face_fluxes_HU_min, 
    face_fluxes_HV_min = maxfloat;

  RealType face_fluxes_H_max, face_fluxes_HU_max, 
    face_fluxes_HV_max = minfloat;
  
  // // Compute min/max for cell values.
  // for (int i = 0; i < m_sim_ptr->var.nb_cells; ++i) {

  //   // in_H, out_H
  //   in_H_min = std::min(m_sim_ptr->var.in_H[i], 
  // 			in_H_min);

  //   in_H_max = std::max(m_sim_ptr->var.in_H[i], 
  // 			in_H_max);

  //   out_H_min = std::min(m_sim_ptr->var.out_H[i], 
  // 			 out_H_min);

  //   out_H_max = std::max(m_sim_ptr->var.out_H[i], 
  // 			 out_H_max);

  //   // in_U, out_U
  //   in_U_min = std::min(m_sim_ptr->var.in_U[i], 
  // 			in_U_min);

  //   in_U_max = std::max(m_sim_ptr->var.in_U[i], 
  // 			in_U_max);

  //   out_U_min = std::min(m_sim_ptr->var.out_U[i], 
  // 			 out_U_min);

  //   out_U_max = std::max(m_sim_ptr->var.out_U[i], 
  // 			 out_U_max);
    
  //   // in_V, out_V
  //   in_V_min = std::min(m_sim_ptr->var.in_V[i], 
  // 			in_V_min);

  //   in_V_max = std::max(m_sim_ptr->var.in_V[i], 
  // 			in_V_max);

  //   out_V_min = std::min(m_sim_ptr->var.out_V[i], 
  // 			 out_V_min);

  //   out_V_max = std::max(m_sim_ptr->var.out_V[i], 
  // 			 out_V_max);
  
  // }

  // for (int i = 0; i < m_sim_ptr->var.nb_faces; ++i) {

  //   // face fluxes.
  //   face_fluxes_H_min = std::min(m_sim_ptr->var.face_fluxes_H[i], 
  // 				 face_fluxes_H_min);

  //   face_fluxes_H_max = std::max(m_sim_ptr->var.face_fluxes_H[i], 
  // 				 face_fluxes_H_max);

  //   face_fluxes_HU_min = std::min(m_sim_ptr->var.face_fluxes_HU[i], 
  // 				  face_fluxes_HU_min);

  //   face_fluxes_HU_max = std::max(m_sim_ptr->var.face_fluxes_HU[i], 
  // 				  face_fluxes_HU_max);

  //   face_fluxes_HV_min = std::min(m_sim_ptr->var.face_fluxes_HV[i], 
  // 				  face_fluxes_HV_min);

  //   face_fluxes_HV_max = std::max(m_sim_ptr->var.face_fluxes_HV[i], 
  // 				  face_fluxes_HV_max);
  // }

  std::cerr << std::fixed << std::setprecision(10);
  std::cerr << "# -------------------------------------------\n";
  std::cerr << "#           Values min/max\n";
  std::cerr << "# -------------------------------------------\n";
  std::cerr << "in_H  " << in_H_min << " " << in_H_max << " || ";
  std::cerr << "out_H " << out_H_min << " " << out_H_max << "\n";
  std::cerr << "in_U  " << in_U_min << " " << in_U_max << " || ";
  std::cerr << "out_U " << out_U_min << " " << out_U_max << "\n";
  std::cerr << "in_V  " << in_V_min << " " << in_V_max << " || ";
  std::cerr << "out_V " << out_V_min << " " << out_V_max << "\n";
  std::cerr << "flux_H " << face_fluxes_H_min 
	    << " " << face_fluxes_H_max << "\n";
  std::cerr << "flux_U " << face_fluxes_HU_min 
	    << " " << face_fluxes_HU_max << "\n";
  std::cerr << "flux_V " << face_fluxes_HV_min 
	    << " " << face_fluxes_HV_max << "\n";
  std::cerr << "# -------------------------------------------\n";

}

OutputTime::OutputTime(Simulation *sim_ptr, std::string &stream_name, 
		       const Timetable &timetable):
  Output(sim_ptr, stream_name, timetable) {};

void OutputTime::Execute() {
  std::cerr << "Iteration " << m_sim_ptr->clock.iter()
	    << ", time: " << m_sim_ptr->clock.time() << std::endl;
}

OutputNorms::OutputNorms(Simulation *sim_ptr, 
			 std::string &stream_name, 
			 const Timetable &timetable,
			 const std::string& name1,
			 const std::string& name2):

  Output(sim_ptr, stream_name, timetable) {
  
  m_name1 = name1;
  m_name2 = name2;

};

void OutputNorms::Execute() {


  // erase the file if it already exists during first simulation
  // iteration.
  if (m_timetable.istart() == m_sim_ptr->clock.iter()) {

    std::ofstream stream(m_stream_name.c_str());

    if (!stream.good())
      std::cerr << "\nWARNING: could not write to file " << m_stream_name.c_str() << "\n";

    else {
      stream << "# Shy | norms between variables " 
	     << m_name1 << " and " << m_name2 << "\n";

    }

  }


  var1 = m_sim_ptr->cell_variables.GetVariable(m_sim_ptr->variables_database, m_name1);

  var2 = m_sim_ptr->cell_variables.GetVariable(m_sim_ptr->variables_database, m_name2);

  std::string cell_volumes_str = "cell_volumes";

  RealType* cell_volumes = m_sim_ptr->cell_variables.GetVariable(m_sim_ptr->variables_database, cell_volumes_str);

  const int nb_cells = m_sim_ptr->m_grid.nb_cells();

  RealType l1_norm = 0.0;
  RealType l2_norm = 0.0;
  RealType l_infinite_norm = 0.0;

  // que pour les mailles internes...

  const int nx = m_sim_ptr->m_grid.nx();
  const int ny = m_sim_ptr->m_grid.ny();

  for (int j = 2; j < ny - 2; ++j) {
    for (int i = 2; i < nx - 2; ++i) {
      
      const int oo = (nx * j) + i;

      l1_norm += fabs(var1[oo] - var2[oo]) * cell_volumes[oo];
      l2_norm += (var1[oo] - var2[oo]) * (var1[oo] - var2[oo]) * cell_volumes[oo];
      l_infinite_norm = std::max(l_infinite_norm, std::abs(var1[oo] - var2[oo]));

    }

  }

  std::cerr << "l1_norm=" << l1_norm << "\n";

  l2_norm = sqrt(l2_norm);

  std::ofstream stream(m_stream_name.c_str(), std::ofstream::app);

  if (!stream.good())
    std::cerr << "\nWARNING: could not write to file " << m_stream_name.c_str() << "\n";
  
  else {

    stream << m_sim_ptr->clock.time() << " ";

    stream << m_sim_ptr->clock.iter() << " ";

    stream << l1_norm << " " << l2_norm << " " << l_infinite_norm << "\n";
    
    stream.close();
  }

}

void OutputNorms::Save(ptree &pt) {
  
  ptree pt_timetable;
  m_timetable.Save(pt_timetable);

  pt.put("filename", m_stream_name);
  pt.put("format", m_io_format);
  pt.put_child("Timetable", pt_timetable);
  pt.put("variable1", m_name1);
  pt.put("variable2", m_name2);
}

 OutputCartesianLine::OutputCartesianLine(Simulation *sim_ptr, 
					  std::string &stream_name, 
					  const Timetable &timetable,
					  const std::string& direction,
					  int id):
   Output(sim_ptr, stream_name, timetable), m_indice(id) {
   
   const int nx = m_sim_ptr->m_grid.nx();
   const int ny = m_sim_ptr->m_grid.ny();

   if ((direction == "x") || (direction == "X")) {
     
     m_direction = X;

     std::cerr << "m_indice=" << m_indice << "\n";

     if ((m_indice < 0) || (m_indice >= nx)) {
       
       std::cerr << "ERROR: in OutputCartesianLine, the number of the line is not between 0 and nx\n";
       assert(0);

     }


   } else if ((direction == "y") || (direction == "Y")) {

     m_direction = Y;

     if ((m_indice < 0) || (m_indice >= ny)) {
       
       std::cerr << "ERROR: in OutputCartesianLine, the number of the line is not between 0 and ny\n";
       assert(0);

     }


   } else {

     std::cerr << "ERROR: in OutputCartesianLine, expected \"x\", \"X\", \"y\", \"Y\" for direction field, got " << direction << "\n";
     assert(0);

   }

 };

void OutputCartesianLine::Execute() {

  std::string filename =  m_stream_name;

  std::stringstream iter_stream;

  iter_stream << std::setfill('0') << std::setw(4) 
	      << m_sim_ptr->clock.iter();
  
  const std::string iter_string = iter_stream.str();
  
  const std::string search_string = "\%i";

  std::string::size_type pos = 0;
  while ((pos = filename.find(search_string, pos)) != std::string::npos) {
    filename.replace(pos, search_string.size(), iter_string);
    ++pos;
  }

  filename += ".txt";

  std::ofstream stream(filename.c_str());


  if (!stream.good())
    std::cerr << "\nWARNING: could not write to file " << m_stream_name.c_str() << "\n";
  
  else {

    stream << "# ";

    for (int id_var = 0; id_var < NB_CELL_VALUES; ++id_var) {
      
      if (cell_variable_attributes[id_var] & WRITTEN)
	stream << cell_variable_names[id_var] << " ";

    }
      
    stream << "\n\n";

    const int nx = m_sim_ptr->m_grid.nx();
    const int ny = m_sim_ptr->m_grid.ny();

    if (m_direction == X) {

      for (int i = 0; i < nx; ++i) {

	const int id = (nx * m_indice) + i;
	
	for (int id_var = 0; id_var < NB_CELL_VALUES; ++id_var) {

	  if (cell_variable_attributes[id_var] & WRITTEN) {
	    
	    stream << m_sim_ptr->cell_variables(id_var)[id] << " ";
	    
	  }
	  
	}
	
	stream << "\n";

      }

    } else if (m_direction == Y) {

      for (int j = 0; j < ny; ++j) {

	const int id = (nx * j) + m_indice;

	for (int id_var = 0; id_var < NB_CELL_VALUES; ++id_var) {

	  if (cell_variable_attributes[id_var] & WRITTEN) {

	    stream << m_sim_ptr->cell_variables(id_var)[id] << " ";
	    
	  }
	
	}

	stream << "\n";

      }

      
    } else {

      assert(0);

    }

    stream.close();
  }

}

void OutputCartesianLine::Save(ptree& pt) {

  //assert(0);

}
