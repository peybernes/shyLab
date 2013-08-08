#include "mpi_boundary_condition.hpp"

#include <cassert>
#include <iostream>

#include "mesh.hpp"
#include "variable_metadata.hpp"
#include "variable_store.hpp"

MpiBoundaryCondition::MpiBoundaryCondition() {
}

void MpiBoundaryCondition::AllocDataBuffers() {
  
  data_to_send.assign(NB_CELL_VALUES * send_cell_ids.size(), 0.0);
  data_recieved.assign(NB_CELL_VALUES * recieve_cell_ids.size(), 0.0);
  
}

void MpiBoundaryCondition::DeallocDataBuffers() {

  // Nothing to do here (thanks to RAII, since we use std::vectors to
  // store our data).

}

void MpiBoundaryCondition::FillSendDataBuffer(const VariableStore& var) {
  
  assert(NB_CELL_VALUES == var.nb_variables());

  const int nb_cells = var.nb_elements();

  const int size = send_cell_ids.size();

  assert(size <= nb_cells);

  for (int id_variable = 0; id_variable < NB_CELL_VALUES; ++id_variable)
    
    if (cell_variable_attributes[id_variable] & COMMUNICATED) {

      for (int i = 0; i < size; ++i) {
      
	const int id = send_cell_ids.at(i);

	data_to_send.at(id_variable * size + i) = var(id_variable)[id];

      }
    }
}
  

void MpiBoundaryCondition::UpdateVariables(VariableStore* var_ptr)
  const {

  assert(NB_CELL_VALUES == var_ptr->nb_variables());

  const int nb_cells = var_ptr->nb_elements();

  const int size = recieve_cell_ids.size();

  assert(size < nb_cells);

  for (int id_variable = 0; id_variable < NB_CELL_VALUES; ++id_variable)

    if (cell_variable_attributes[id_variable] & COMMUNICATED) {
    
      for (int i = 0; i < size; ++i) {
      
	const int id = recieve_cell_ids.at(i);

	(*var_ptr)(id_variable)[id]
	  = data_recieved.at(id_variable * size + i);

      }

    }
}

void MpiBoundaryCondition::Print(std::ostream* os_ptr) const {

  *os_ptr << "send cell ids size: " 
	  << send_cell_ids.size() 
	  << "\nsend cell ids: [ ";

    for (int i = 0; i < send_cell_ids.size(); ++i)
      *os_ptr << send_cell_ids.at(i) << " ";

  *os_ptr << " ]\n\n";

  *os_ptr << "recieve cell ids size: " 
	  << recieve_cell_ids.size() 
	  << "\nrecieve cell ids: [ ";

    for (int i = 0; i < recieve_cell_ids.size(); ++i)
      *os_ptr << recieve_cell_ids.at(i) << " ";

  *os_ptr << " ]\n";


}

void MpiBoundaryCondition::SwitchToLocalNumbering(const Mesh& mesh) {

  for (int i = 0; i < send_cell_ids.size(); ++i) {

    const int send_cell_global_id = send_cell_ids.at(i);

    if (mesh.unique_to_local_cell_ids().count(send_cell_global_id) == 1) {

      send_cell_ids.at(i) = mesh.unique_to_local_cell_ids().at(send_cell_global_id);

    }

    else {

      assert(0);

    }

  }

  for (int i = 0; i < recieve_cell_ids.size(); ++i) {

    const int recieve_cell_global_id = recieve_cell_ids.at(i);

    if (mesh.unique_to_local_cell_ids().count(recieve_cell_global_id) ==1) {

      recieve_cell_ids.at(i) = mesh.unique_to_local_cell_ids().at(recieve_cell_global_id);

    }

    else {

      assert(0);

    }

  }

}

#ifdef HAVE_MPI
void CommunicateAllMpiBC(MPI_Comm communicator, 
			 std::vector<MpiBoundaryCondition>* bc_ptr,
			 std::vector<MPI_Request>* send_handles_ptr,
			 std::vector<MPI_Request>* recv_handles_ptr) {

  int nb_processes = -1;
  int process_rank = -1;

  MPI_Comm_size(communicator, &nb_processes);
  MPI_Comm_rank(communicator, &process_rank);

  assert(nb_processes == bc_ptr->size());
  assert(send_handles_ptr->size() == nb_processes);
  assert(recv_handles_ptr->size() == nb_processes);

  // Prepost MPI recieves.
  for (int i = 0; i < nb_processes; ++i) {

    const int recieve_count = bc_ptr->at(i).data_recieved.size();

    void* recieve_buffer = 
      (recieve_count == 0 ? NULL : &(bc_ptr->at(i).data_recieved.at(0)));

    const int recieve_tag = i;

    MPI_Irecv(recieve_buffer, recieve_count, MPI_DOUBLE, 
	      i, recieve_tag, communicator, 
	      &(recv_handles_ptr->at(i)));

  }

  // If I understand well, we should interleave computations between
  // the recieves and sends (maybe mpi buffer data packing is
  // sufficient ?), *and also* between the send and waitall. To be
  // tested...
  for (int i = 0; i < nb_processes; ++i) {
    
    const int send_count = bc_ptr->at(i).data_to_send.size();

    void* send_buffer =
      (send_count == 0 ? NULL : &(bc_ptr->at(i).data_to_send.at(0)));

    const int send_tag = process_rank;
    
    MPI_Isend(send_buffer, send_count, MPI_DOUBLE, i, send_tag,
	      communicator, &(send_handles_ptr->at(i)));

  }

}

void FinishAllMpiBC(int nb_communications,
		    const std::vector<MPI_Request>& send_handles,
		    const std::vector<MPI_Request>& recv_handles) {

  assert(0 < nb_communications);
  assert(send_handles.size() == nb_communications);
  assert(recv_handles.size() == nb_communications);

  std::vector<MPI_Status> are_sends_finished(nb_communications);
  std::vector<MPI_Status> are_recieves_finished(nb_communications);

  // Casts deemed necessary to respect MPI API.

  MPI_Request* send_requests = 
    const_cast<MPI_Request*>(&(send_handles.at(0)));
  
  MPI_Waitall(nb_communications, 
	      send_requests,
  	      &(are_sends_finished.at(0)));

  MPI_Request* recv_requests = 
    const_cast<MPI_Request*>(&(recv_handles.at(0)));

  MPI_Waitall(nb_communications, 
	      recv_requests,
  	      &(are_recieves_finished.at(0)));

  
}
#endif // HAVE_MPI
