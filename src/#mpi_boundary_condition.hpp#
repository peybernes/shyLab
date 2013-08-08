#ifndef BOUNDARY_CONDITION_HPP
#define BOUNDARY_CONDITION_HPP

#ifdef HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI

#include <iosfwd>
#include <vector>

#include "variable_store.hpp"

class VariableStore;
class Mesh;

class MpiBoundaryCondition {
public:
  MpiBoundaryCondition();
  void AllocDataBuffers();
  void DeallocDataBuffers();
  void FillSendDataBuffer(const VariableStore& var);
  void UpdateVariables(VariableStore* var_ptr) const;
  void Print(std::ostream* os_ptr) const;
  void SwitchToLocalNumbering(const Mesh& mesh);
  std::vector<int> send_cell_ids;
  std::vector<int> recieve_cell_ids;
  std::vector<RealType> data_to_send;
  std::vector<RealType> data_recieved;
private:
};

#ifdef HAVE_MPI
void CommunicateAllMpiBC(MPI_Comm communicator, 
			 std::vector<MpiBoundaryCondition>* bc_ptr,
			 std::vector<MPI_Request>* sends_handles_ptr,
			 std::vector<MPI_Request>* recvs_handles_ptr);

void FinishAllMpiBC(int nb_communications,
		    const std::vector<MPI_Request>& send_handles,
		    const std::vector<MPI_Request>& recv_handles);

#endif // HAVE_MPI

#endif // BOUNDARY_CONDITION_HPP
