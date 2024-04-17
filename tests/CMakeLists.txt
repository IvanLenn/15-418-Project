add_executable(lpseq main.cpp)
target_link_libraries(lpseq lp_seq)

add_executable(flowseq flow.cpp)
target_link_libraries(flowseq lp_seq)
target_link_libraries(flowseq flow_seq)

add_executable(lppar main_mpi.cpp)
target_include_directories(lppar PRIVATE ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(lppar ${MPI_CXX_LIBRARIES})
target_link_libraries(lppar lp_par)