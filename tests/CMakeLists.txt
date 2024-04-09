add_executable(lp_seq_test main.cpp)
target_link_libraries(lp_seq_test linear_programming_seq)

add_executable(lp_test main_mpi.cpp)
target_include_directories(lp_test PRIVATE ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(lp_test ${MPI_CXX_LIBRARIES})