add_executable(lpseq main.cpp)
target_link_libraries(lpseq linear_programming_seq)

add_executable(lppar main_mpi.cpp)
target_include_directories(lppar PRIVATE ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(lppar ${MPI_CXX_LIBRARIES})
target_link_libraries(lppar linear_programming_par)