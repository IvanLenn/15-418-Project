add_library(lp_par1 lp_par1.cpp)
target_link_libraries(lp_par1 ${MPI_CXX_LIBRARIES})
target_include_directories(lp_par1 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MPI_CXX_INCLUDE_PATH})