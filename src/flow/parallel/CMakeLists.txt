add_library(flow_par flow_par.cpp)
target_link_libraries(flow_par lp_par1)
target_link_libraries(flow_par ${MPI_CXX_LIBRARIES})
target_include_directories(flow_par PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MPI_CXX_INCLUDE_PATH})