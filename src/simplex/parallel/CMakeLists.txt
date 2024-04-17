add_library(lp_par
    lp_par.h
    lp_par1.cpp
    lp_factory.cpp
)
target_link_libraries(lp_par ${MPI_CXX_LIBRARIES})
target_include_directories(lp_par PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${MPI_CXX_INCLUDE_PATH})