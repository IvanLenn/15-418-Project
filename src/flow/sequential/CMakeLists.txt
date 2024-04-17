add_library(flow_seq flow_seq.cpp)
target_link_libraries(flow_seq lp_seq)
target_include_directories(flow_seq PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})