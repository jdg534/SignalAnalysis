cmake_minimum_required(VERSION 3.10.0)

file(GLOB_RECURSE SignalAnalysisSource *.cpp)
file(GLOB_RECURSE SignalAnalysisHeaders *.h)
add_library(SignalAnalysis ${SignalAnalysisHeaders} ${SignalAnalysisSource})

set(SignalAnalysisHeaderDir ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "" FORCE)
