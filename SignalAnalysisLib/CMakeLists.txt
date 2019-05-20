cmake_minimum_required(VERSION 3.10.0)

file(GLOB_RECURSE SignalAnalysisSource *.cpp)
file(GLOB_RECURSE SignalAnalysisHeaders *.h)
add_library(SignalAnalysis ${SignalAnalysisHeaders} ${SignalAnalysisSource})