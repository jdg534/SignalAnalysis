cmake_minimum_required(VERSION 3.10.0)

project (SignalAnalysis)

# Find gnu plot
# find_package(Gnuplot REQUIRED)
# not needed that was just reading output text file from the test exe

add_subdirectory(SignalAnalysisLib)
add_subdirectory(TestExe)
