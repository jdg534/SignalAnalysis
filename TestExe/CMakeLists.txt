cmake_minimum_required(VERSION 3.10.0)

file(GLOB_RECURSE SignalAnalysisTestExeSource *.h *.cpp)
add_executable(SignalAnalysisLibTestExe ${SignalAnalysisTestExeSource})

target_link_libraries(SignalAnalysisLibTestExe SignalAnalysis)

target_include_directories(SignalAnalysisLibTestExe PUBLIC ${SignalAnalysisHeaderDir})
