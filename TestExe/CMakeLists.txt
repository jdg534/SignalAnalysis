cmake_minimum_required(VERSION 3.10.0)

file(GLOB_RECURSE TestExeSource *.h *.cpp)
add_executable(TestExe ${TestExeSource})