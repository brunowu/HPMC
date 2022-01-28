add_executable(getMemSize.exe getMemSize.c)
target_link_libraries (getMemSize.exe PUBLIC sparsematrix)
target_include_directories(getMemSize.exe PUBLIC ${CMAKE_SOURCE_DIR}/SparseMatrix)

