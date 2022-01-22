add_executable(SpMV.exe SpMV.c)
target_link_libraries (SpMV.exe PUBLIC sparsematrix)
target_include_directories(SpMV.exe PUBLIC ${CMAKE_SOURCE_DIR}/SparseMatrix)

