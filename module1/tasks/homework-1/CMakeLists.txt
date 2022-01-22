add_executable(spmv.exe SpMV.c)
target_link_libraries( spmv.exe PUBLIC sparsematrix )
target_include_directories(spmv.exe PUBLIC ${CMAKE_SOURCE_DIR}/SparseMatrix)

