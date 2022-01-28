add_executable(spmv_dev.exe SpMV_dev.c)
target_link_libraries( spmv_dev.exe PUBLIC sparsematrix )
target_include_directories(spmv_dev.exe PUBLIC ${CMAKE_SOURCE_DIR}/SparseMatrix)

add_executable(spmv.exe SpMV.c)
target_link_libraries( spmv.exe PUBLIC sparsematrix )
target_include_directories(spmv.exe PUBLIC ${CMAKE_SOURCE_DIR}/SparseMatrix)

