add_executable(mkl_spmv_csr.exe mkl_spmv_csr.c)
target_link_libraries (mkl_spmv_csr.exe PUBLIC sparsematrix BLAS::BLAS )
target_include_directories(mkl_spmv_csr.exe PUBLIC ${CMAKE_SOURCE_DIR}/SparseMatrix)

add_executable(mkl_spmv_coo.exe mkl_spmv_coo.c)
target_link_libraries (mkl_spmv_coo.exe PUBLIC sparsematrix BLAS::BLAS )
target_include_directories(mkl_spmv_coo.exe PUBLIC ${CMAKE_SOURCE_DIR}/SparseMatrix)

