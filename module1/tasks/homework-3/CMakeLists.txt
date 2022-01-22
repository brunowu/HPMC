add_executable(cusparse_spmv_csr.exe cusparse_spmv_csr.c)
target_link_libraries(cusparse_spmv_csr.exe PRIVATE CUDA::cudart CUDA::cusparse sparsematrix)
target_include_directories(mkl_spmv_csr.exe PUBLIC ${CMAKE_SOURCE_DIR}/SparseMatrix)

add_executable(cusparse_spmv_coo.exe cusparse_spmv_coo.c)
target_link_libraries(cusparse_spmv_coo.exe PRIVATE CUDA::cudart CUDA::cusparse sparsematrix)
target_include_directories(mkl_spmv_coo.exe PUBLIC ${CMAKE_SOURCE_DIR}/SparseMatrix)
