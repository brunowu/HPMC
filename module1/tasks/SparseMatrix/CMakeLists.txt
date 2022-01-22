set(SOURCE
	mmio.c
	loadmm.c
	cmp.c
	dense.c
	coo.c
	csr.c
	ell.c
	dia.c
)

add_library(sparsematrix ${SOURCE})
if(OpenMP_FOUND)
	target_compile_definitions(sparsematrix PRIVATE USE_OPENMP=1)
	target_link_libraries( sparsematrix OpenMP::OpenMP_C)	
endif()

target_include_directories(sparsematrix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
