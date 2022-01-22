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

target_include_directories(sparsematrix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})