set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/SRC )

set( ALLAUX   
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/superlu_timer.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dclock.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sp_ienv.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/lsame.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/xerbla.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/util.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pmemory.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/qrnzcnt.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/await.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/get_perm_c.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/mmd.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/colamd.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sp_coletree.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_scheduler.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sp_colorder.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_mark_busy_descends.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_pruneL.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_super_bnd_dfs.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_relax_snode.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/heap_relax_snode.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_synch.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_finalize.c
     )
#SET_SOURCE_FILES_PROPERTIES(await.c PROPERTIES COMPILE_FLAGS -O0) #enable it for Cray
SET_SOURCE_FILES_PROPERTIES(superlu_timer.c PROPERTIES COMPILE_FLAGS -O0)

set( DZAUX     
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlamch.c
     )
SET_SOURCE_FILES_PROPERTIES(${DZAUX} PROPERTIES COMPILE_FLAGS -O0)

set( DLUSRC2  
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dreadhb.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dreadrb.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dmatgen.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgssv.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgssvx.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgstrs.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgsrfs.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgscon.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlacon.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlangs.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgsequ.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlaqgs.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dpivotgrowth.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdmemory.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdutil.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dmyblas2.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_init.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_thread.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_thread_init.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_thread_finalize.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_factor_snode.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_snode_dfs.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_snode_bmod.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_panel_dfs.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_panel_bmod.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_copy_to_ucol.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_pivotL.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_column_dfs.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_column_bmod.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_bmod1D.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_bmod2D.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_bmod1D_mv2.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_bmod2D_mv2.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsp_blas2.c
        ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsp_blas3.c
    )

add_definitions( -DUSE_VENDOR_BLAS )
add_definitions( -DAdd_ )
add_definitions( -D_LONGINT )
add_definitions( -w )
add_definitions( -D__OPENMP ) #be careful to choose the correct MACH for SuperLU_MT

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions( -DNO_TIMER )
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")


###############################################################
add_library(superlu_mt STATIC ${DLUSRC2} ${ALLAUX} ${DZAUX} )
target_link_libraries(superlu_mt ${BLAS_LIBRARIES} )
# install(TARGETS superlu DESTINATION libs )
