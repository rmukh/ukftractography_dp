#-----------------------------------------------------------------------------
set(MODULE_NAME UKFTractography)

#-----------------------------------------------------------------------------
# Check if float precision type is specified and use it if it is
if(UKF_USE_FLOAT)
  add_definitions("-DUKF_USE_FLOAT=1")
endif()

set(MODULE_INCLUDE_DIRECTORIES
  ${Teem_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../common
  ${CMAKE_CURRENT_SOURCE_DIR}/../ukf
  ${Eigen_INCLUDE_DIR}
  ${SphericalRidgelets_INCLUDE_DIR}
  )

set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${TEEM_LIB}
  ${ZLIB_LIBRARIES}
  ${VTK_LIBRARIES}
  ${Eigen_LIBRARIES}
  ${SphericalRidgelets_LIBRARY}
  ${OPENMP}
  UKFBase
  )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME UKFTractography
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ConvertVTK
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
add_subdirectory(Testing)
