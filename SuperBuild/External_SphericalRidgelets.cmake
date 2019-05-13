
set(proj SphericalRidgelets)

# Set dependency list
set(${proj}_DEPENDENCIES ITK Eigen)

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  unset(SphericalRidgelets_DIR CACHE)
  find_package(SphericalRidgelets REQUIRED)
endif()

# Sanity checks
if(DEFINED SphericalRidgelets_DIR AND NOT EXISTS ${SphericalRidgelets_DIR})
  message(FATAL_ERROR "SphericalRidgelets_DIR variable is defined but corresponds to nonexistent directory")
endif()

if(NOT DEFINED SphericalRidgelets_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

  set(EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS)

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG
    "42b1b6255161b3cdc90de4421b1684647647b555"
  )
  
  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY 
    "${git_protocol}://github.com/research-enthusiast/spherical_ridgelets.git"
  )

  set(EP_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
  set(EP_BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
  set(EP_INSTALL_DIR ${CMAKE_BINARY_DIR}/${proj}-install)

  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY "${${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY}"
    GIT_TAG "${${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG}"
    SOURCE_DIR ${EP_SOURCE_DIR}
    BINARY_DIR ${EP_BINARY_DIR}
    INSTALL_DIR ${EP_INSTALL_DIR}
    CMAKE_CACHE_ARGS
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      # Not needed -DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
      ${CMAKE_PROJECT_INCLUDE_EXTERNAL_PROJECT_ARG}
      -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
      -DITK_DIR:PATH=${ITK_DIR}
      -DEigen3_DIR:PATH=${Eigen_DIR}
      -DSPH_DIR:PATH=${EP_SOURCE_DIR}/Ridgelets/SphericalRidgelets
      ${EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS}
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )

  set(SphericalRidgelets_DIR ${EP_BINARY_DIR})

  #-----------------------------------------------------------------------------
  # Launcher setting specific to build tree

  # library paths
  set(${proj}_LIBRARY_PATHS_LAUNCHER_BUILD ${SphericalRidgelets_DIR}/bin/<CMAKE_CFG_INTDIR>)
  mark_as_superbuild(
    VARS ${proj}_LIBRARY_PATHS_LAUNCHER_BUILD
    LABELS "LIBRARY_PATHS_LAUNCHER_BUILD" "PATHS_LAUNCHER_BUILD"
    )

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDENCIES})
endif()

mark_as_superbuild(
  VARS SphericalRidgelets_DIR:PATH
  LABELS "FIND_PACKAGE"
  )
