set(proj zlib)

# Set dependency list
set(${proj}_DEPENDENCIES "")

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

# Fix for CMP0074
if (POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  unset(zlib_DIR CACHE)
  find_package(ZLIB REQUIRED)
  set(ZLIB_INCLUDE_DIR ${ZLIB_INCLUDE_DIRS})
  set(ZLIB_LIBRARY ${ZLIB_LIBRARIES})
endif()

# Sanity checks
if(DEFINED zlib_DIR AND NOT EXISTS ${zlib_DIR})
  message(FATAL_ERROR "zlib_DIR variable is defined but corresponds to nonexistent directory")
endif()

if(NOT DEFINED zlib_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

  set(EP_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
  set(EP_BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
  set(EP_INSTALL_DIR ${CMAKE_BINARY_DIR}/${proj}-install)

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY
    "https://github.com/madler/zlib.git"
    QUIET
    )

  ExternalProject_SetIfNotDefined(
    ${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG
    "21767c654d31d2dccdde4330529775c6c5fd5389"
    QUIET
    )

  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY "${${CMAKE_PROJECT_NAME}_${proj}_GIT_REPOSITORY}"
    GIT_TAG "${${CMAKE_PROJECT_NAME}_${proj}_GIT_TAG}"
    SOURCE_DIR ${EP_SOURCE_DIR}
    BINARY_DIR ${EP_BINARY_DIR}
    INSTALL_DIR ${EP_INSTALL_DIR}
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CMAKE_CACHE_ARGS
      ## CXX should not be needed, but it a cmake default test
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
#      -DZLIB_MANGLE_PREFIX:STRING=slicer_zlib_
      -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
    )

  ExternalProject_GenerateProjectDescription_Step(${proj})

  set(zlib_DIR ${EP_INSTALL_DIR} CACHE PATH "zlib dir" FORCE)
  set(ZLIB_ROOT ${zlib_DIR} CACHE PATH "zlib root" FORCE)
  set(ZLIB_INCLUDE_DIR ${zlib_DIR}/include CACHE PATH "zlib include dir" FORCE)
  if(WIN32)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") #to make compatible with msys2 gcc build
      set(ZLIB_LIBRARY ${zlib_DIR}/lib/libzlib.a CACHE FILEPATH "zlib library" FORCE)
    else()
      set(ZLIB_LIBRARY ${zlib_DIR}/lib/zlib.lib CACHE FILEPATH "zlib library" FORCE)
    endif()
  else()
    set(ZLIB_LIBRARY ${zlib_DIR}/lib/libz.a CACHE FILEPATH "zlib library" FORCE)
  endif()
else()
  # The project is provided using zlib_DIR, nevertheless since other project may depend on zlib,
  # let's add an 'empty' one
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDENCIES})
endif()

mark_as_superbuild(
  VARS
    ZLIB_INCLUDE_DIR:PATH
    ZLIB_LIBRARY:FILEPATH
    ZLIB_ROOT:PATH
  LABELS "FIND_PACKAGE"
  )

ExternalProject_Message(${proj} "ZLIB_INCLUDE_DIR:${ZLIB_INCLUDE_DIR}")
ExternalProject_Message(${proj} "ZLIB_LIBRARY:${ZLIB_LIBRARY}")
if(ZLIB_ROOT)
  ExternalProject_Message(${proj} "ZLIB_ROOT:${ZLIB_ROOT}")
endif()
