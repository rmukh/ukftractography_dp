get_filename_component(INPUT ${CMAKE_CURRENT_SOURCE_DIR}/../Data/Input REALPATH)
get_filename_component(BASELINE ${CMAKE_CURRENT_SOURCE_DIR}/../Data/Baseline REALPATH)
set(TESTING_RESULTS_DIRECTORY "${CMAKE_BINARY_DIR}/Testing/Temporary")
file(MAKE_DIRECTORY ${TESTING_RESULTS_DIRECTORY}) ## Need to generate the directory for writing the temporary test outputs to

set(CLP UKFTractography)

include_directories(${CMAKE_INSTALL_PREFIX}/include)
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

add_executable( ${CLP}Test CompareFibers.cc )
target_link_libraries( ${CLP}Test ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${SphericalRidgelets_LIBRARY})
add_dependencies(${CLP}Test ${CLP})
set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

# remove output files before running tests
macro(RMTestFile CleanupTestName FileToRemove)
  add_test(NAME ${CleanupTestName}
    COMMAND ${CMAKE_COMMAND} -E remove ${FileToRemove}
    )
endmacro()

###############################################################################
# UKF Slicer CLI Test
###############################################################################
if (${UKFTractography_BUILD_SLICER_EXTENSION})
    if (MSVC) # need release subdir on Windows :|
      set(REL_VAR "/Release")
    else()
      set(REL_VAR "")
    endif()

  set(UKF_UNITTEST_LIB_PATHS "--additional-module-paths;${CMAKE_BINARY_DIR}/${Slicer_QTSCRIPTEDMODULES_LIB_DIR}${REL_VAR};${CMAKE_BINARY_DIR}/${Slicer_QTLOADABLEMODULES_LIB_DIR}${REL_VAR};${CMAKE_BINARY_DIR}/${Slicer_CLIMODULES_LIB_DIR}${REL_VAR}")

    set(UKF_TEST_CLI_PYTHON "${CMAKE_CURRENT_SOURCE_DIR}/slicer_cli_test.py")

    set(testname UKF_Slicer_CLI)
    add_test(
      NAME ${testname}
      COMMAND ${Slicer_LAUNCHER_EXECUTABLE}
      ${UKF_UNITTEST_LIB_PATHS}
      --python-code "data_dir = '${INPUT}'; exec(open('${UKF_TEST_CLI_PYTHON}').read(), globals())"
      )
    set_tests_properties(${testname} PROPERTIES TIMEOUT 180)

endif()
