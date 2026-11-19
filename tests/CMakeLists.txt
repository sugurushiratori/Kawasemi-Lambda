add_executable(test_config test_config.cpp)
target_link_libraries(test_config PRIVATE lambda_core)
add_test(NAME config COMMAND test_config)

add_executable(test_solver test_solver.cpp)
target_link_libraries(test_solver PRIVATE lambda_core)
target_compile_definitions(test_solver PRIVATE
  LAMBDA_REFERENCE_CSV="${CMAKE_CURRENT_SOURCE_DIR}/reference/chandrasekhar_roots.csv")
add_test(NAME chandrasekhar COMMAND test_solver)

add_executable(test_convergence test_convergence.cpp)
target_link_libraries(test_convergence PRIVATE lambda_core)
add_test(NAME convergence COMMAND test_convergence)

# end-to-end: the example inputs must run
file(GLOB _examples ${PROJECT_SOURCE_DIR}/examples/*.yaml)
foreach(_ex ${_examples})
  get_filename_component(_name ${_ex} NAME_WE)
  add_test(NAME example_${_name} COMMAND lambdaDrop ${_ex} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()
add_test(NAME cli_legacy COMMAND lambdaDrop 1e-3 2 0 32 32)

add_executable(test_eigenfunction test_eigenfunction.cpp)
target_link_libraries(test_eigenfunction PRIVATE lambda_core)
add_test(NAME eigenfunction COMMAND test_eigenfunction)
