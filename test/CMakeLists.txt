add_executable(main_test
  main_test.c
  linear_algebra_test.c
  geometry_test.c
  kinematic_chain_test.c
)

target_link_libraries(main_test
  dyn2b
  ${CHECK_LIBRARY}
)

add_test(main_test
  ${CMAKE_CURRENT_BINARY_DIR}/main_test
)