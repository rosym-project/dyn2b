add_library(dyn2b
  dyn2b/linear_algebra.c
  dyn2b/geometry.c
  dyn2b/kinematic_chain.c
)

target_link_libraries(dyn2b m)


add_library(dyn2b_example
  example/solver_state.c
  example/robots/one_dof.c
  example/robots/two_dof.c
)

target_link_libraries(dyn2b_example dyn2b)


add_executable(fpk_algorithm example/fpk_algorithm.c)
target_link_libraries(fpk_algorithm dyn2b_example)


install(
  TARGETS dyn2b
  EXPORT dyn2b-targets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
