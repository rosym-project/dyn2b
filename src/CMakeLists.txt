add_library(dyn2b SHARED
  dyn2b/linear_algebra.c
  dyn2b/geometry.c
  dyn2b/mechanics.c
  dyn2b/kinematic_chain.c

  dyn2b/geometry_nbx.c
  dyn2b/kinematic_chain_nbx.c
)

target_link_libraries(dyn2b m)


add_library(dyn2b_example SHARED
  example/chain_iterator.c
  example/solver_state.c
  example/robots/one_dof.c
  example/robots/two_dof.c
)

target_link_libraries(dyn2b_example dyn2b)


add_executable(fpk_algorithm example/fpk_algorithm.c)
target_link_libraries(fpk_algorithm dyn2b_example)

add_executable(fpk_algorithm_iter example/fpk_algorithm_iter.c)
target_link_libraries(fpk_algorithm_iter dyn2b_example)

add_executable(fvk_algorithm example/fvk_algorithm.c)
target_link_libraries(fvk_algorithm dyn2b_example)

add_executable(fak_algorithm example/fak_algorithm.c)
target_link_libraries(fak_algorithm dyn2b_example)

add_executable(ab_algorithm example/ab_algorithm.c)
target_link_libraries(ab_algorithm dyn2b_example)


install(
  TARGETS dyn2b
  EXPORT dyn2b-targets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install(
  TARGETS dyn2b_example
  EXPORT dyn2b_example-targets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
