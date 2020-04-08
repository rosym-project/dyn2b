add_library(dyn2b
  dyn2b/linear_algebra.c
  dyn2b/geometry.c
)

target_link_libraries(dyn2b m)