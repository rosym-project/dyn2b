if(CHECK_INCLUDE_DIR)
  # Already in cache, be silent
  set(CHECK_FIND_QUIETLY TRUE)
endif(CHECK_INCLUDE_DIR)

find_path(CHECK_INCLUDE_DIR NAMES check.h)

# Look for the library.
find_library(CHECK_LIBRARY NAMES check)

if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.4)
  # handle the QUIETLY and REQUIRED arguments and set CHECK_FOUND to TRUE if
  # all listed variables are TRUE
  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(Check "Please install 'check' and 'check-devel' packages" CHECK_LIBRARY CHECK_INCLUDE_DIR)
endif(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.4)

if(CHECK_FOUND)
  set(CHECK_LIBRARIES ${CHECK_LIBRARY})
else(CHECK_FOUND)
  set(CHECK_LIBRARIES)
endif(CHECK_FOUND)

mark_as_advanced(CHECK_INCLUDE_DIR)
mark_as_advanced(CHECK_LIBRARY)