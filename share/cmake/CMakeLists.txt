configure_file(triqstools-config.cmake.in triqstools-config.cmake @ONLY)
configure_file(triqstools-config-version.cmake.in triqstools-config-version.cmake @ONLY)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/triqstools-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/triqstools-config-version.cmake
  DESTINATION lib/cmake/triqstools
)

install(EXPORT triqstools-targets DESTINATION lib/cmake/triqstools)
