# header-only library
vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO p12tic/libsimdpp
    REF 9dac213d6965bf57c7accf9dda6d16c3bed0e3ac
    SHA512 d90af5d34bd8f52e1a5917e17ad2b77d8cb2a50a93af45f2dc361ed4df276c9ac6d21f0e646a8c38235227eb10544388e5f3b476a8dface6acebb440b665c1fd
    HEAD_REF master
)

vcpkg_configure_cmake(
    SOURCE_PATH ${SOURCE_PATH}
    PREFER_NINJA
)

vcpkg_install_cmake()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")

# Include architecture detection
configure_file("${SOURCE_PATH}/cmake/SimdppMultiarch.cmake" "${SOURCE_PATH}/libsimdpp-config.cmake" COPYONLY)
file(APPEND "${SOURCE_PATH}/libsimdpp-config.cmake" "\nset(LIBSIMDPP_INCLUDE_DIRS \"${CURRENT_PACKAGES_DIR}/include/libsimdpp-2.1\")\n")
file(INSTALL "${SOURCE_PATH}/libsimdpp-config.cmake" DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}")

# Handle copyright
file(INSTALL ${SOURCE_PATH}/LICENSE_1_0.txt DESTINATION ${CURRENT_PACKAGES_DIR}/share/${PORT} RENAME copyright)

