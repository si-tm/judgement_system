################################################################################

set(link_outputs)

foreach(x ${NUPACK_PYTHON_FILES})
    link_file(${CMAKE_CURRENT_SOURCE_DIR}/python/${x}
              ${CMAKE_CURRENT_BINARY_DIR}/nupack/${x})
endforeach()

foreach(x dna04.json rna95.json rna06.json rna99.json dna04-nupack3.json rna95-nupack3.json rna99-nupack3.json)
    link_file(${CMAKE_CURRENT_SOURCE_DIR}/parameters/${x}
              ${CMAKE_CURRENT_BINARY_DIR}/nupack/parameters/${x})
endforeach()

foreach(x ${NUPACK_PYTHON_TEST_FILES})
    link_file(${CMAKE_CURRENT_SOURCE_DIR}/test/python/${x}
              ${CMAKE_CURRENT_BINARY_DIR}/nupack/test/${x})
endforeach()

foreach(x ${REBIND_PYTHON_FILES})
    link_file(${CMAKE_CURRENT_SOURCE_DIR}/external/rebind/${REBIND_PYTHON_ROOT}/${x}
              ${CMAKE_CURRENT_BINARY_DIR}/nupack/${x})
endforeach()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/package/setup.py
               ${CMAKE_CURRENT_BINARY_DIR}/setup.py)

link_file(${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt
          ${CMAKE_CURRENT_BINARY_DIR}/LICENSE.txt)
# file(WRITE ${BUILD_DIR}/build.sh "echo 'Installing nupack module'\n$PYTHON -m pip install --no-deps --ignore-installed .\n\n")

add_custom_target(copy-python ALL DEPENDS ${link_outputs})
set(NUPACK_PYTHON_OUTPUTS ${link_outputs} CACHE STRING "list of python output files" FORCE)

################################################################################

if(NUPACK_BUILD_CXX)
    rebind_module(nupack-python cpp nupack-bind)
else()
    set(NUPACK_BIND_LIBRARY)
    if(NUPACK_BIND_LIBRARY)
        set(NUPACK_BIND_PATH ${NUPACK_BIND_LIBRARY})
        message("-- Using binding library ${NUPACK_BIND_PATH}")
    else()
        find_library(NUPACK_BIND_PATH NAMES nupack-bind)
        if(NUPACK_BIND_PATH)
            message("-- Found binding library ${NUPACK_BIND_PATH}")
        else()
            message(FATAL_ERROR "libnupack-bind cannot be found")
        endif()
    endif()
    add_library(nupack-bind-import INTERFACE)
    target_link_libraries(nupack-bind-import INTERFACE ${NUPACK_BIND_PATH})
    rebind_module(nupack-python cpp nupack-bind-import)
endif()

set_target_properties(nupack-python PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/nupack) # CXX_VISIBILITY_PRESET hidden)
add_dependencies(nupack-python copy-python)
add_library(nupack::python ALIAS nupack-python)
nupack_cxx_target(nupack-python)

################################################################################
