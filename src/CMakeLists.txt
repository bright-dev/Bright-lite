INCLUDE_DIRECTORIES(${CYCLUS_CORE_INCLUDE_DIRS} ${CYCLUS_CORE_TEST_INCLUDE_DIR})

if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif()

# add lib
ADD_LIBRARY(brightlite burnupcalc.cpp origenBuilder.cpp)
INSTALL(TARGETS brightlite
    LIBRARY DESTINATION lib/cyclus/reactor_facility
    COMPONENT brightlite
    )

install_cyclus_standalone("ReactorFacility" "reactor_facility" "reactor")
TARGET_LINK_LIBRARIES(ReactorFacility brightlite)
