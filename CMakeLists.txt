INCLUDE_DIRECTORIES(${CYCLUS_CORE_INCLUDE_DIRS} ${CYCLUS_CORE_TEST_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif()

# add lib
ADD_LIBRARY(brightlite burnupcalc.cpp origenBuilder.cpp)
INSTALL(TARGETS brightlite
    LIBRARY DESTINATION lib/cyclus
    COMPONENT brightlite
    )
use_cyclus("Brightlite" "reactor_facility")
use_cyclus("Brightlite" "fuelfab_facility")
use_cyclus("Brightlite" "reprocess_facility")
install_cyclus_module("Brightlite" "")
TARGET_LINK_LIBRARIES(Brightlite brightlite)
TARGET_LINK_LIBRARIES(Brightlite_unit_tests brightlite)
#SET_TARGET_PROPERTIES(Brightlite brightlite PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)

