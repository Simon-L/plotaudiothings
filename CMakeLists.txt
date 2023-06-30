cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)

project(
    plotenv
    VERSION 0.1.0
    DESCRIPTION "Short description"
    HOMEPAGE_URL "https://example.com/"
    LANGUAGES CXX
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare library ----

add_library(
    plotenv_lib OBJECT
    source/lib.cpp
)

target_include_directories(
    plotenv_lib ${warning_guard}
    PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>"
)

target_compile_features(plotenv_lib PUBLIC cxx_std_17)

# ---- Declare executable ----

add_subdirectory(sst-basic-blocks)

add_executable(plotenv_exe source/main.cpp)
add_executable(plotenv::exe ALIAS plotenv_exe)

set_property(TARGET plotenv_exe PROPERTY OUTPUT_NAME plotenv)

target_compile_features(plotenv_exe PRIVATE cxx_std_17)

target_link_libraries(plotenv_exe PRIVATE plotenv_lib)
target_link_libraries(plotenv_exe PRIVATE sst-basic-blocks)
target_include_directories(plotenv_exe PUBLIC dr_libs)
target_include_directories(plotenv_exe PUBLIC include)

add_custom_target(run
    DEPENDS plotenv::exe
    COMMAND plotenv::exe
    COMMAND ../audiowaveform -i out.wav -o out.png -w 800 -h 700 -e 0.3)
    # COMMAND xdg-open out.png)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Developer mode ----

if(NOT plotenv_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of plotenv"
  )
endif()

include(cmake/dev-mode.cmake)
