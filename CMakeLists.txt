cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)

project(
    plotaudiothings
    VERSION 0.1.0
    DESCRIPTION "Plot envlopes and filter response!"
    HOMEPAGE_URL ""
    LANGUAGES CXX
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare executable ----

add_executable(plotaudiothings_exe source/main.cpp)
add_executable(plotaudiothings::exe ALIAS plotaudiothings_exe)

set_property(TARGET plotaudiothings_exe PROPERTY OUTPUT_NAME plotaudiothings)

target_compile_features(plotaudiothings_exe PRIVATE cxx_std_17)

add_subdirectory(JUCE)
add_subdirectory(chowdsp_utils)

target_link_libraries(plotaudiothings_exe PUBLIC
    juce::juce_audio_utils
    juce::juce_dsp
    juce::juce_gui_basics
    juce::juce_gui_extra
    chowdsp::chowdsp_dsp_utils
    chowdsp::chowdsp_sources
    chowdsp::chowdsp_filters
    chowdsp::chowdsp_visualizers
    chowdsp::chowdsp_plugin_utils
)

add_subdirectory(sst-basic-blocks)

target_link_libraries(plotaudiothings_exe PRIVATE sst-basic-blocks)
target_link_libraries(plotaudiothings_exe PRIVATE curl)
target_include_directories(plotaudiothings_exe PUBLIC chowdsp_wdf/include)
target_include_directories(plotaudiothings_exe PUBLIC dr_libs)
target_include_directories(plotaudiothings_exe PUBLIC /usr/include/gtk-3.0
    /usr/include/glib-2.0
    /usr/include/harfbuzz
    /usr/include/cairo
    /usr/include/gdk-pixbuf-2.0
    /usr/include/atk-1.0
    /usr/include/webkitgtk-4.0
    /usr/include/libsoup-2.4
    /usr/lib/x86_64-linux-gnu/glib-2.0/include
    /usr/include/pango-1.0)
target_include_directories(plotaudiothings_exe PUBLIC sciplot)
target_include_directories(plotaudiothings_exe PUBLIC include)

add_custom_target(run
    DEPENDS plotaudiothings::exe
COMMAND plotaudiothings::exe ../options.json)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Developer mode ----

if(NOT plotaudiothings_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of plotaudiothings"
  )
endif()

include(cmake/dev-mode.cmake)
