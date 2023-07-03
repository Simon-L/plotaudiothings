install(
    TARGETS plotaudiothings_exe
    RUNTIME COMPONENT plotaudiothings_Runtime
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()
