install(
    TARGETS plotenv_exe
    RUNTIME COMPONENT plotenv_Runtime
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()
