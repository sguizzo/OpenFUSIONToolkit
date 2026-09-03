oft_add_exe( hartmann.F90 )
install(FILES oft_hartmann.in DESTINATION examples/MUG/hartmann COMPONENT examples)
configure_file( oft_hartmann.in ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
