# Declare your game executable target
add_executable(game)
set_target_properties(game PROPERTIES OUTPUT_NAME "${CMAKE_PROJECT_NAME}")




# set the target source files
file(GLOB_RECURSE MY_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
target_sources(game PRIVATE ${MY_SOURCES})


# Link the engine library to the game executable
target_link_libraries(game PRIVATE LittleEngine)

# Include directories if needed (usually engine's include is enough via target_link_libraries)
target_include_directories(game PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(game PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/gameLayer/")

