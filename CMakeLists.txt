cmake_minimum_required(VERSION 3.16)
project(TemplateGame)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



# ----------------- TODO: move this to configuration in CMakeSettings.json
# Build options
set(PRODUCTION_BUILD OFF CACHE BOOL "Make this a production build" FORCE)
set(PLATFORM "GLFW" CACHE STRING "Choose platform: GLFW or SDL" FORCE)
# -----------------

# Add engine as a library
add_subdirectory(LittleEngine)

# Add game as an executable that links to engine
add_subdirectory(game)






if(PRODUCTION_BUILD)
	# setup the ASSETS_PATH macro to be in the root folder of your exe
	target_compile_definitions(game PUBLIC RESOURCES_PATH="./resources/") 
	# remove the option to debug asserts.
	target_compile_definitions(game PUBLIC PRODUCTION_BUILD) 
	target_compile_definitions(LittleEngine PUBLIC PRODUCTION_BUILD)

	
	#may give problems on linux
	set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

else()
	# This is useful to get an ASSETS_PATH in your IDE during development
	target_compile_definitions(game PUBLIC RESOURCES_PATH="${CMAKE_CURRENT_SOURCE_DIR}/game/resources/")
	target_compile_definitions(game PUBLIC DEVELOPMENT_BUILD) 
	target_compile_definitions(LittleEngine PUBLIC DEVELOPMENT_BUILD) 
 
	# enable IMGUI
	target_compile_definitions(game PUBLIC ENABLE_IMGUI)
	target_compile_definitions(LittleEngine PUBLIC ENABLE_IMGUI)

	
	set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)

	# does not work with nija i think
	## For MSVC compiler only
	#if (MSVC)
	#	# Add these flags only for Debug configuration
	#	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /ZI /Od")
	#endif()

endif()



if(MSVC) # If using the VS compiler...

	target_compile_definitions(game PUBLIC _CRT_SECURE_NO_WARNINGS)

	
	set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>$<$<CONFIG:Release>:Release>")


	add_compile_options(/arch:AVX2) #make sure SIMD optimizations take place
	
if(PRODUCTION_BUILD)
	set_target_properties(game PROPERTIES LINK_FLAGS "/SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup") #no console
	endif()
endif()

#force remove unicode
if (WIN32)
	target_compile_options(game PRIVATE /UUNICODE /U_UNICODE)
endif()


