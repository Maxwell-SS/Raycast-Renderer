cmake_minimum_required(VERSION 3.2)
set(CMAKE_CXX_STANDARD 20)
project(RaycastRenderer)

# delete folders
execute_process(COMMAND ${CMAKE_COMMAND} -E rm -rf res)

# copying shaders and resources file into build folder
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/res DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

include_directories(
	${PROJECT_SOURCE_DIR}/src
	${PROJECT_SOURCE_DIR}/libs/SDL2/build/include/SDL2
)

add_subdirectory(${PROJECT_SOURCE_DIR}/libs/SDL2)

# making main executable
add_executable(${PROJECT_NAME} src/main.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE SDL2-static)