project("ADBench_Clad" CXX)

add_library("Clad" MODULE)

SET(clad_install_root "/home/parth/clubs/open-source/clad/old/builds/build-15/inst")

target_sources("Clad" PRIVATE "${CMAKE_SOURCE_DIR}/src/cpp/shared/utils.cpp")
target_sources("Clad" PRIVATE "ba/ba.cpp")
target_sources("Clad" PRIVATE "ba/ba_grad.cpp")
target_sources("Clad" PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/CladBA.cpp")
target_include_directories("Clad" PRIVATE ${clad_install_root}/include)
add_compile_definitions("Clad" CLAD_NO_NUM_DIFF)