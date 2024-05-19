project("ADBench_CladLatestTBR" CXX)

add_library("CladLatestTBR" MODULE)

SET(clad_install_root "/home/mad-scientist/clubs/open-source/clad/formal/builds/build-llvm-16/inst")

target_sources("CladLatestTBR" PRIVATE "${CMAKE_SOURCE_DIR}/src/cpp/shared/utils.cpp")
target_sources("CladLatestTBR" PRIVATE "ba/ba_grad.cpp")
target_sources("CladLatestTBR" PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/CladLatestTBRBA.cpp")
target_sources("CladLatestTBR" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/gmm/gmm_grad_tbr.cpp")
target_sources("CladLatestTBR" PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/CladLatestTBRGMM.cpp")

target_include_directories("CladLatestTBR" PRIVATE ${clad_install_root}/include)
target_include_directories("CladLatestTBR" PRIVATE ${CMAKE_SOURCE_DIR}/src/cpp/shared)
# add_compile_definitions("CladLatestTBR" CLAD_NO_NUM_DIFF)