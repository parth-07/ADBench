project("ADBench_CladLatest" CXX)

add_library("CladLatest" MODULE)

SET(clad_install_root "/home/mad-scientist/clubs/open-source/clad/formal/builds/build-llvm-16/inst")

target_sources("CladLatest" PRIVATE "${CMAKE_SOURCE_DIR}/src/cpp/shared/utils.cpp")
target_sources("CladLatest" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/ba/ba_grad.cpp")
target_sources("CladLatest" PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/CladLatestBA.cpp")
target_sources("CladLatest" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/gmm/gmm_grad.cpp")
target_sources("CladLatest" PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/CladLatestGMM.cpp")

target_include_directories("CladLatest" PRIVATE ${clad_install_root}/include)
target_include_directories("CladLatest" PRIVATE ${CMAKE_SOURCE_DIR}/src/cpp/shared)
# add_compile_definitions("CladLatest" CLAD_NO_NUM_DIFF)