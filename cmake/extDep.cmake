cmake_minimum_required(VERSION 3.10.3)

message("-- Fetching HLSLib")
# External Dependencies
## HLSLib - finds Intel OpenCL for FPGA installation
include(FetchContent)
FetchContent_Declare(
  hlslib
  GIT_REPOSITORY      https://github.com/definelicht/hlslib.git
  GIT_TAG             v1.4.2
)
FetchContent_MakeAvailable(hlslib)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${hlslib_SOURCE_DIR}/cmake)



