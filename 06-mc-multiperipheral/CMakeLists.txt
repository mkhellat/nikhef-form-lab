# Real CMake build for Topic 6 C gates (numeric + kinc2-driver).
# Portable across Linux / macOS / Windows. Extended-precision digit-loss
# tests self-skip when long double is not wider than double (MSVC).
#
#   cmake -S . -B build
#   cmake --build build
#   ctest --test-dir build --output-on-failure

cmake_minimum_required(VERSION 3.16)
project(nikhef_form_lab_mumu LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(CTest)
enable_testing()

set(MUMU_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(NUMERIC_DIR ${MUMU_ROOT}/numeric)
set(KINC2_DIR ${MUMU_ROOT}/kinc2)
set(DRIVER_DIR ${MUMU_ROOT}/kinc2-driver)

# libm is separate on Unix; folded into the CRT on MSVC.
set(MUMU_MATH_LIBS "")
if(NOT MSVC)
  list(APPEND MUMU_MATH_LIBS m)
endif()

# ---- numeric (header-only evaluators + gates) ----

add_executable(check_stable_c ${NUMERIC_DIR}/check_stable_c.c)
target_include_directories(check_stable_c PRIVATE ${NUMERIC_DIR})
target_link_libraries(check_stable_c PRIVATE ${MUMU_MATH_LIBS})
add_test(NAME check_stable_c COMMAND $<TARGET_FILE:check_stable_c>)

add_executable(digit_loss_mij ${NUMERIC_DIR}/digit_loss_mij.c)
target_include_directories(digit_loss_mij PRIVATE ${NUMERIC_DIR})
target_link_libraries(digit_loss_mij PRIVATE ${MUMU_MATH_LIBS})
add_test(NAME digit_loss_mij COMMAND $<TARGET_FILE:digit_loss_mij>)

add_executable(digit_loss_naive ${NUMERIC_DIR}/digit_loss_naive.c)
target_include_directories(digit_loss_naive PRIVATE ${NUMERIC_DIR})
target_link_libraries(digit_loss_naive PRIVATE ${MUMU_MATH_LIBS})

add_executable(mc_digit_loss ${NUMERIC_DIR}/mc_digit_loss.c)
target_include_directories(mc_digit_loss PRIVATE ${NUMERIC_DIR})
target_link_libraries(mc_digit_loss PRIVATE ${MUMU_MATH_LIBS})

# ---- vendored kinc2 (static lib) ----

add_library(kinc2 STATIC
  ${KINC2_DIR}/orient.c
  ${KINC2_DIR}/pickin.c
  ${KINC2_DIR}/mapt.c
  ${KINC2_DIR}/mapw.c
  ${KINC2_DIR}/mapla.c
  ${KINC2_DIR}/ipow.c
  ${KINC2_DIR}/iipow.c
  ${KINC2_DIR}/ranf.c
  ${KINC2_DIR}/gamgam.c
  ${KINC2_DIR}/peripp.c
  ${KINC2_DIR}/pi0.c
)
target_include_directories(kinc2 PUBLIC ${KINC2_DIR})
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  target_compile_options(kinc2 PRIVATE -Wno-unused-parameter)
endif()

# ---- kinc2-driver gates ----

add_executable(check_eemm_gamgam ${DRIVER_DIR}/check_eemm_gamgam.c)
target_include_directories(check_eemm_gamgam PRIVATE ${NUMERIC_DIR})
target_link_libraries(check_eemm_gamgam PRIVATE kinc2 ${MUMU_MATH_LIBS})
add_test(NAME check_eemm_gamgam COMMAND $<TARGET_FILE:check_eemm_gamgam>)

add_executable(digit_loss_eemm_gamgam ${DRIVER_DIR}/digit_loss_eemm_gamgam.c)
target_include_directories(digit_loss_eemm_gamgam PRIVATE ${NUMERIC_DIR})
target_link_libraries(digit_loss_eemm_gamgam PRIVATE kinc2 ${MUMU_MATH_LIBS})
add_test(NAME digit_loss_eemm_gamgam COMMAND $<TARGET_FILE:digit_loss_eemm_gamgam>)
