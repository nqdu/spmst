include(GoogleTest)

add_executable(unit_tests
    test_bilinear.cpp
    test_gps2dist.cpp
    test_parallel_tools.cpp
    test_io.cpp
    test_invparam.cpp
)
target_link_libraries(unit_tests share GTest::gtest_main)

gtest_discover_tests(unit_tests)
