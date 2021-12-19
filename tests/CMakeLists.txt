add_executable(test_sss "test_sss.cpp")
target_link_libraries(test_sss PRIVATE volppm)

add_executable(test_glass_sss "test_glass_sss.cpp")
target_link_libraries(test_glass_sss PRIVATE volppm)