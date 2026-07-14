add_library(numcpp STATIC
    src/cfin/letsberational.cpp
    src/cfin/baw.cpp
    src/stats/univariate.cpp
)
target_link_libraries(numcpp PUBLIC Eigen3::Eigen)
target_include_directories(numcpp PUBLIC include)
