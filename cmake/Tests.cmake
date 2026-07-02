include(CTest)

set(TESTS

    testnumcpp-stats-dist-normal tests/stats/distributions/normal.cpp
    testnumcpp-stats-dist-gamma tests/stats/distributions/gamma.cpp
    testnumcpp-stats-dist-student tests/stats/distributions/tstudent.cpp
    testnumcpp-stats-covmat tests/stats/covmat.cpp
    testnumcpp-gquad-glaguerre tests/gquad/gausslaguerre.cpp
    testnumcpp-optim-nraphson tests/optim/nraphson.cpp
    testnumcpp-optim-bisection tests/optim/bisection.cpp
    testnumcpp-optim-secant tests/optim/secant.cpp
    testnumcpp-optim-brent tests/optim/brent.cpp
    testnumcpp-cfin-bs tests/cfin/bs.cpp
    testnumcpp-cfin-baw tests/cfin/baw.cpp
    testnumcpp-cfin-bachelier tests/cfin/bachelier.cpp
    #testnumcpp-cfin-bssim tests/cfin/bssim.cpp
    testnumcpp-cfin-heston tests/cfin/heston.cpp
    #testnumcpp-cfin-hestonqesim tests/cfin/hestonqesim.cpp
    testnumcpp-cfin-nss tests/cfin/nss.cpp
    testnumcpp-interpolation tests/interpolation/test.cpp
    testnumcpp-reg-ols tests/regression/ols.cpp
)

while(TESTS)
    list(POP_FRONT TESTS test_name test_file)

    add_executable(${test_name} ${test_file})
    target_link_libraries(${test_name} PRIVATE numcpp)

    add_test(NAME ${test_name} COMMAND ${test_name})
endwhile()