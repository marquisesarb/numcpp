include(CTest)

set(TESTS

    testnumcpp-stats-dist-normal tests/stats/distributions/normal.cpp
    testnumcpp-stats-dist-gamma tests/stats/distributions/gamma.cpp
    testnumcpp-stats-dist-student tests/stats/distributions/tstudent.cpp
    testnumcpp-stats-covmat tests/stats/covmat.cpp
    testnumcpp-stats-univariate tests/stats/univariate.cpp
    testnumcpp-gquad-glaguerre tests/gquad/gausslaguerre.cpp
    testnumcpp-optim-nraphson tests/optim/nraphson.cpp
    testnumcpp-optim-bisection tests/optim/bisection.cpp
    testnumcpp-optim-secant tests/optim/secant.cpp
    testnumcpp-optim-brent tests/optim/brent.cpp
    testnumcpp-reg-ols tests/regression/ols.cpp
    testnumcpp-interpolation tests/interpolation/test.cpp
    testnumcpp-solvers-tridiagthomas tests/solvers/tridiagthomas.cpp

    testnumcpp-cfin-nss tests/cfin/nss.cpp
    testnumcpp-cfin-ssvi tests/cfin/ssvi.cpp

    testnumcpp-cfin-blackscholes-cf tests/cfin/blackscholes/cf.cpp
    testnumcpp-cfin-blackscholes-mcbs tests/cfin/blackscholes/mcbs.cpp
    testnumcpp-cfin-blackscholes-mclv tests/cfin/blackscholes/mclv.cpp
    testnumcpp-cfin-blackscholes-baw tests/cfin/blackscholes/baw.cpp
    testnumcpp-cfin-blackscholes-pdesolver-vanilla tests/cfin/blackscholes/pdesolver/vanilla.cpp

    testnumcpp-cfin-bachelier-cf tests/cfin/bachelier/cf.cpp

    testnumcpp-cfin-heston-cf tests/cfin/heston/cf.cpp
    testnumcpp-cfin-heston-qemc tests/cfin/heston/qemc.cpp

)

while(TESTS)
    list(POP_FRONT TESTS test_name test_file)

    add_executable(${test_name} ${test_file})
    target_link_libraries(${test_name} PRIVATE numcpp)

    add_test(NAME ${test_name} COMMAND ${test_name})
endwhile()