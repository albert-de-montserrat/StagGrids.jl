using StagGrids, Test

function run_tests()
    f = filter(x -> contains(x, "test_"), readdir(pwd()))
    map(include, f)
end

run_tests()

