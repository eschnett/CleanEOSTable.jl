using CleanEOSTable
using Base.Test
using HDF5

mktempdir() do path
    println("Downloading...")
    download("https://stellarcollapse.org/EOS/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2")
    println("Unpacking...")
    run(`pwd`)
    run(`ls`)
    run(`ls $path`)
    run(`bunzip2 $path/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2`)

    tablename = joinpath(path, "$path/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5")
    new_tablename = joinpath(path, "$path/LS220_234r_136t_50y_analmu_20091212_SVNr26_corrected.h5")
    cleanEOSTable(tablename, new_tablename)

    println("Comparing tables...")
    table = h5open(tablename, "r")
    table_logenergy = read(table, "logenergy")
    new_table = h5open(new_tablename, "r")
    new_table_logenergy = read(new_table, "logenergy")

    @test count(new_table_logenergy .!= table_logenergy) == 741
end
