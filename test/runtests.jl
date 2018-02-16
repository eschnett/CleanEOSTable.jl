using CleanEOSTable
using Base.Test
using HDF5

mktempdir() do path
    println("Downloading...")
    download("https://stellarcollapse.org/EOS/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2")
    println("Unpacking...")
    run(`bunzip2 LS220_234r_136t_50y_analmu_20091212_SVNr26.h5.bz2`)

    tablename = joinpath(path, "LS220_234r_136t_50y_analmu_20091212_SVNr26.h5")
    new_tablename = joinpath(path, "LS220_234r_136t_50y_analmu_20091212_SVNr26_corrected.h5")
    cleanEOSTable(tablename, new_tablename)

    @test count(new_table_logenergy .!= table_logenergy) == 741
end
