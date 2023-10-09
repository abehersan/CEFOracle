using CEFOracle
using DataFrames
using StatsPlots


ion = single_ion("Yb3")
bfactors = blm_dframe(Dict("B20"=>0.5622, "B40"=>1.6087e-5,
                           "B60"=>6.412e-7, "B66"=>-8.324e-6))

xyzw = SOPHE_xyzw(2, pi/2, 10, true)
ES = collect(-12:0.1:12)
BS = collect(0.0:0.5:12)
TS = collect(1.5:0.5:100)


function performance_overview()
    println("INS XSECTION POWDER CALCULATION, $(length(ES))-step energy scan")
    display(
        @benchmark cef_neutronxsection_powder!(ion, bfactors, DataFrame(T=1.5,
                    Q=2.55, EN=ES))
    )
    println()

    println("INS XSECTION SC CALCULATION, $(length(ES))-step energy scan")
    display(
    @benchmark cef_neutronxsection_crystal!(ion, bfactors,
                DataFrame(T=1.5, Qx=0.0, Qy=1.0, Qz=2.0, EN=ES,
                    Bx=0.0, By=0.0, Bz=0.0))
    )
    println()

    println("MAGNETIC MOMENT POWDER CALCULATION, $(size(xyzw)[1])-orientation powder average, $(length(BS))-magnetic field values")
    display(
        @benchmark cef_magneticmoment_powder!(ion, bfactors,
            DataFrame(T=1.5, B=BS), units="ATOMIC", xyzw=xyzw)
    )

    println("MAGNETIC MOMENT SC CALCULATION, $(length(BS))-magnetic field values")
    display(
        @benchmark cef_magneticmoment_crystal!(ion, bfactors,
            DataFrame(T=1.5, Bx=0.0, By=0.0, Bz=BS), units="ATOMIC")
    )

    println("SUSCEPTIBILITY POWDER CALCULATION, $(size(xyzw)[1])-orientation powder average, $(length(TS))-temperature values")
    display(
        @benchmark cef_susceptibility_powder!(ion, bfactors,
            DataFrame(B=0.01, T=TS), units="CGS", xyzw=xyzw)
    )

    println("SUSCEPTIBILITY SC CALCULATION, $(length(TS))-temperature values")
    display(
        @benchmark cef_susceptibility_crystal!(ion, bfactors,
            DataFrame(T=TS, Bx=0.0, By=0.0, Bz=0.01), units="CGS")
    )

    println("ENTROPY SC CALCULATION, $(length(TS))-temperature values")
    display(
        @benchmark cef_entropy!(ion, bfactors, DataFrame(T=TS), units="SI")
    )
end


performance_overview()


"""
results on my HP notebook - maximum calculation time expected per set is well
under a second and can be possibly lowered if the powder positions are set
to a minimum number that still allows the computation to converge


INS XSECTION POWDER CALCULATION, 241-step energy scan
BenchmarkTools.Trial: 60 samples with 1 evaluation.
 Range (min … max):  60.590 ms … 169.506 ms  ┊ GC (min … max):  9.84% … 4.97%
 Time  (median):     82.443 ms               ┊ GC (median):    10.25%
 Time  (mean ± σ):   84.814 ms ±  20.528 ms  ┊ GC (mean ± σ):   9.74% ± 1.98%

  █▂  ▂          ▂▂       ▂
  ██▅██▁▁██▁▅▁▅▁▅██▅▅▁███▅█▅▅▅▅█▁▁▁▁█▅▅▁▅▁▅▁▁▁▅▅▅▅▁▁▁▅▁▁▅▅▁▁▁▅ ▁
  60.6 ms         Histogram: frequency by time          125 ms <

 Memory estimate: 154.56 MiB, allocs estimate: 851802.

INS XSECTION SC CALCULATION, 241-step energy scan
BenchmarkTools.Trial: 62 samples with 1 evaluation.
 Range (min … max):  54.496 ms … 134.858 ms  ┊ GC (min … max): 12.18% … 7.56%
 Time  (median):     81.211 ms               ┊ GC (median):    10.19%
 Time  (mean ± σ):   82.354 ms ±  18.719 ms  ┊ GC (mean ± σ):   9.97% ± 2.34%

          █                   ▆
  ▄▁▁▁▄██▆█▄▄▆▁██▁▄▁▆▁▄▆█▄▁▄▆▄█▁▆▁▄▄▁▄▄▁▁▄▄▁█▁▁▁▁▁▁▄▁▁▁▁▁▁▁▁▄▄ ▁
  54.5 ms         Histogram: frequency by time          131 ms <

 Memory estimate: 154.64 MiB, allocs estimate: 852070.

MAGNETIC MOMENT POWDER CALCULATION, 100-orientation powder average, 25-magnetic field values
BenchmarkTools.Trial: 25 samples with 1 evaluation.
 Range (min … max):  158.602 ms … 346.883 ms  ┊ GC (min … max): 11.15% … 5.73%
 Time  (median):     192.961 ms               ┊ GC (median):    10.07%
 Time  (mean ± σ):   200.234 ms ±  41.010 ms  ┊ GC (mean ± σ):   9.97% ± 1.47%

  █ ▁     ▁      ▄
  █▁█▆▆▁▁▆█▁▆▆▁▆▆█▁▁▆▆▆▁▁▁▆▆▁▁▁▁▁▁▁▆▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▆ ▁
  159 ms           Histogram: frequency by time          347 ms <

 Memory estimate: 314.49 MiB, allocs estimate: 667649.
MAGNETIC MOMENT SC CALCULATION, 25-magnetic field values
BenchmarkTools.Trial: 2590 samples with 1 evaluation.
 Range (min … max):  1.001 ms … 10.896 ms  ┊ GC (min … max):  0.00% …  0.00%
 Time  (median):     1.293 ms              ┊ GC (median):     0.00%
 Time  (mean ± σ):   1.923 ms ±  1.210 ms  ┊ GC (mean ± σ):  10.95% ± 15.60%

  ▂▄█▄▂▁ ▁   ▁▁▁▁   ▁▂▁▁▂▁
  ███████████████▇▇███████▇▆▇██▇▇██▇▇▆▆▅▄▁▃▁▅▅▅▄▁▁▃▃▃▁▃▃▃▄▁▃ █
  1 ms         Histogram: log(frequency) by time     7.24 ms <

 Memory estimate: 3.27 MiB, allocs estimate: 7103.
SUSCEPTIBILITY POWDER CALCULATION, 100-orientation powder average, 198-temperature values
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.437 s …   1.498 s  ┊ GC (min … max): 9.58% … 9.02%
 Time  (median):     1.463 s              ┊ GC (median):    9.32%
 Time  (mean ± σ):   1.465 s ± 24.972 ms  ┊ GC (mean ± σ):  9.35% ± 0.30%

  █                    █    █                             █
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.44 s         Histogram: frequency by time         1.5 s <

 Memory estimate: 2.41 GiB, allocs estimate: 13757552.
SUSCEPTIBILITY SC CALCULATION, 198-temperature values
BenchmarkTools.Trial: 302 samples with 1 evaluation.
 Range (min … max):   7.839 ms … 49.962 ms  ┊ GC (min … max):  0.00% … 0.00%
 Time  (median):     15.588 ms              ┊ GC (median):    12.22%
 Time  (mean ± σ):   16.548 ms ±  6.428 ms  ┊ GC (mean ± σ):   8.24% ± 7.71%

   ▃▂▁▂▁ █▅ ▃ ▃ ▁  ▂ ▁▁  ▁       ▂
  ▄██████████▇█▄█▄▇█▇███▇█▅▄▇▄█▅▇█▃▇▄▅▆▄▆▅▄▃▄▃▃▁▁▁▁▃▃▁▁▃▃▄▃▁▃ ▄
  7.84 ms         Histogram: frequency by time        34.6 ms <

 Memory estimate: 24.72 MiB, allocs estimate: 137777.
ENTROPY SC CALCULATION, 198-temperature values
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  132.800 μs …   5.191 ms  ┊ GC (min … max): 0.00% … 85.71%
 Time  (median):     192.000 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   294.077 μs ± 317.641 μs  ┊ GC (mean ± σ):  9.44% ±  9.12%

  ▇█▅▅▄▄▄▄▅▅▃▁▁                                                 ▂
  █████████████████▇▇██▇▆▄▅▅▄▁▁▁▄▁▁▁▃▃▁▃▁▃▁▁▃▃▁▃▆▃▃▁▁▁▃▁▃▄▄▁▄▅▅ █
  133 μs        Histogram: log(frequency) by time       2.22 ms <

 Memory estimate: 525.81 KiB, allocs estimate: 3542.
"""