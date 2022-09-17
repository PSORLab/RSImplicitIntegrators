using Plots, DynamicBoundsBase, DynamicBoundspODEsDiscrete, LaTeXStrings, Measures
using DataFrames, CSV, DifferentialEquations

const DBpODEsDiscrete = DynamicBoundspODEsDiscrete
const DBB = DynamicBoundsBase
const DBDisc = DynamicBoundspODEsDiscrete
LINEWIDTH = 3
use_test_factory = false