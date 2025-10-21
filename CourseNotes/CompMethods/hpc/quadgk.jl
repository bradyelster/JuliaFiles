using QuadGK
using SpecialFunctions

f(x) = exp(-x^4)

quadgk(f, 0, Inf)[1]
(1 / 4) * gamma(1 / 4)

