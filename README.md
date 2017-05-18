
# ApproXD.jl: Function Approximation in X Dimensions in Julia

[![Build Status](https://travis-ci.org/michakraus/ApproXD.jl.png?branch=master)](https://travis-ci.org/michakraus/ApproXD.jl)
[![Coverage Status](https://coveralls.io/repos/github/michakraus/ApproXD.jl/badge.svg?branch=master)](https://coveralls.io/github/michakraus/ApproXD.jl?branch=master)


[![2D approximation](https://dl.dropboxusercontent.com/u/109115/BSplines.jl/approx.png)]()

This package implements high-dimensional approximation of real-valued functions in Julia:

```julia
x = [1:10.0]
y = [1:10.0]
z = [1:10.0]
fun(x,y,z) = x^2 + y^(0.5) + (x-z)^2

# what is fun(4.5,pi,8.01) ?
```

[Documentation is available on readthedocs](http://approxdjl.readthedocs.org/en/latest/index.html).
