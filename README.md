![license](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 
[![Travis build status](https://travis-ci.org/Jean-Romain/lidRplugins.svg?branch=master)](https://travis-ci.com/Jean-Romain/lidRplugins)
[![Codecov test coverage](https://codecov.io/gh/Jean-Romain/lidRplugins/branch/master/graph/badge.svg)](https://codecov.io/gh/Jean-Romain/lidRplugins?branch=master)

This package contains functions and algorithms to extend the [lidR](https://github.com/Jean-Romain/lidR) package (versions >= 3.1). These functions or algorithms are not yet or will not be included in the `lidR` package either because they are:

- :microscope: **Experimental** and not supported by a peer-reviewed and accessible publications.
- :zap: **Non suitable for `lidR`**  usually because they are not sufficiently efficient. 
- :warning: **Not tested enought** and I'm not sure they are sufficiently robust.
- :octocat: **Require extra packages** available on github but not on CRAN

This package will NOT be submitted on CRAN and must be installed from github. It depends on `lidR (>= 3.1.0)` and should be seen as a laboratory with more or less interesting content inside.   

 
```r
remotes::install_github("Jean-Romain/lidRplugins")
```

To install the package from github make sure you have a working development environment.

* **Windows**: Install [Rtools.exe](https://cran.r-project.org/bin/windows/Rtools/).  
* **Mac**: Install `Xcode` from the Mac App Store.
* **Linux**: Install the R development package, usually called `r-devel` or `r-base-dev`
