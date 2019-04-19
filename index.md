The rflexscan package provides functions and classes to analyze spatial count data using the flexible spatial scan statistic developed by Tango and Takahashi (2005) on [R](https://www.r-project.org/).

rflexscanは疾病集積性（disease clustering）の検討をするための統計解析をR上で行うことができるパッケージです。集積性の検定には、Kulldorff’s Circular Scan法、Tango and Takahashi’s Flexible Scan法の2つによって検討できます。
このパッケージはWindows版FleXScanを元に開発しています。

# Install
```r
# Install development version from GitHub:
install.packages("devtools")
devtools::install_github("tkhrotn/rflexscan")

# Install from CRAN (in preparation!)
# install.packages("rflexscan")
```


# References
 * Tango T. and Takahashi K. (2005). A flexibly shaped spatial scan statistic for detecting clusters, International Journal of Health Geographics 4:11.
 * Takahashi K, Yokoyama T and Tango T. (2010). FleXScan v3.1: Software for the Flexible Scan Statistic. National Institute of Public Health, Japan.
   <https://sites.google.com/site/flexscansoftware/home>
