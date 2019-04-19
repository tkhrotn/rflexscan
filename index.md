rflexscan is an [R](https://www.r-project.org/) package for analysing spatial count data using the flexible spatial scan statistic developed by Tango and Takahashi (2005) and Kulldorff's circular spatial scan statistic (1997).
This package is designed for any of the following interrelated purposes:

1. To evaluate reported spatial disease clusters, to see if they are statistically significant.
2. To test whether a disease is randomly distributed over space.
3. To perform geographical surveillance of disease, to detect areas of significantly high rates.

rflexscan is developed based on the FleXScan 3.1.2 for Windows developed by Kunihiko Takahashi, Tetsuji Yokoyama and Toshiro Tango.


rflexscanは疾病集積性 (disease clustering) の検討をするための統計解析をR上で行うことができるパッケージです。
以下のような用途に使用できます。

1. 報告された疾患集積性を評価し、それらが統計的に有意であるかどうかを確かめる。
2. 疾患が空間的にランダムに分布しているかどうかを検定する。
3. 感染症のサーベイランス、著しく高い割合の地域を検出する。

rflexscanは高橋邦彦、横山徹爾、丹後俊郎によって開発されたWindows版FleXScan 3.1.2を元に開発しています。

# Install インストール
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
