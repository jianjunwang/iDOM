
# iDOM

Statistical analysis of dissolved organic matter characterized by
high-resolution mass spectrometry

## About

iDOM is a multifunctional tool that facilitates basic analyses, such as
the calculation of molecular traits, the assignment of molecular
classes, and the evaluation of chemical diversity and dissimilarity. It
also includes functions for advanced analyses to quantify the assembly
processes of DOM assemblages, the effect of molecular dark matter on DOM
molecular interactions, and the associations between DOM molecules and
microbial taxa.

We expect that iDOM will serve as a comprehensive pipeline for DOM
statistical analyses and bridge the gap between chemical
characterization and ecological interpretation.

Package maintainer: Jianjun Wang (<jjwang@niglas.ac.cn>)

Developers: Ang Hu, Fanfan Meng

## Installation

You can install the development version of **iDOM** from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jianjunwang/iDOM")
```

Installing packages in R from zip source. You may have downloaded
**iDOM** in ***zip*** or ***tar.gz*** format. In order to install the
package from a local zip file, you just need to call the
install.packages function with arguments repos = NULL and type =
“source”.

``` r
install.packages("file_path/package_file_name.extension",repos = NULL, type = "source")
```

## Key functions in iDOM package

[commProc](https://pubs.acs.org/doi/abs/10.1021/acs.est.2c01432):
Quantifies the relative influences of deterministic and stochastic
processes in the assembly of DOM assemblages.

[iDME](https://pubs.acs.org/doi/10.1021/acs.est.2c05052): Quantifies the
effect of molecular dark matter on DOM assemblages by constructing
co-occurrence networks based on the presence and absence of dark matter.

[H2](https://doi.org/10.1101/2021.08.12.456177): Calculates the
network-level specialization of all interacting trophic levels in
DOM-microbe bipartite networks, including full, negative and positive
networks.

[iCER](https://www.nature.com/articles/s41467-024-44813-2): Calculates
the indicator of compositional-level environmental response to assess
the thermal response of DOM.

## References

A Hu, K-S Jang, A J Tanentzap, W Zhao, J T Lennon, J Liu, M Li, J
Stegen, M Choi, Y Lu, X Feng, and J Wang. 2024. Thermal responses of
dissolved organic matter under global change. ***Nature
communications***. <https://www.nature.com/articles/s41467-024-44813-2>

A Hu, M Choi, A J Tanentzap, J Liu, K-S Jang, J T Lennon, Y Liu, J
Soininen, X Lu, Y Zhang, J Shen, and J Wang. 2022. Ecological networks
of dissolved organic matter and microorganisms under global change.
***Nature communications***.
<https://www.nature.com/articles/s41467-022-31251-1>

A Hu, F Meng, A J Tanentzap, K-S Jang, and J Wang. 2023. Dark Matter
Enhances Interactions within Both Microbes and Dissolved Organic Matter
under Global Change. ***Environmental Science & Technology***.
<https://pubs.acs.org/doi/10.1021/acs.est.2c05052>

A Hu, K-S Jang, F Meng, J Stegen, A J Tanentzap, M Choi, J T Lennon, J
Soininen, and J Wang. 2022. Microbial and Environmental Processes Shape
the Link between Organic Matter Functional Traits and Composition.
***Environmental Science & Technology***.
<https://pubs.acs.org/doi/abs/10.1021/acs.est.2c01432>

F Meng, A Hu, K-S Jang, and J Wang. 2024. iDOM: Statistical analysis of
dissolved organic matter based on high-resolution mass spectrometry.
***Research Square***. <https://doi.org/10.21203/rs.3.rs-4660944/v1>

F Meng, A Hu, K-S Jang, and J Wang. 2025. iDOM: Statistical analysis of 
dissolved organic matter characterized by high‐resolution mass spectrometry. 
***mLife***. <https://doi.org/10.1002/mlf2.70002>.

