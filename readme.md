Multivariate mixed membership modeling: Inferring domain-specific risk profiles
===============================================================================

This tutorial shows how to reproduce simulation results in Section 6
[**Multivariate mixed membership modeling: Inferring domain-specific
risk profiles**](https://arxiv.org/pdf/1901.05191.pdf).

As first step, download and install the package <tt>MMM</tt> by running
the following code (or using <tt>devtools</tt>)

``` r
R CMD build MMM
R CMD INSTALL MMM -l your_path
```

Once the R package has been correctly installed, the following
<tt>.md</tt> files provide instruction and R code to reproduce Figure
2–5 and Tables 1 and 2 of the paper, and figure S2–S5 and Tables S1 and
S2, in the Supplementary Material.

-   The file [`simulated_data.md`](simulated_data.md) contains the R
    code to generate the data for all the proposed scenarios.

-   The file [`simulation_MMM.md`](simulation_MMM.md) contains
    instruction on how to estimate the proposed MMM model on the
    simulated data.

-   The file [`simulation_MM.md`](simulation_MM.md) contains instruction
    on how to estimate MM models based on R package <tt>mixedMem</tt>.

-   The file [`simulation_MM_NIMBLE.md`](simulation_MM_NIMBLE.md)
    contains instruction on how to estimate MM models using the software
    <tt>NIMBLE</tt>.

-   The file [`plots_and_tables.md`](plots_and_tables.md) contains code
    to produce Figures 3 to 5, and Tables 1 and 2.

-   The file [`plots_and_tables_NIMBLE.md`](plots_and_tables_NIMBLE.md)
    contains code to produce Figures S1 to S5, and Tables S1 and S2 of
    the Supplementary Material.
