# Analysis of Circumcision Data
Sean Pascoe
2024-03-27

- [Overview + Highlights](#overview--highlights)
  - [Bigger Plots](#bigger-plots)
  - [Specific Insights](#specific-insights)
  - [Remaining Questions/Concerns](#remaining-questionsconcerns)
  - [Conclusions + Future Analysis](#conclusions--future-analysis)
- [In Depth Analysis](#in-depth-analysis)
  - [Virion Particle count](#virion-particle-count)
  - [Mast Cell Count Analysis](#mast-cell-count-analysis)
  - [Immunofluorescence Cell
    Distribution](#immunofluorescence-cell-distribution)
    - [Adding Back virion penetration
      data](#adding-back-virion-penetration-data)

<details class="code-fold">
<summary>Code</summary>

``` r
library(tidyverse)
library(DT)
```

</details>

# Overview + Highlights

This analysis analyzes \[\]

In this sort of analysis, our expectation would be that shaft samples
look the same between circumcised and uncircumcised individuals, and
that the glans samples react to virus differently and/or have different
resident immune cell populations.

- There were x,y,z differences that we found

- These differences were/were not the same in shaft/glans

- We did not detect differences between glans/shaft or
  circumcised/uncircumcised for the distance that virions penetrated

## Bigger Plots

## Specific Insights

## Remaining Questions/Concerns

## Conclusions + Future Analysis

- A lot of the emmeans models perform kinda poorly, and there isn’t much
  we can do about this I feel like

# In Depth Analysis

## Virion Particle count

## Mast Cell Count Analysis

## Immunofluorescence Cell Distribution

We also had data from an immunofluorescence experiment, in which these
samples of glans/shaft had been stained with CD3, CD4, and CCR10.
Approximate cell types based on these markers are delineated below:

<details class="code-fold">
<summary>Code</summary>

``` r
tribble(
  ~ "IF Expression", ~ "Expected Cell Type",
  "CD3+CD4+CCR10+", "Th22 (Tissue Homing T cell)",
  "CD3+CCR10+", "CD8+ Cutaneous Lymphoid Antigen T cell",
  "CD4+CCR10+", "Langerhans skin cell",
  "CD3+CD4+", "CD4+ T cell (likely Th1, Th2, Th9, or Th17)",
  "CD3+", "Other T cell",
  "CD4+", "HIV susceptible innate immune cells (macrophages, mast cells, etc)",
  "CCR10+", "Melanocytes"
) |>
  gt::gt() |>
  gt::tab_header(title = "Expected IF Cell Types")
```

</details>
<div id="kpvgtyupby" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#kpvgtyupby table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}
&#10;#kpvgtyupby thead, #kpvgtyupby tbody, #kpvgtyupby tfoot, #kpvgtyupby tr, #kpvgtyupby td, #kpvgtyupby th {
  border-style: none;
}
&#10;#kpvgtyupby p {
  margin: 0;
  padding: 0;
}
&#10;#kpvgtyupby .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}
&#10;#kpvgtyupby .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}
&#10;#kpvgtyupby .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}
&#10;#kpvgtyupby .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}
&#10;#kpvgtyupby .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#kpvgtyupby .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#kpvgtyupby .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}
&#10;#kpvgtyupby .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}
&#10;#kpvgtyupby .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}
&#10;#kpvgtyupby .gt_column_spanner_outer:first-child {
  padding-left: 0;
}
&#10;#kpvgtyupby .gt_column_spanner_outer:last-child {
  padding-right: 0;
}
&#10;#kpvgtyupby .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}
&#10;#kpvgtyupby .gt_spanner_row {
  border-bottom-style: hidden;
}
&#10;#kpvgtyupby .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}
&#10;#kpvgtyupby .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}
&#10;#kpvgtyupby .gt_from_md > :first-child {
  margin-top: 0;
}
&#10;#kpvgtyupby .gt_from_md > :last-child {
  margin-bottom: 0;
}
&#10;#kpvgtyupby .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}
&#10;#kpvgtyupby .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#kpvgtyupby .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}
&#10;#kpvgtyupby .gt_row_group_first td {
  border-top-width: 2px;
}
&#10;#kpvgtyupby .gt_row_group_first th {
  border-top-width: 2px;
}
&#10;#kpvgtyupby .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#kpvgtyupby .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}
&#10;#kpvgtyupby .gt_first_summary_row.thick {
  border-top-width: 2px;
}
&#10;#kpvgtyupby .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#kpvgtyupby .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#kpvgtyupby .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}
&#10;#kpvgtyupby .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}
&#10;#kpvgtyupby .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}
&#10;#kpvgtyupby .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}
&#10;#kpvgtyupby .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#kpvgtyupby .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#kpvgtyupby .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}
&#10;#kpvgtyupby .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}
&#10;#kpvgtyupby .gt_left {
  text-align: left;
}
&#10;#kpvgtyupby .gt_center {
  text-align: center;
}
&#10;#kpvgtyupby .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}
&#10;#kpvgtyupby .gt_font_normal {
  font-weight: normal;
}
&#10;#kpvgtyupby .gt_font_bold {
  font-weight: bold;
}
&#10;#kpvgtyupby .gt_font_italic {
  font-style: italic;
}
&#10;#kpvgtyupby .gt_super {
  font-size: 65%;
}
&#10;#kpvgtyupby .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}
&#10;#kpvgtyupby .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}
&#10;#kpvgtyupby .gt_indent_1 {
  text-indent: 5px;
}
&#10;#kpvgtyupby .gt_indent_2 {
  text-indent: 10px;
}
&#10;#kpvgtyupby .gt_indent_3 {
  text-indent: 15px;
}
&#10;#kpvgtyupby .gt_indent_4 {
  text-indent: 20px;
}
&#10;#kpvgtyupby .gt_indent_5 {
  text-indent: 25px;
}
</style>

| Expected IF Cell Types |                                                                    |
|------------------------|--------------------------------------------------------------------|
| IF Expression          | Expected Cell Type                                                 |
| CD3+CD4+CCR10+         | Th22 (Tissue Homing T cell)                                        |
| CD3+CCR10+             | CD8+ Cutaneous Lymphoid Antigen T cell                             |
| CD4+CCR10+             | Langerhans skin cell                                               |
| CD3+CD4+               | CD4+ T cell (likely Th1, Th2, Th9, or Th17)                        |
| CD3+                   | Other T cell                                                       |
| CD4+                   | HIV susceptible innate immune cells (macrophages, mast cells, etc) |
| CCR10+                 | Melanocytes                                                        |

</div>

Here, looking at the number of cells per sample isn’t the most
interesting, as this is likely just a result of how the tissue was cut.
What *would* be interesting would be to look at if the
***distribution***, i.e. how far any particular cell type was from the
beginning of the epithelial layer, of these cells changed as a result of
the virion infection.

\[INSERT ANALYSIS\]

### Adding Back virion penetration data

We can also then go back to correlate this with virion penetration data:

<details class="code-fold">
<summary>Code</summary>

``` r
load(file = "figures/plots/correlation_plot_1.rda")

correlation_plot_1
```

</details>

![](foreskin_analysis_files/figure-commonmark/correlation-plot-one-1.png)

Here, we can also zoom in on the tip of the biopsy– as the plot above
shows, virions aren’t really penetrating further than 20um. If we just
zoom into this area, we can see what % of each cell type is within our
“viral range” (here defined as 30um from the epithelium), both in the
case where virus was present (4h) and was not (Non infected). *Does the
presence of virions affect the % of immune cell types at the tip?*

<details class="code-fold">
<summary>Code</summary>

``` r
load(file = "figures/plots/correlation_plot_2.rda")

correlation_plot_2
```

</details>

![](foreskin_analysis_files/figure-commonmark/correlation-plot-two-1.png)

Probably not. Maybe something weird (again) happening with the CD4
single positive cells; not sure how relevant this is though and may just
reflect a subset of these cells that we don’t have markers to define.
