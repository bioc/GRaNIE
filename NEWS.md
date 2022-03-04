# GRaNIE 0.15-0.17 (2021-12-13)

## Major changes

- all enrichment analyses have been extended and improved, we added additional ontologies (KEGG, DO, and Reactome), more information in the resulting summary plots
- all plotting functions have been homogenized and extended, PDF width and height can now be set for all exported plot functions. Also, the possibility to not to a PDF but instead to the currently active graphics device is possible. Lastly, setting a different filename is finally possible. Collectively, this provides ultimate flexibility for customizing file names, the output devices used and PDF sizes
- we added a function to build the eGRN network that can be parameterized and that allows future developmemt more easily. Now, network-specific parameters can be changed, such as whether loops should be allowed
- we removed the GRaNIEdev package, the development now happens in a separate branch rather than a different package
- we added Leiden clustering for community clustering (see https://www.nature.com/articles/s41598-019-41695-z for a comparison with louvain)
- extensive vignette updates

## Bug fixes

- various minor bug fixes

## Minor changes

- changed the object structure slightly (graph slot and structure within the stats$enrichment slot)


# GRaNIE 0.9-0.14 (2021-12-13)

## Major changes

- major overhaul and continuous work on peak-gene QC plots
- the *filterData* functions has now more filter parameters, such as filtering for CV. Also, all filters uniformly have a *min* and *max* filter.
- integrated network statistics and various enrichment analyses
- handling of edge cases and rare events in various functions
- packages have been renamed to *GRaNIE* as basename (before: *GRN*)

## Bug fixes

- various minor bug fixes

## Minor changes

- changed the object structure slightly and moved some gene and peak annotation data (such as mean, CV) to the appropriate annotation slot


# GRaNIE 0.8 (2021-05-07)

## Major changes

- improved PCA plotting, PCA plots are now produced for both raw and normalized data
- new filters for the function *filterGRaNIEAndConnectGenes* (*peak_gene.maxDistance*) as well as more flexibility how to adjust the peak-gene raw p-values for multiple testing (including the possibility to use IHW - experimental)
- new function *plotDiagnosticPlots_TFPeaks* for plotting (this function was previously called only internally, but is now properly exported), in analogy to *plotDiagnosticPlots_peakGene*

## Bug fixes

- various minor bug fixes (PCA plotting, compatibility when providing pre-normalized data)

## Minor changes

- changed the object structure slightly and cleaned the config slot, for example
- some functions have been added / renamed to make the workflow more clear and streamlined, see Vignette for details
- some default parameters changed

# GRaNIE 0.7 (2021-03-12)

## Major changes

- improved PCA plotting, also works for pre-normalized counts now when provided as input originally
- more flexibility for data normalization
- homogenized wordings, function calls and workflow clarity, removed unnecessary warnings when plotting peak-gene diagnostic plots, added more R help documentation
- added IHW (Independent Hypothesis Weighting) as a multiple testing procedure for peak-gene p-values in addition to now allowing all methods that are supported by p.adjust

## Bug fixes

- various minor bug fixes

## Minor changes


# GRaNIE 0.6 (2021-02-09)

## Major changes

- significant speed improvements for the peak-FDR calculations and subsequent plotting
- TF-peak diagnostic plots now also show negatively correlated TF-peak statistics irrespective of whether they have been filtered out in the object / pipeline. This may be useful for diagnostic purposes to check whether excluding them is a sensible choice and to confirm the numbers are low

## Bug fixes

- Numbers for connections per correlation bin in the TF-peak diagnostic plots were wrong as they did not correctly differentiate between the different connection types in case multiple ones had been specified (e.g., expression and TF activity). This has been fixed.

## Minor changes


# GRaNIE 0.5 (2021-02-02)

first published package version
