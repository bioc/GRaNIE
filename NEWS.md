# GRaNIE 1.1.14 to 1.1.20 (2022-12-13)

## Major changes
- major object changes and optimizations, particularly related to storing the count matrices in an optimized and simpler format. In short, the count matrices are now stored either as normal or sparse matrices, depending on the amount of zeros present. In addition, only the counts after normalization are saved, the raw counts before applying normalization are not stored anymore. If no normalization is wished by the user, as before, the "normalized" counts are equal to the raw counts. `GRaNIE` is now more readily applicable for larger analyses and single-cell analysis even though we just started actively optimizing for it, so we cannot yet recommend applying our framework in a single-cell manner. Older GRN objects are automatically changed internally when executing the major functions upon the first invocation.
- various Documentation and R help updates
- the function `generateStatsSummary` now doesnt alter the stored filtered connections in the object anymore. This makes its usage more intuitive and it can be used anywhere in the workflow.
- removed redundant `biomaRt` calls in the code. This saves time and makes the code less vulnerable to timeout issues caused by remote services
- due to the changes described above, the function `plotPCA_all` now can only plot the normalized counts and not the raw counts anymore (except when no normalization is wanted)
- the GO enrichments are now also storing, for each GO term, the ENSEMBL IDs of the genes that were found in the foreground. This facilitates further exploration of the enrichment results.

## Minor changes
- many small changes in the code

# GRaNIE 1.1.12 and 1.1.13 (2022-09-13)

## Major changes

- many Documentation and R help updates, the [Package Details Vignette](https://grp-zaugg.embl-community.io/GRaNIE/articles/GRaNIE_packageDetails.html#installation) is online
- The workflow vignette is now improved: better figure resolution, figure aspect ratios are optimized, and a few other changes
- the eGRN graph structure as built by `build_eGRN_graph()` in the `GRaNIE` object is now reset whenever the function `filterGRNAndConnectGenes()` is successfully executed to make sure that enrichment functions etc are not using an outdated graph structure. 
- the landing page of the website has been extended and overhauled
- removed some dependency packages and moved others into `Suggests` to lower the installation burden of the package. In addition, removed `topGO` from the `Depends` section (now in `Suggests`) and removed `tidyverse` altogether (before in `Depends`). Detailed explanations when and how the packages listed under `Suggests` are needed can now be found in the new [Package Details Vignette](https://grp-zaugg.embl-community.io/GRaNIE/articles/GRaNIE_packageDetails.html#installation) and are clearly given to the user when executing the respective functions
- major updates to the function `getGRNConnections`, which now has more arguments allowing a more fine-tuned and rich retrieval of eGRN connections, features and feature metadata
- a new function `add_featureVariation` to quantify and interpret multiple sources of biological and technical variation for features (TFs, peaks, and genes) in a GRN object, see the R help for more information 
- `filterGRNAndConnectGenes` now doesnt include feature metadata columns to save space in the result data frame that is created. The help has been updated to make clear that `getGRNConnections` includes these features now.

## Minor changes
- small changes in the GRN object structure, moved `GRN@data$TFs@translationTable` to `GRN@annotation@TFs`. All exported functions run automatically a small helper function to make this change for any GRN object automatically to adapt to the new structure
- many small changes in the code, updated argument checking, and preparing rigorous unit test inclusion
- internally renaming the (recently changed / renamed) gene type `lncRNA` from `biomaRt` to `lincRNA` to be compatible with older versions of `GRaNIE`


# GRaNIE 1.1.X (2022-05-31)

## Minor changes

- added the argument *maxWidth_nchar_plot* to all functions that plot enrichments, and changed the default from 100 to 50. 

## Bug fixes

- fixed a small bug that resulted in the enrichment plots to ignore the value of *maxWidth_nchar_plot*

# GRaNIE 0.99.X (2022-04-26)

## Major changes

- Bioconductor acceptance: this version is the final version for the Bioconductor 3.15 release branch
- full inclusion of the GRN visualization
- extensive vignette updates
- added the possibility to print only particular output pages for all plot functions

## Bug fixes

- various minor bug fixes

## Minor changes

- various minor changes


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
