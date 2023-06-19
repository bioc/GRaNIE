# GRaNIE 1.5.1 (2023-06-19)

- version jump due to new Bioconductor development cycle

## New features and stability improvements
- we replaced `biomaRt` for the full genome annotation retrieval in `addData` with a different approach that is more reliable, as we had more and more issues with `biomaRt` in the recent past. While using the old `biomaRt` approach is still an option, the default is now to use the `AnnotationHub` package from Bioconductor. This makes GRaNIE overall more stable and less reliant on `biomaRt` due to the strict timeouts and query size restrictions.

# GRaNIE 1.3.36 (2023-05-09)

## New features and vignette updates
- a new correlation method has been implemented that replaces the "old" robust method (via `addRobustRegression`) that was available as an experimental feature until now. It is implemented in the `WGCNA` package and called *biweight midcorrelation* or short *bicor*, a robust type of correlation based on medians that can be used as an alternative to Spearman correlation. Biweight midcorrelation has been shown to be more robust in evaluating similarity in gene expression networks and is often used for weighted correlation network analysis. In addition, this new correlation type can now also be selected for the TF-peak correlations, which was not possible before. Lastly, the code has been cleaned and simplified and all instances of `addRobustRegression` have been removed and replaced by a new third option `bicor` (in addition to Pearson and Spearman, as before) for the `corMethod` argument in multiple functions that support this feature. All vignettes have been updated accordingly.
- in `addData`, when a `DESeq2` size factor normalization is selected by the user, it is now explicitly checked whether enough genes are available that contain no 0 values on which the size factor normalization is based on. If this is not the case (the default and hard-coded limit is currently set to a minimum of 100), an error is thrown. This becomes particularly relevant for single-cell derived data with a high fraction of 0s, and prevents a normalization based on very few genes and improves the error messages that `DESeq2` throws otherwise for an improved user experience.


# GRaNIE 1.3.35 (2023-05-08)

## Bugfixes and stability improvements
- improved the stability of the `biomaRt` call, which did not work as originally intended in case of temporary connection failures. Now, calls to `biomaRt` are attempted up to 40 times to increase the chances of not suffering from connection issues. Also, the approach to deal with `BiocParallel` failures has been changed.

## Paper acceptance and publication update
- **We are particularly excited to announce that the GRaNIE and GRaNPA paper is finally out! See [here for the final article in Molecular Systems Biology](https://www.embopress.org/doi/full/10.15252/msb.202311627).**

# GRaNIE 1.3.34 (2023-03-29)

## Bugfixes
- many small bugfixes and other small improvements to homogenize the user experience due to the usage of systematic unit tests


# GRaNIE 1.3.33 (2023-03-02-2023-03-06)

## New features and vignette updates
- we provide two new functions with this update:
    1. `getGRNSummary()` that summarizes a `GRN` object and returns a named list, which can be used to compare different `GRN` objects ore easily among each other, for example.
    2. `plotCorrelations()` for scatter plots of the underlying data for either TF-peak, peak-gene or TF-gene pairs. This can be useful to visualize specific TF-peak, peak-gene or TF-gene pairs to investigate the underlying data and to judge the reasonability of the inferred connection.
    
- methods vignette updates

# GRaNIE 1.3.31-1.3.32 (2023-03-02-2023-03-06)

## Bugfixes
- various small bugfixes that were accidentally introduced in the latest change from using the `TF.ID` instead of `TF.name` column as unique TF identifier


# GRaNIE 1.3.26-1.3.30 (2023-03-02-2023-03-06)

## New features and vignette updates
- added two new supported genomes: `rn6`/`rn7` and `dm6` for the rat and the Drosophila (fruit fly) genome, respectively
- **added preliminary support for a new, alternative way of how to import TF and TFBS data into `GRaNIE`. We now additionally offer a more user-friendly way by making it possible to directly use the `JASPAR2022` database. You do not need any custom files anymore for this approach! [See the Package vignette for more details](https://grp-zaugg.embl-community.io/GRaNIE/articles/GRaNIE_packageDetails.html#input_TF).**

## Bugfixes
- fixed a regression bug in `addConnections_TF_peak` (`Column `peak.GC.class` doesn't exist.`) that was caused due to the recent GC modifications

# GRaNIE 1.3.25 (2023-02-20)

## New features and vignette updates
- additional significant methods vignette updates
- updates and clarifications for the workflow vignette
- a new QC plot for `plotDiagnosticPlots_TFPeaks` (and indirectly in `addConnections_TF_peak` when `plotDiagnosticPlots = TRUE`) on page 1 that shows the total number of connections for real and background TF-peak links as calculated and stored in the `GRN` object, stratified by TF-peak FDR and correlation bin. This is a similar plot as we show in the paper and helps comparing foreground and background.

## Improvements
- speed improvements for `plotDiagnosticPlots_TFPeaks` (and indirectly in `addConnections_TF_peak` when `plotDiagnosticPlots = TRUE`) when `plotAsPDF = FALSE`

## Bugfixes
- fixed a bug that  only occurred in `addConnections_TF_peak` when using `useGCCorrection = TRUE`

# GRaNIE 1.3.24 (2023-02-20)

## New features and vignette updates
- significant methods vignette updates that help clarifying methods details


# GRaNIE 1.3.22-1.3.23 (2023-02-09-2023-02-15)

## Minor changes
- Small workflow vignette updates 

## Bug fixes
- we were informed that newer versions of `dplyr` (1.1.0) changed their default behavior for the function `if_else` when `NULL` is involved, which caused an error. We changed the implementation to accommodate for that and now avoid `dplyr::if_else` and use base R `ifelse` instead.

# GRaNIE 1.3.18-1.3.21 (2023-02-01-2023-02-07)

## Minor changes
- Small vignette updates and fixing typos / improved wording

## Bug fixes
- due to a change from USCS that affected `GenomeInfoDb::getChromInfoFromUCSC("hg38")` (see [here](https://github.com/Bioconductor/GenomeInfoDb/issues/82) for more details), the minimum required version of `GenomeInfoDb` had to be increased to `1.34.8`. If you have troubles installing at least this version, we recommend updating to the newest Bioconductor version 3.16 or (without warranties) use the following line to manually install the newest version directly from GitHub outside of Bioconductor (not recommended): `BiocManager::install("Bioconductor/GenomeInfoDb)"` 
- small change in `addData()` so that peak IDs are stored with the same name in the object in case the user-provided peak IDs have the format `chr:start:end` as opposed to the required `chr:start-end`. `filterData()` otherwise incorrectly discarded all peaks because of the ID mismatch caused by the two different formats.
- fixed a rare edge case in `filterGRNAndConnectGenes()` that caused an error when 0 TF-peak connections were found beforehand

# GRaNIE 1.3.17 (2023-01-26)

## New features
- We are excited to announce that we added a new vignette for how to use `GRaNIE` for single-cell data! We plan to update it regularly with new information. [Check it out here!](https://grp-zaugg.embl-community.io/GRaNIE/articles/GRaNIE_singleCell_eGRNs.html)

# GRaNIE 1.3.16 (2023-01-24)

## New features
- significant updated to the package details vignette
- revisited and improved the internal logging and object history. The time when a function was called is now added to the list name, which allows the storage of multiple instances of the same function.
- new parameter in `addData()`: `geneAnnotation_customHost` to specify a custom host and overriding the default and previously hard-coded hostname when retrieving gene annotation data via `biomaRt`.
- the function `getGRNConnections()` can now also include the various additional metadata for all `type` parameters and not only the default type `all.filtered`.

# GRaNIE 1.3.15 (2023-01-20)

## Bug fixes
- fixed an error that appeared in rare cases when a chromosome name from either peak or RNA data could not be found in `biomaRt` such as `GL000194.1`. Peaks from chromosomes with irretrievable lengths are now automatically discarded.
- significant updates to the package details vignette

# GRaNIE 1.3.13-1.3.14 (2023-01-20)

## New features
- the function `plotDiagnosticPlots_peakGene()` (which is also called indirectly from `addConnections_peak_gene()` when setting `plotDiagnosticPlots = TRUE`) now stores the plot data for the QC plots from the first page into the GRN object. It is stored in `GRN@stats$peak_genes` 
- the columns of the result table from `getGRNConnections()` are now explained in detail in the R help, and we reference this from the Vignette and other places
- various significant Vignette updates

## Bug fixes
- optimized the column names for the function `getGRNConnections()`, which now does not return duplicate columns for particular cases anymore
- improved printing in the log for the function `filterData()` and `addData()`
- the `loadExampleObject()` function has been optimized and should now force download an example object when requesting it.
- the package version as stored in the GRN object now works correctly.

## Minor changes
- further code cleaning in light of the `tidyselect` changes in version 1.2.0 to eliminate deprecated warnings
- the default gene types for `addConnections_peak_gene()` and `plotDiagnosticPlots_peakGene()` have been homogenized and changed to `list(c("all"), c("protein_coding"))`. Before, the default was `list(c("protein_coding", "lincRNA"))`, but we decided to now split this into two separate lists: Once for all genes irrespective of the gene type and once for only protein-coding genes. As before, `lincRNA` or other gene types can of course still be selected and chosen.
- various minor changes


# GRaNIE 1.3.12 (2022-12-22)

## Bug fixes
- bug fix in `plotCommunitiesEnrichment()` that was introduced due to the `tidyselect` 1.2.0 changes

## Minor changes
- further code cleaning in light of the `tidyselect` changes in version 1.2.0 to eliminate deprecated warnings

# GRaNIE 1.3.11 (2022-12-16)

## Major changes
- the default URL for the example `GRN` object in `loadExampleObject()` had to be changed due to changes in the IT infrastructure. The new stable default URL is now \url{https://git.embl.de/grp-zaugg/GRaNIE/-/raw/master/data/GRN.rds}, in the same Git repository that provides `GRaNIE` outside of Bioconductor.

## Bug fixes
- fixing bugs introduced due to the tidyverse 1.2.0 related code cleaning
- other bugfix accidentally introduced in the previous commits

# GRaNIE 1.3.10 (2022-12-15)

## Bug fixes
- revisited the import of TADs and made the code more error-prone and fixed some bugs related to TADs. Importing TADs now works again as before.

## Minor changes
- code cleaning in light of the `tidyselect` changes in version 1.2.0 to eliminate deprecated warnings

## New features
- new argument for `addConnections_peak_gene()`: `TADs_mergeOverlapping`. See the R help for more details.

# GRaNIE 1.3.9 (2022-12-14)

## New features
- new argument for `addConnections_peak_gene()`: `shuffleRNACounts`. See the R help for more details.

## Minor changes
- first round of code cleaning in light of the `tidyselect` changes in version 1.2.0 to eliminate deprecated warnings

# GRaNIE 1.3.4-1.3.8 (2022-12-06)

## Major changes
- the `topGO` package is now  required package and not optional anymore. The reasoning for this is that the standard vignette should run through with the default arguments, and `GO` annotation is the default ontology so `topGO` is needed for this. Despite this package still being optional from a strict workflow point of view, we feel this is a better way and improves user friendliness by not having to install another package in the middle of the workflow.

## Minor changes
- in `initializeGRN()`, the `objectMetadata` argument is now checked whether it contains only atomic elements, and an error is thrown if this is not the case. As this list is not supposed to contain real data, checking this prevents the print(GRN) function to unnecessarily print the whole content of the provided object metadata, thereby breaking the original purpose.

## New features
- `addTFBS()` got two more arguments to make it more flexible. Now, it is possible to specify the file name of the translation table to be used via the argument `translationTable`, which makes it more flexible than the previously hard-coded name `"translationTable.csv`. In addition, the column separator for this file can now be specified via the argument `translationTable_sep`
- Overlapping TFBS data with the peak is now more error-tolerant and does not error out in case that some chromosome or contig names from the TFBS BED files contain elements the size of which cannot be retrieved online. This was the case for some contig names with the suffix `decoy`, for example. If such elements are found, a warning is now thrown and they are ignored as they are usually not wanted anyway.
- in case a GRN objects contains 0 connections (e..g, because of too strict filtering), subsequent functions as well as the `print` function now give a more user-friendly warning / error message.

# GRaNIE 1.1.22-1.3.3 (2022-11-29)

## New features
- additional normalization schemes have been implemented, including GC-aware normalization schemes for peaks, and existing normalization methods have been renamed for clarity. See `?addData` for details.
- further reduced the package burden; the large genome annotation packages are now more or less fully optional and only needed when a GC-aware normalization has been chosen or when additional peak annotation is wanted. However, in contrast to before, none of these annotation packages are strictly required anywhere anymore. The vignettes have been updated accordingly.

## Minor changes
- various small changes in the code
- vignette updates


# GRaNIE 1.1.14-1.1.21 (2022-11-13)

## Major changes and new features
- major object changes and optimizations, particularly related to storing the count matrices in an optimized and simpler format. In short, the count matrices are now stored either as normal or sparse matrices, depending on the amount of zeros present. In addition, only the counts after normalization are saved, the raw counts before applying normalization are not stored anymore. If no normalization is wished by the user, as before, the "normalized" counts are equal to the raw counts. `GRaNIE` is now more readily applicable for larger analyses and single-cell analysis even though we just started actively optimizing for it, so we cannot yet recommend applying our framework in a single-cell manner. Older GRN objects are automatically changed internally when executing the major functions upon the first invocation.
- various Documentation and R help updates
- the function `generateStatsSummary()` now doesnt alter the stored filtered connections in the object anymore. This makes its usage more intuitive and it can be used anywhere in the workflow.
- removed redundant `biomaRt` calls in the code. This saves time and makes the code less vulnerable to timeout issues caused by remote services
- due to the changes described above, the function `plotPCA_all()` now can only plot the normalized counts and not the raw counts anymore (except when no normalization is wanted)
- the GO enrichments are now also storing, for each GO term, the ENSEMBL IDs of the genes that were found in the foreground. This facilitates further exploration of the enrichment results.

## Minor changes
- many small changes in the code

# GRaNIE 1.1.12 and 1.1.13 (2022-09-13)

## Major changes and new features

- many Documentation and R help updates, the [Package Details Vignette](https://grp-zaugg.embl-community.io/GRaNIE/articles/GRaNIE_packageDetails.html#installation) is online
- The workflow vignette is now improved: better figure resolution, figure aspect ratios are optimized, and a few other changes
- the eGRN graph structure as built by `build_eGRN_graph()` in the `GRaNIE` object is now reset whenever the function `filterGRNAndConnectGenes()` is successfully executed to make sure that enrichment functions etc are not using an outdated graph structure. 
- the landing page of the website has been extended and overhauled
- removed some dependency packages and moved others into `Suggests` to lower the installation burden of the package. In addition, removed `topGO` from the `Depends` section (now in `Suggests`) and removed `tidyverse` altogether (before in `Depends`). Detailed explanations when and how the packages listed under `Suggests` are needed can now be found in the new [Package Details Vignette](https://grp-zaugg.embl-community.io/GRaNIE/articles/GRaNIE_packageDetails.html#installation) and are clearly given to the user when executing the respective functions
- major updates to the function `getGRNConnections()`, which now has more arguments allowing a more fine-tuned and rich retrieval of eGRN connections, features and feature metadata
- a new function `add_featureVariation()` to quantify and interpret multiple sources of biological and technical variation for features (TFs, peaks, and genes) in a GRN object, see the R help for more information 
- `filterGRNAndConnectGenes()` now doesnt include feature metadata columns to save space in the result data frame that is created. The help has been updated to make clear that `getGRNConnections()` includes these features now.

## Minor changes
- small changes in the GRN object structure, moved `GRN@data$TFs@translationTable` to `GRN@annotation@TFs`. All exported functions run automatically a small helper function to make this change for any GRN object automatically to adapt to the new structure
- many small changes in the code, updated argument checking, and preparing rigorous unit test inclusion
- internally renaming the (recently changed / renamed) gene type `lncRNA` from `biomaRt` to `lincRNA` to be compatible with older versions of `GRaNIE`


# GRaNIE 1.1.X (2022-05-31)

## New features

- added the argument *maxWidth_nchar_plot* to all functions that plot enrichments, and changed the default from 100 to 50. 

## Bug fixes

- fixed a small bug that resulted in the enrichment plots to ignore the value of *maxWidth_nchar_plot*

# GRaNIE 0.99.X (2022-04-26)

## Major changes and new features

- Bioconductor acceptance: this version is the final version for the Bioconductor 3.15 release branch
- full inclusion of the GRN visualization
- extensive vignette updates
- added the possibility to print only particular output pages for all plot functions

## Bug fixes

- various minor bug fixes

## Minor changes

- various minor changes


# GRaNIE 0.15-0.17 (2021-12-13)

## Major changes and new features

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

## Major changes and new features

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

## Major changes and new features

- improved PCA plotting, PCA plots are now produced for both raw and normalized data
- new filters for the function `filterGRaNIEAndConnectGenes()` (`peak_gene.maxDistance`) as well as more flexibility how to adjust the peak-gene raw p-values for multiple testing (including the possibility to use IHW - experimental)
- new function `plotDiagnosticPlots_TFPeaks()` for plotting (this function was previously called only internally, but is now properly exported), in analogy to `plotDiagnosticPlots_peakGene()`

## Bug fixes

- various minor bug fixes (PCA plotting, compatibility when providing pre-normalized data)

## Minor changes

- changed the object structure slightly and cleaned the config slot, for example
- some functions have been added / renamed to make the workflow more clear and streamlined, see Vignette for details
- some default parameters changed

# GRaNIE 0.7 (2021-03-12)

## Major changes and new features

- improved PCA plotting, also works for pre-normalized counts now when provided as input originally
- more flexibility for data normalization
- homogenized wordings, function calls and workflow clarity, removed unnecessary warnings when plotting peak-gene diagnostic plots, added more R help documentation
- added IHW (Independent Hypothesis Weighting) as a multiple testing procedure for peak-gene p-values in addition to now allowing all methods that are supported by p.adjust

## Bug fixes

- various minor bug fixes

## Minor changes


# GRaNIE 0.6 (2021-02-09)

## Major changes and new features

- significant speed improvements for the peak-FDR calculations and subsequent plotting
- TF-peak diagnostic plots now also show negatively correlated TF-peak statistics irrespective of whether they have been filtered out in the object / pipeline. This may be useful for diagnostic purposes to check whether excluding them is a sensible choice and to confirm the numbers are low

## Bug fixes

- Numbers for connections per correlation bin in the TF-peak diagnostic plots were wrong as they did not correctly differentiate between the different connection types in case multiple ones had been specified (e.g., expression and TF activity). This has been fixed.

## Minor changes


# GRaNIE 0.5 (2021-02-02)

first published package version
