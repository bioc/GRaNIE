test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("addConnections_TF_peak", {
  options(timeout = 500)
  GRN = loadExampleObject(forceDownload = TRUE)
  
  GRN = addConnections_TF_peak(GRN, connectionTypes = c("expression"), outputFolder = ".", forceRerun = TRUE)
  expect_s4_class(GRN, "GRN")
  
  expect_error(addConnections_TF_peak(GRN, connectionTypes = c("bla"), outputFolder = ".", forceRerun = TRUE), 
               regexp = "Assertion on")
  
  expect_error(addConnections_TF_peak(GRN, connectionTypes = c("expression", "TF_activity"), outputFolder = ".", forceRerun = TRUE))
  
  expect_error(addConnections_TF_peak(GRN, connectionTypes = c("expression", "TF_activity"), outputFolder = ".", forceRerun = TRUE, removeNegativeCorrelation = c(FALSE)), 
               regexp = "Assertion on")
  
  # Make sure TF activity is added now
  # TODO: Add TF activity here
  # GRN = addConnections_TF_peak(GRN, connectionTypes = c("expression", "TF_activity"), removeNegativeCorrelation = c(FALSE, TRUE), forceRerun = TRUE)
  # expect_s4_class(GRN, "GRN")
  # 
  # GRN = addConnections_TF_peak(GRN, connectionTypes = c("expression", "TF_activity"), removeNegativeCorrelation = c(TRUE, FALSE), forceRerun = TRUE)
  # expect_s4_class(GRN, "GRN")
  # 
  # GRN = addConnections_TF_peak(GRN, connectionTypes = c("TF_activity"), removeNegativeCorrelation = c(TRUE), forceRerun = TRUE)
  # expect_s4_class(GRN, "GRN")
  

})

test_that("AR wrapper", {
    
    options(timeout = 500)
    GRN = loadExampleObject()

    expect_error(AR_classification_wrapper(GRN, significanceThreshold_Wilcoxon = 1.5, outputFolder = "."),  regexp = "Assertion on")
    
    # Run with TF activity and expression
    GRN = AR_classification_wrapper(GRN, significanceThreshold_Wilcoxon = 0.05,
                                            plot_minNoTFBS_heatmap = 100, plotDiagnosticPlots = TRUE,
                                            outputFolder = ".", 
                                            deleteIntermediateData = TRUE, forceRerun = TRUE)
    
    # What happens if TF activity is defined but addConnections_TF_peak with only expression?
    
    
})
