test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("addConnections_TF_peak", {
  GRN = loadExampleObject()
  
  GRN = addConnections_TF_peak(GRN, connectionTypes = c("expression"), forceRerun = TRUE)
  expect_s4_class(GRN, "GRN")
  
  expect_error(addConnections_TF_peak(GRN, connectionTypes = c("bla"), forceRerun = TRUE), 
               regexp = "Assertion on")
  
  expect_error(addConnections_TF_peak(GRN, connectionTypes = c("expression", "TF_activity"), forceRerun = TRUE))
  
  expect_error(addConnections_TF_peak(GRN, connectionTypes = c("expression", "TF_activity"), forceRerun = TRUE, removeNegativeCorrelation = c(FALSE)), 
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
    
    GRN = loadExampleObject()

    expect_error(AR_classification_wrapper(GRN, significanceThreshold_Wilcoxon = 1.5), regexp = "Assertion on")
    
    # Run with TF activity and expression
    GRN = AR_classification_wrapper(GRN, significanceThreshold_Wilcoxon = 0.05,
                                            plot_minNoTFBS_heatmap = 100, plotDiagnosticPlots = TRUE,
                                            deleteIntermediateData = TRUE, forceRerun = TRUE)
    
    # What happens if TF activity is defined but addConnections_TF_peak with only expression?
    
    
})
