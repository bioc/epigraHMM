test_that("check output from epigraHMM object (bam)",{

    if ("genomationData" %in% rownames(installed.packages())) {
        bamFiles <- system.file("extdata",
                                "wgEncodeBroadHistoneH1hescCtcfStdAlnRep1.chr21.bam",
                                package="genomationData")

        colData <- data.frame(condition = 'CTCF', replicate = 1)

        object <- epigraHMMDataSetFromBam(bamFiles = bamFiles,
                                          colData = colData,
                                          genome = 'hg19',
                                          windowSize = 25000,
                                          gapTrack = TRUE,
                                          blackList = TRUE)
        # Check genome
        expect_true(all(genome(object)=='hg19'))

        # Check window size
        expect_true(sum(width(object)==25000)>=(nrow(object)-1))

        # Check genome via GRanges
        objectOne <- epigraHMMDataSetFromBam(bamFiles = bamFiles,
                                             colData = colData,
                                             genome = rowRanges(object)[1],
                                             windowSize = 25000,
                                             gapTrack = TRUE,
                                             blackList = TRUE)

        expect_equal(assay(objectOne)[1],assay(object)[1])
        expect_equal(width(objectOne),25000)

        # Check discards
        objectTwo <- epigraHMMDataSetFromBam(bamFiles = bamFiles,
                                             colData = colData,
                                             genome = rowRanges(object)[1:3],
                                             windowSize = 25000,
                                             gapTrack = rowRanges(object)[2],
                                             blackList = rowRanges(object)[3])

        expect_equal(assay(objectTwo)[1],assay(objectOne)[1])
        expect_equal(width(objectTwo),25000)
    } else{
        expect_true(TRUE)
    }
})
