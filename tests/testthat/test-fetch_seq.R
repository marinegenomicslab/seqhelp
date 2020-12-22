test_that("fasta access is correct", {
  fa <- system.file("testdata", "tiny.fasta", package="seqhelp")
  tinyseq <- fetch_seq(fasta = fa, seq_name = "seq1", seq_start = 1, seq_end = 3)

  expect_equal(tinyseq, "ATC")
})
