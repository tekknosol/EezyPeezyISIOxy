#!/usr/bin/env Rscript

# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.

# targets::tar_make()
# targets::tar_make_clustermq(workers = 300, reporter = "summary") # nolint

targets::tar_make(
  seconds_meta_append = 15,
  seconds_reporter = 0.5
)

#targets::tar_make_clustermq(workers = 32, reporter = "silent")

