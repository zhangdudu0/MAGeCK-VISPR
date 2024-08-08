# Change Log

All notable changes to this project will be documented in this file.


## Changed

## [0.5.6] - 2020-12-04
### Changed
- Allowing multiple sgRNAs to be annotated in one library
- Separate additional-rra-mle-parameter into two entries (additional-mle-parameter and additional-rra-parameter) in config.yaml.

## [0.5.5] - 2020-09-04
### Changed
- Add indicator to whether turn on sgRNA annotation
- A BED file indicating the log fold change for each sgRNA is generated, When annotation is turned on and either 1) it's an RRA experiment or 2) a day0 label is specified
- Simplify rules for vispr for paired-end read support

### Added
- Add annotation support for MAGeCK-VISPR to add custom values in the bed file
- Add copy number variation (CNV) correction support

## [0.5.3] - 2016-09-08
### Changed
- Prefer mle rule over rra if design matrix is present.


## [0.5.1] - 2016-03-10
### Added
- It is now possible to activate batch effect correction via combat/sva. The mechanism is described in the example config file that is obtained via "mageck-vispr init".

