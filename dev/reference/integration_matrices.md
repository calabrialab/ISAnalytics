# Example of imported multi-quantification integration matrices.

The data was obtained manually by simulating real research data.

## Usage

``` r
data("integration_matrices")
```

## Format

Data frame with 1689 rows and 8 columns

- chr:

  The chromosome number (as character)

- integration_locus:

  Number of the base at which the viral insertion occurred

- strand:

  Strand of the integration

- GeneName:

  Symbol of the closest gene

- GeneStrand:

  Strand of the closest gene

- CompleteAmplificationID:

  Unique sample identifier

- seqCount:

  Value of the sequence count quantification

- fragmentEstimate:

  Value of the fragment estimate quantification
