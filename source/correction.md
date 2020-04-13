
# Usage 2: RNase R effect correction

When you have both RNase R treated and untreated samples, CIRIquant can estimate the before-treatment expression levels of circRNAs detected in RNase R data.

In order to remove RNase R treatment effect, two steps are needed:

1. Run CIRIquant with RNase R treated sample.
2. Run CIRIquant with untreaded total RNA sample, specific `--RNaseR` option using the output gtf file in Step1

Then, CIRIquant will output estimated expression levels of circRNAs detected in RNaseR data, and the header lines will include additional information of RNase R treatment effciency.

## Example usage

```bash
# Step1. Run CIRIquant with RNase R treated data
CIRIquant --config ./hg19.yml \
          -1 ./RNaseR_treated_1.fq.gz \
          -2 ./RNaseR_treated_2.fq.gz \
          --no-gene \
          -o ./RNaseR_treated \
          -p RNaseR_treated \
          -t 6

# Step2. Run CIRIquant with untreated total RNA
CIRIquant --config ./hg19.yml \
          -1 ./TotalRNA_1.fq.gz \
          -2 ./TotalRNA_2.fq.gz \
          -o ./TotalRNA \
          -p TotalRNA \
          -t 6 \
          --RNaseR ./RNaseR_treated/RNaseR_treated.gtf
```