
# Usage 2: RNase R effect correction

In order to remove effect for RNase R treatment, two steps of programs are needed

1. Run CIRIquant with RNase R treated sample
2. Use output gtf file in Step1 and run CIRIquant with `--RNaseR` option using output gtf in previous step

The output is in the same format as normal run, however the header line is appended with additional 
information of RNase R treatment

## Example usage

```
CIRIquant --config ./chr1.yml -1 ./test_1.fq.gz -2 ./test_2.fq.gz --no-gene -o ./test -p test -t 6
```