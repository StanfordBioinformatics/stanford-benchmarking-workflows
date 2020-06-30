Last Update: May 17, 2019

This is the readme file for the validation scripts: benchmarking_truth_set.py and verify_variants.py.

PLEASE NOTE THAT THESE SCRIPTS USE PYTHON 2.7, THE SUPPORT FOR WHICH WILL STOP BY 2020!!

########################################################################################################################

    benchmarking_truth_set.py

########################################################################################################################

This script produces a tab delimited text file that contains the benchmarking metrics for NIST samples.

The required inputs are the 'number of bases' file and four input files per NIST sample.  The four input files per NIST
sample are in a subdirectory of the input directory.  See directory structure below.

The number of bases file has three columns.  The first column is case name.  The second column is the number of bases in
the whole exome.  The third column is the number bases in the coding exons.  The number of bases file must be prefixed
with 'number_of_bases'.

The input files for each NIST sample must have one file with each of the following suffixes:
    1. CodingExons_indelSizeDistribution.txt
    2. CodingExons.extended.csv
    3. WholeExomeRegions_indelSizeDistribution.txt
    4. WholeExomeRegions.extended.csv

For each sample folder, files (1) and (2) listed above can be found in 'indelSizeDistribution_CodingExons_HappyResults' and 'vcfComparison_by_Happy_CodingExons' sub folders respectively.

Similarly, for each sample folder, files (3) and (4) listed above can be found in 'indelSizeDistribution_WholeExomeRegions_HappyResults' and 'vcfComparison_by_Happy_WholeExomeRegions' sub folders respectively.

The input directory MUST have the following structure.  Replace '*' with the rest of the file name.

Input files:
    input_dir/
        number_of_bases*.txt
        NA12878/
            *CodingExons.extended.csv
            *CodingExons_indelSizeDistribution.txt
            *WholeExomeRegions.extended.csv
            *WholeExomeRegions_indelSizeDistribution.txt
        NA24143/
            *CodingExons.extended.csv
            *CodingExons_indelSizeDistribution.txt
            *WholeExomeRegions.extended.csv
            *WholeExomeRegions_indelSizeDistribution.txt
        NA24149/
            *CodingExons.extended.csv
            *CodingExons_indelSizeDistribution.txt
            *WholeExomeRegions.extended.csv
            *WholeExomeRegions_indelSizeDistribution.txt
        NA24385/
            *CodingExons.extended.csv
            *CodingExons_indelSizeDistribution.txt
            *WholeExomeRegions.extended.csv
            *WholeExomeRegions_indelSizeDistribution.txt
        NA24631/
            *CodingExons.extended.csv
            *CodingExons_indelSizeDistribution.txt
            *WholeExomeRegions.extended.csv
            *WholeExomeRegions_indelSizeDistribution.txt


The columns of the output file are Case, Number of bases, Truth total, TP, FP, FN, TN = TotalBases - (TP + FN + FP),
TotalNegative  = TN + FP, NPA = TN/(Total Negative), Precision and Recall.

Output file:
    Final_benchmarking_metrics_YYYY-MM-DD.txt


usage: benchmarking_truth_set.py [-h] [-i INPUT] [-o OUTPUT]

optional arguments:
  -h, --help                    Show this help message and exit
  -i INPUT, --input INPUT       The input directory; default is current directory
  -o OUTPUT, --output OUTPUT    The output directory; default is current directory


########################################################################################################################

    verify_variants.py

########################################################################################################################

This script verifies that the expected variants in the truth file are in the corresponding sample VCF file.  If all
the variants in the truth file are found in the VCF file, then this VCF file passes validation.  Otherwise, the VCF file
fails validation.  This script produces screen output.

The required inputs are either a GVCF or VCF file and a truth file.

Input files:
    *.recal.annotated.final.g.vcf (single sample * = sample ID)
    *.recal.annotated.final.vcf (trio with * = family ID)
    expected_variants/*_Truth.txt


Output to screen:
    The VCF file name
    The truth file name
    The set of variants in the truth file
    The intersect of variants in the truth file and the VCF file
    Pass | Fail


usage: verify_variants.py [-h] [-d] [-t TRUTH] [-v VCF]

optional arguments:
  -h, --help                Show this help message and exit
  -t TRUTH, --truth TRUTH   The truth file
  -v VCF, --vcf VCF         The GVCF or VCF file
