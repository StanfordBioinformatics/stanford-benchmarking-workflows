#!/usr/bin/python

"""
This script will verify that the expected variants in the truth file are in the corresponding VCF file.  If all the
variants in the truth file are found in the VCF file, then this VCF file passes validation.  Otherwise, the VCF file
fails validation.
"""

import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--truth", help="The truth file")
parser.add_argument("-v", "--vcf", help="The vcf file")

args = parser.parse_args()
test_file = args.truth
vcf_file = args.vcf


def get_test_variants():
    """
    This function returns a tuple of a set of tuples of the test variants and the chromosome.

    :return: Set of tuples of (chrom, pos, ref alt, filter) and the chromosome
    :rtype: tuple of set and str
    """
    test_variant_set = set()
    chromosome = ''

    with open(test_file, 'r') as test_obj:
        for line in test_obj:
            if not re.search('^#', line):
                line_items = line.strip().split('\t')

                chrom = line_items[0]
                pos = line_items[1]
                ref = line_items[2]
                alt = line_items[3]
                vcf_filter = line_items[4]

                chromosome = chrom + '\t'
                var = (chrom, pos, ref, alt, vcf_filter)
                test_variant_set.add(var)
            else:
                continue

    return test_variant_set, chromosome


def get_vcf_variants(chromosome):
    """
    This function returns a set of variants with the specified chromosome from the VCF file.

    :param chromosome: The variants' chromosome
    :type chromosome: str

    :return: Set of tuples of (chrom, pos, ref alt, filter)
    :rtype: set
    """
    vcf_variant_set = set()

    with open(vcf_file, 'r') as vcf_obj:
        for line in vcf_obj:
            if re.search('^' + chromosome, line) and not re.search('^#', line):
                line_items = line.strip().split('\t')

                vcf_filter = line_items[6]

                if vcf_filter == 'PASS':
                    chrom = line_items[0]
                    pos = line_items[1]
                    ref = line_items[3]
                    alt = line_items[4]

                    var = (chrom, pos, ref, alt, vcf_filter)
                    vcf_variant_set.add(var)
                else:
                    continue
            else:
                continue

    return vcf_variant_set


########################################################################################################################
#
#   MAIN
#
########################################################################################################################
def main():
    """
    This is the main function.  It prints 'Pass' to the screen if the test variants are found in the VCF variants;
    otherwise, it prints 'Fail'.

    :rtype: void
    """
    test_variant_set, chromosome = get_test_variants()

    vcf_variant_set = get_vcf_variants(chromosome)

    print 'VCF File:', vcf_file.split('/')[-1]
    print 'Truth File:', test_file.split('/')[-1]

    if test_variant_set.issubset(vcf_variant_set):
        print 'Truth Set: ', test_variant_set
        print 'Intersect:', test_variant_set.intersection(vcf_variant_set)
        print 'Pass\n'
    else:
        print 'Truth Set: ', test_variant_set
        print 'Intersect:', test_variant_set.intersection(vcf_variant_set)
        print 'Fail\n'


main()
