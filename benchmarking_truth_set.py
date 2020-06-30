#!/usr/bin/python

"""
This script creates a file with the benchmarking metrics for the NIST sample truth sets.
"""

import argparse
import datetime
import os
import sys

create_date = str(datetime.date.today())

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="The input directory; default is current directory", default='.')
parser.add_argument("-o", "--output", help="The output directory; default is current directory", default='.')

args = parser.parse_args()
input_dir = args.input
output_dir = args.output

output_file = 'Final_benchmarking_metrics_' + create_date + '.txt'


def get_tn(num_bases, tp, fp, fn):
    """
    This functions returns the number of true negatives.

    :param num_bases: Number of bases
    :type num_bases: int
    :param tp: Number of true positives
    :type tp: int
    :param fp: Number of false positives
    :type fp: int
    :param fn: Number of false negatives
    :type fn: int

    :return: Number of true negatives
    :rtype: int
    """
    return num_bases - (tp + fp + fn)


def get_total_negative(fp, tn):
    """
    This functions returns the number of total negatives.

    :param fp: Number of false positives
    :type fp: int
    :param tn: Number of true negatives
    :type tn: int

    :return: Number of total negative
    :rtype: int
    """
    return fp + tn


def get_npa(fp, tn):
    """
    This functions returns the negative percent agreement.

    :param fp: Number of false positives
    :type fp: int
    :param tn: Number of true negatives
    :type tn: int

    :return: The negative percent agreement, NPA = (tn) / (fp+tn)
    :rtype: float
    """
    if fp + tn == 0:
        return 'NaN'
    else:
        return format(round(100 * (float(tn) / (fp + tn)), 2), '.2f')


def get_precision(tp, fp):
    """
    This function returns the precision for indels by size.

    :param tp: Number of true positives
    :type tp: int
    :param fp: Number of false positives
    :type fp: int

    :return: The precision, (tp) / (tp+fp)
    :rtype: float
    """
    if tp + fp == 0:
        return 'NaN'
    else:
        return format(round(100 * (float(tp) / (tp + fp)), 2), '.2f')


def get_recall(tp, fn):
    """
    This function returns the recall for indels by size.

    :param tp: Number of true positives
    :type tp: int
    :param fn: Number of false negatives
    :type fn: int

    :return: The recall, (tp) / (tp+fn)
    :rtype: float
    """
    if tp + fn == 0:
        return 'NaN'
    else:
        return format(round(100 * (float(tp) / (tp + fn)), 2), '.2f')


def get_indel_by_size(path, filename, case_name, num_bases):
    """
    This function generates the output data from the indel by size text file.

    :param path: The path to the input TXT file.
    :type path: str
    :param filename: The file name
    :type filename: str
    :param case_name: The case name
    :type case_name: str
    :param num_bases: Number of bases
    :type num_bases: int

    :return: A list of indel size values for the output file.
    :rtype: list
    """
    indels_11_20 = list()
    indels_21_50 = list()

    tp_1_10_list = list()
    fp_1_10_list = list()
    fn_1_10_list = list()

    with open(os.path.join(path, filename), 'r') as infile_obj:  # indel txt file
        for line in infile_obj:
            if line.startswith('I') or line.startswith('51'):
                continue
            else:
                line_items = line.split('\t')
                tp = int(line_items[2])
                fp = int(line_items[3])
                fn = int(line_items[4])

                if line.startswith('1\t'):
                    tp_1_10_list.append(tp)
                    fp_1_10_list.append(fp)
                    fn_1_10_list.append(fn)

                elif line.startswith('2 - 5\t'):
                    tp_1_10_list.append(tp)
                    fp_1_10_list.append(fp)
                    fn_1_10_list.append(fn)

                elif line.startswith('6 - 10\t'):
                    tp_1_10_list.append(tp)
                    fp_1_10_list.append(fp)
                    fn_1_10_list.append(fn)

                elif line.startswith('11 - 20\t'):
                    truth_total = tp + fn
                    tn = get_tn(num_bases, tp, fp, fn)
                    total_negative = get_total_negative(fp, tn)
                    npa = get_npa(fp, tn)
                    precision = get_precision(tp, fp)
                    recall = get_recall(tp, fn)

                    indels_11_20 = [case_name + ' - Indels 11 - 20', str(num_bases),
                                    str(truth_total), str(tp), str(fp), str(fn), str(tn),
                                    str(total_negative), str(npa), str(precision), str(recall)]

                elif line.startswith('21 - 50\t'):
                    truth_total = tp + fn
                    tn = get_tn(num_bases, tp, fp, fn)
                    total_negative = get_total_negative(fp, tn)
                    npa = get_npa(fp, tn)
                    precision = get_precision(tp, fp)
                    recall = get_recall(tp, fn)

                    indels_21_50 = [case_name + ' - Indels 21 - 50', str(num_bases),
                                    str(truth_total), str(tp), str(fp), str(fn), str(tn),
                                    str(total_negative), str(npa), str(precision), str(recall)]
                else:
                    continue

        # add indel size 1, 2-5, and 6-10
        tp_1_10 = sum(tp_1_10_list)
        fp_1_10 = sum(fp_1_10_list)
        fn_1_10 = sum(fn_1_10_list)

        # do the calculations for indel size 1-10
        truth_total = tp_1_10 + fn_1_10
        tn = get_tn(num_bases, tp_1_10, fp_1_10, fn_1_10)
        total_negative = get_total_negative(fp_1_10, tn)
        npa = get_npa(fp_1_10, tn)
        precision = get_precision(tp_1_10, fp_1_10)
        recall = get_recall(tp_1_10, fn_1_10)

        indels_1_10 = [case_name + ' - Indels 1 - 10', str(num_bases),
                       str(truth_total), str(tp_1_10), str(fp_1_10), str(fn_1_10), str(tn),
                       str(total_negative), str(npa), str(precision), str(recall)]

    return [indels_1_10, indels_11_20, indels_21_50]


def get_column_indexes(column_name_items):
    """
    This function returns the indexes for the columns of interest from the CSV file.

    :param column_name_items: List of column names
    :type column_name_items: list

    :return: Column index for 'TRUTH.TOTAL', 'QUERY.TP', 'QUERY.FP', 'TRUTH.FN', 'METRIC.Precision' and 'METRIC.Recall'
    :rtype: list
    """
    total_index = column_name_items.index('TRUTH.TOTAL')
    tp_index = column_name_items.index('QUERY.TP')
    fp_index = column_name_items.index('QUERY.FP')
    fn_index = column_name_items.index('TRUTH.FN')
    precision_index = column_name_items.index('METRIC.Precision')
    recall_index = column_name_items.index('METRIC.Recall')

    return [total_index, tp_index, fp_index, fn_index, precision_index, recall_index]


def get_csv_data(indexes, data_items, case_name, num_bases):
    """
    This function returns a list of the data from the CSV file that will go into the output file.

    :param indexes: Column indexes for TRUTH.TOTAL, QUERY.TP, QUERY.FP, TRUTH.FN, METRIC.Precision and METRIC.Recall
    :type indexes: list
    :param data_items: Line items from the CSV file
    :type data_items: list
    :param case_name: The Truth Set name
    :type case_name: str
    :param num_bases: Number of bases in the Truth Set
    :type num_bases: int

    :return: The data from the CSV file that will go into the output file
    :rtype: list
    """
    total_index, tp_index, fp_index, fn_index, precision_index, recall_index = indexes
    truth_total = int(float(data_items[total_index]))  # truth.total
    tp = int(float(data_items[tp_index]))  # query.tp
    fp = int(float(data_items[fp_index]))  # query.fp
    fn = int(float(data_items[fn_index]))  # truth.fn
    precision = round(float(data_items[precision_index]) * 100, 2)  # metric.precision
    recall = round(float(data_items[recall_index]) * 100, 2)  # metric.recall

    tn = get_tn(num_bases, tp, fp, fn)
    total_negative = get_total_negative(fp, tn)
    npa = get_npa(fp, tn)

    return [case_name, str(num_bases),
            str(truth_total), str(tp), str(fp), str(fn),
            str(tn), str(total_negative), str(npa), str(precision), str(recall)]


def get_indel_and_snp(path, filename, case_name, num_bases):
    """
    This function returns a tuple of lists of the InDel and SNP output data from the CSV file.

    :param path: The path to the input CSV file.
    :type path: str
    :param filename: The file name
    :type filename: str
    :param case_name: The case name
    :type case_name: str
    :param num_bases: Number of bases
    :type num_bases: int

    :return: A tuple of the list of values for the indels and the list of values for the snps
    :rtype: tuple
    """
    with open(os.path.join(path, filename), 'r') as infile_obj:
        column_name_items = list()
        indel_items = list()
        snp_items = list()

        for line in infile_obj:
            if line.strip().startswith(',METRIC.'):  # column header line
                column_name_items = line.strip().split(',')
            elif line.startswith('Locations.INDEL,'):
                indel_items = line.strip().split(',')
            elif line.startswith('Locations.SNP,'):
                snp_items = line.strip().split(',')
            else:
                continue

        indexes = get_column_indexes(column_name_items)

        indel = get_csv_data(indexes, indel_items, case_name, num_bases)
        snp = get_csv_data(indexes, snp_items, case_name, num_bases)

    return indel, snp


def create_base_num_dicts(bases_file):
    """
    This function creates two dictionaries from the number of bases file, which has three columns.  The first column is
    case name, the second column is the number of bases in the whole exome, and the third column is the number bases in
    the coding exon.

    :param bases_file: The file name
    :type bases_file: str

    :return: Two dictionaries: one has key = case and value = number of bases in the whole exome; the other has
    key = case and value = number of bases in coding exon
    :rtype: tuple
    """
    num_bases_whole_exome_dict = dict()
    num_bases_coding_exons_dict = dict()

    with open(bases_file, 'r') as num_bases_file_obj:
        for line in num_bases_file_obj:
            if not line.startswith('#'):
                line_items = line.strip().split('\t')
                case = line_items[0]
                whole_exome = line_items[1]
                coding_exon = line_items[2]

                num_bases_whole_exome_dict[case] = int(whole_exome)
                num_bases_coding_exons_dict[case] = int(coding_exon)
            else:
                continue

    return num_bases_whole_exome_dict, num_bases_coding_exons_dict


def verify_required_files_exists(files_list, case):
    required = ['WholeExomeRegions.extended.csv', 'WholeExomeRegions_indelSizeDistribution.txt',
                'CodingExons.extended.csv', 'CodingExons_indelSizeDistribution.txt']

    required_set = set(required)

    suffixes_set = set()

    for filename in files_list:
        if filename.endswith('WholeExomeRegions.extended.csv'):
            suffixes_set.add('WholeExomeRegions.extended.csv')
        elif filename.endswith('WholeExomeRegions_indelSizeDistribution.txt'):
            suffixes_set.add('WholeExomeRegions_indelSizeDistribution.txt')
        elif filename.endswith('CodingExons.extended.csv'):
            suffixes_set.add('CodingExons.extended.csv')
        elif filename.endswith('CodingExons_indelSizeDistribution.txt'):
            suffixes_set.add('CodingExons_indelSizeDistribution.txt')
        else:
            continue

    missing = required_set - suffixes_set

    if missing:
        print 'Error:'
        if len(missing) > 1:
            print 'The following files are missing in {}:'.format(case)
        else:
            print 'The following file is missing in {}:'.format(case)

        for item in missing:
            print '\t', item

        sys.exit(1)


def create_output(whole_exome_indel_list, whole_exome_snp_list, coding_exons_indel_list, coding_exons_snp_list):
    """
    This function creates the output file.

    :param whole_exome_indel_list: A list the values for a single output line.
    :type whole_exome_indel_list: list
    :param whole_exome_snp_list: A list the values for a single output line.
    :type whole_exome_snp_list: list
    :param coding_exons_indel_list: A list the values for a single output line.
    :type coding_exons_indel_list: list
    :param coding_exons_snp_list: A list the values for a single output line.
    :type coding_exons_snp_list: list

    :rtype: void
    """
    with open(os.path.join(output_dir, output_file), 'w') as outfile_obj:
        header_columns = ['Case', 'Number of bases', 'Truth total', 'TP', 'FP', 'FN',
                          'TN = TotalBases - (TP + FN + FP)', 'TotalNegative  = TN + FP', 'NPA = TN/(Total Negative)',
                          'Precision', 'Recall']

        # SNPs Whole Exome
        outfile_obj.write('\tBenchmarking SNPs Whole Exome\n')
        outfile_obj.write('\t'.join(header_columns) + '\n')

        for snp in whole_exome_snp_list:
            outfile_obj.write('\t'.join(snp) + '\n')

        # INDELs Whole Exome
        outfile_obj.write('\tBenchmarking INDELs Whole Exome\n')
        outfile_obj.write('\t'.join(header_columns) + '\n')

        for indels in whole_exome_indel_list:
            for indel in indels:
                outfile_obj.write('\t'.join(indel) + '\n')

        # SNPs  Coding Exons
        outfile_obj.write('\tBenchmarking SNPs Coding Exons\n')
        outfile_obj.write('\t'.join(header_columns) + '\n')

        for snp in coding_exons_snp_list:
            outfile_obj.write('\t'.join(snp) + '\n')

        # INDELs Coding Exons
        outfile_obj.write('\tBenchmarking INDELs Coding Exons\n')
        outfile_obj.write('\t'.join(header_columns) + '\n')

        for indels in coding_exons_indel_list:
            for indel in indels:
                outfile_obj.write('\t'.join(indel) + '\n')

        # screen output
        print 'Output file created. It can be found at', os.path.join(output_dir, output_file)


########################################################################################################################
#
#   MAIN
#
########################################################################################################################
def main():
    """
    This is the main function. It iterates over all the the subdirectories in the input directory with names starting
    with 'NA' or 'HuRef'.  All CSV and TXT files with 'WholeExome' or 'CodingExons' in the file name are processed.

    A benchmarking metrics file is created.

    :rtype: void
    """
    whole_exome_indel_list = list()
    whole_exome_snp_list = list()
    coding_exons_indel_list = list()
    coding_exons_snp_list = list()
    num_bases_whole_exome_dict = dict()
    num_bases_coding_exons_dict = dict()

    prefixes = [filename[:15] for filename in os.listdir(input_dir)]

    if 'number_of_bases' not in prefixes:
        print 'Error:'
        print 'The number of bases file is missing.'
        sys.exit(1)

    for filename in os.listdir(input_dir):
        if filename.startswith('number_of_bases') and filename.endswith('txt'):
            num_bases_file = os.path.join(input_dir, filename)
            num_bases_whole_exome_dict, num_bases_coding_exons_dict = create_base_num_dicts(num_bases_file)
        else:
            continue

    for case_dir_file in os.listdir(input_dir):
        if case_dir_file.startswith('NA') or case_dir_file.startswith('HuRef'):
            files_list = os.listdir(os.path.join(input_dir, case_dir_file))

            case_name = case_dir_file

            verify_required_files_exists(files_list, case_name)

            whole_exome_indel = list()
            coding_indel = list()

            for casefilename in os.listdir(os.path.join(input_dir, case_dir_file)):
                path = os.path.join(input_dir, case_dir_file)

                if casefilename.endswith('WholeExomeRegions.extended.csv'):
                    num_bases = num_bases_whole_exome_dict[case_name]
                    whole_exome_indel, snp = get_indel_and_snp(path, casefilename, case_name, num_bases)
                    whole_exome_snp_list.append(snp)

                elif casefilename.endswith('CodingExons.extended.csv'):
                    num_bases = num_bases_coding_exons_dict[case_name]
                    coding_indel, snp = get_indel_and_snp(path, casefilename, case_name, num_bases)
                    coding_exons_snp_list.append(snp)

                else:  # skip all other files
                    continue

            for casefilename in os.listdir(os.path.join(input_dir, case_dir_file)):
                path = os.path.join(input_dir, case_dir_file)

                if casefilename.endswith('WholeExomeRegions_indelSizeDistribution.txt'):
                    num_bases = num_bases_whole_exome_dict[case_name]
                    indel_by_size = get_indel_by_size(path, casefilename, case_name, num_bases)
                    whole_exome_indel_list.append([whole_exome_indel] + indel_by_size)

                elif casefilename.endswith('CodingExons_indelSizeDistribution.txt'):
                    num_bases = num_bases_coding_exons_dict[case_name]
                    indel_by_size = get_indel_by_size(path, casefilename, case_name, num_bases)
                    coding_exons_indel_list.append([coding_indel] + indel_by_size)

                else:  # skip all other files
                    continue
        else:
            continue

    create_output(whole_exome_indel_list, whole_exome_snp_list, coding_exons_indel_list, coding_exons_snp_list)


main()
