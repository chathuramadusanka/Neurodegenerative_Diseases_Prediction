import gc
import sys
import time
from typing import DefaultDict

from summary.retreve_data import data_retrive


def refence_seq_pattern(lowest_seq_combination, initial_case, table_name, reference_table):
    # process flow 1

    # retrieve reference frequency from  retreve_data.py
    sequence_list, reference_dict = data_retrive(reference_table)
    # sequence_list = ["protein_id", "123456789x"]

    # retrieve data from the compare protein sequences database tables only to pass into the other method
    seq_comp_list, seq_com_dict = data_retrive(table_name)
    seq_comp_list = seq_comp_list[0:]  # in case if you need to control the list you can change the attributes inside[]
    # print(seq_comp_list)
    # seq_com_dict = {'id_01' : '54671234512345', "id_02": "01234565789x"}

    # get only the amyloid beta sequence
    ref_seq = sequence_list[0]
    print(ref_seq)
    ref_length = len(ref_seq)  # get the length of the beta amyloid sequence
    ref_indexes = (ref_length)  # total number of indexes thinking it is a list
    # print(ref_length, ref_indexes)

    # initializing the lowest sequence combination
    # means that what is the minimum number of characters should be there in any combination
    low_com = lowest_seq_combination + 1

    # crating and initializing case variables
    inital_case = initial_case
    max_cases = ref_indexes  # total number of indexes assing to max cases
    cases_to_run = (max_cases - (low_com - 1))  # total number fo cases need to run
    sequence_grab_list = []  # all the sequence combinations will be added
    all_count_details_list = []  # all count details and collection of all the dictionaries list

    animation = ["[■□□□□□□□□□]", "[■■□□□□□□□□]", "[■■■□□□□□□□]", "[■■■■□□□□□□]", "[■■■■■□□□□□]", "[■■■■■■□□□□]",
                 "[■■■■■■■□□□]", "[■■■■■■■■□□]", "[■■■■■■■■■□]", "[■■■■■■■■■■]"]

    for cases in range(ref_length - inital_case):
        iteration_count = 0  # iteration count is initialize to 0 at outer loop
        back_iteration = inital_case  # back iteration is reusing form its initial case to  bottom words

        time.sleep(0.2)
        sys.stdout.write("\r" + animation[cases % len(animation)])
        sys.stdout.flush()

        if inital_case <= cases_to_run:
            for row in range(inital_case + 1):
                if iteration_count <= inital_case:
                    # print("case_%s" %inital_case,inital_case, iteration_count)

                    sequence_grab = ref_seq[iteration_count:(ref_indexes - back_iteration)]
                    # print("case_%s" %inital_case,inital_case, iteration_count, sequence_grab)
                    sequence_grab_list.append(sequence_grab)

                    # passing data to compare_seq method and returning summery details and getting output from
                    # process flow 2 : output is a dictionary
                    case_details = ("case_%s" % inital_case)
                    all_count_details = compare_sequences(sequence_grab, low_com, seq_comp_list, seq_com_dict,
                                                          case_details)

                    # adding this loop case details and sequence to the dictionary
                    # this addition is no longer required because it has been already addend in the previous method
                    # returned details
                    # all_count_details["Ref_seq_combination"] = ["case_%s" %inital_case, sequence_grab]

                    # adding all_count_details to the all_count_details_list when iterating
                    all_count_details_list.append(all_count_details)

                    iteration_count += 1
                    back_iteration -= 1
        inital_case += 1
        gc.collect()

    return all_count_details_list


# defining method to compare sequences
def compare_sequences(grabed_sequene, low_com, seq_comp_list, compare_seq_dict, case_details):
    # process flow 2
    # reference sequence combination that pass form reference_seq_pattern method
    grabed_sequene = grabed_sequene
    grabed_sequene_rev = grabed_sequene[::-1]
    len_grab = len(grabed_sequene)

    # this is the small part of the sequence that get form the grabbed_sequence
    sequence_part1 = grabed_sequene[0:low_com - 1]

    # getting the reverse complement of the sequence_part1
    sequence_part2 = sequence_part1[::-1]

    # this dict items contains both protein sequnce id and sequence itself
    comp_dict_items = compare_seq_dict.items()

    # initialising match frequency to 0
    match_frequency = []
    match_frequency_reverse = []
    column_name_list = []
    case_details_list = []
    grabed_sequene_list = []

    # getting extracted protein column name list as x and sequence as y  in separate variables
    for column_name, seq_comp in comp_dict_items:

        count = 0
        reversed_count = 0

        # check whether the sequence_part1 is equal to the given / selected_sequence
        if seq_comp.find(sequence_part1) != -1:
            # get the sequence match count to the specific reference sequence
            count = count_occurrences(seq_comp, grabed_sequene)

        # check whether the sequence_part2 is equal to the given / selected sequence
        if seq_comp.find(sequence_part2) != -1:
            reversed_count = count_occurrences(seq_comp, grabed_sequene_rev)
            count = reversed_count + count

        # adding each values to list throughout the iteration
        match_frequency.append(count)
        match_frequency_reverse.append(reversed_count)
        column_name_list.append(column_name)
        case_details_list.append(case_details)
        grabed_sequene_list.append(grabed_sequene)

    # all the results in one variable which is a dictionary  # this needs complex procedure when extracting results
    # last_result = dict(zip(column_name_list, zip(match_frequency, match_frequency_reverse)))

    # using this method we can simply create a variable independent dictionary
    last_result = DefaultDict(dict)
    for case, grb_seq, col_name, matc_fre, rev_fre in zip(case_details_list, grabed_sequene_list, column_name_list,
                                                          match_frequency, match_frequency_reverse):
        last_result["id"]["case details"] = case
        last_result["id"]["Reference_sequence_segment"] = grb_seq
        last_result["id"][col_name] = matc_fre
        # last_result["id"][column_name] = rev_fre

    # print(column_name_list)
    # print(seq_comp_list[0:2])
    # print(seq_comp_list)

    # print(grabed_sequene, len_grab, low_com, sequence_part1, sequence_part2)
    # print("original", match_frequency, "reverse", match_frequency_reverse)
    gc.collect()
    return last_result
    # retun to process flow 1


def count_occurrences(full_sequence, grab_sequence):
    # process flow 3

    # Initialize count and start to 0
    count = 0
    start = 0

    # Search through the full_sequence till
    # we reach the end of it
    while start < len(full_sequence):

        # Check if a grab_sequence is present from
        # 'start' position till the end
        pos = full_sequence.find(grab_sequence, start)

        if pos != -1:
            # If a grab_sequence is present, move 'start' to
            # the next position from start of the grab_sequence
            start = pos + 1

            # Increment the count
            count += 1
        else:
            # If no further grab_sequence is present
            break
    # return the value of count
    return count
    # return to process flow 2


"""all_count_details_dictionary = refence_seq_pattern(1)

print(len(all_count_details_dictionary))
for values_list in range(len(all_count_details_dictionary)):
    print(all_count_details_dictionary[values_list])"""
