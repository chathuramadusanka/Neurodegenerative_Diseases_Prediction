from sqlite3.dbapi2 import Row
from summary.compare import refence_seq_pattern
import pandas as pd
import csv


def generate_summary(lowest_num_of_residues, starting_case, database_tabel_name, reference_table_name,
                     result_saving_file_name):
    # from plot import extract_csv_gen_plot

    # calling method to get the summary details of the protein sequence matching
    # the fist parameter is to indicate the lowest needed character combination during the sequence analysis
    # the second parameter is the database table name where the protein sequences that need to compare
    # for amyloid beta description is ""original_seq"" # for human tau "human_tau"
    # all_count_details_dictionary = refence_seq_pattern(6, 750, "humaproteome", "human_tau")

    all_count_details_dictionary = refence_seq_pattern(lowest_num_of_residues, starting_case, database_tabel_name,
                                                       reference_table_name)

    # initializing csv file path
    csv_path = 'results/' + result_saving_file_name + '.csv'

    # creating dataframe using pandas jason normalizer and passing variables; used to extract data from nested
    # dictionary df = pd.json_normalize(all_count_details_dictionary, record_path=None)

    # creating pandas dataframe using concat mehtod to extract data from dictionary
    df = pd.concat([pd.DataFrame(l) for l in all_count_details_dictionary], axis=1).T

    # saving the dataframe to the csv file
    df.to_csv(csv_path, index=True)
    print("successfully saved")
