import pandas as pd
# hydrophobicity, pI value and instability index
import peptides
# from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP

# for plotting sequences using plotly
import plotly.express as px

"""# defining plotting function
def plot(table):
    plt.rcParams["figure.figsize"] = [20, 15]
    table.plot(kind="bar")
    plt.title("Total Matching Frequency vs Protein Sequence Segments")
    plt.xlabel("Protein Sequence Segment")
    plt.ylabel("Matching Frequency")
    plt.show()"""


# defining analysis function to analyse all the extracted results through the datasets
def analysis_common(file_1_name, file_2_name, file_3_name, file_4_name):
    # loading CSV files into dataframes =============================================================
    try:
        df_1 = pd.read_csv('results/summary/short_summary/' + file_1_name + '.csv')  # abeta ==> proteome
    except:
        df_1 = pd.DataFrame()

    try:
        df_2 = pd.read_csv('results/summary/short_summary/' + file_2_name + '.csv')  # tau ==> proteome
    except:
        df_2 = pd.DataFrame()

    try:
        df_3 = pd.read_csv('results/summary/short_summary/' + file_3_name + '.csv')  # abeta ==> neuro
    except:
        df_3 = pd.DataFrame()

    try:
        df_4 = pd.read_csv('results/summary/short_summary/' + file_4_name + '.csv')  # tau ==> neuro
    except:
        df_4 = pd.DataFrame()

    # ====================================================================================================

    # get all intersecting protein id s with all the imported dataframes
    compare_1 = set(df_1['protein_id']).intersection(df_2['protein_id'])  # abeta ==> proteome + tau ==> proteome
    compare_2 = set(df_3['protein_id']).intersection(df_1['protein_id'])  # abeta ==> proteome + abeta ==> neuro
    compare_3 = set(df_2['protein_id']).intersection(df_4['protein_id'])  # tau ==> proteome + tau ==> neuro
    compare_4 = set(df_3['protein_id']).intersection(df_4['protein_id'])  # abeta ==> neuro + tau ==> neuro
    compare_5 = set(df_1['protein_id']).intersection(df_4['protein_id'])  # abeta ==> proteome + tau ==> neuro
    compare_6 = set(df_2['protein_id']).intersection(df_3['protein_id'])  # tau ==> proteome + abeta ==> neuro

    # print(compare_1)

    # filtering out all the details that contain intersected protein ID's
    data_1 = pd.concat([df_1.query('protein_id in @compare_1'), df_2.query('protein_id in @compare_1')],
                       ignore_index=True)
    # removing duplicates in dataframe
    data_1_1 = data_1.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    data_2 = pd.concat([df_3.query('protein_id in @compare_2'), df_1.query('protein_id in @compare_2')],
                       ignore_index=True)
    data_1_2 = data_2.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    data_3 = pd.concat([df_2.query('protein_id in @compare_3'), df_4.query('protein_id in @compare_3')],
                       ignore_index=True)
    data_1_3 = data_3.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    data_4 = pd.concat([df_3.query('protein_id in @compare_4'), df_4.query('protein_id in @compare_4')],
                       ignore_index=True)
    data_1_4 = data_4.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    data_5 = pd.concat([df_1.query('protein_id in @compare_5'), df_4.query('protein_id in @compare_5')],
                       ignore_index=True)
    data_1_5 = data_5.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    data_6 = pd.concat([df_2.query('protein_id in @compare_6'), df_3.query('protein_id in @compare_6')],
                       ignore_index=True)
    data_1_6 = data_6.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    compare_7 = set(data_1['protein_id']).intersection(data_4['protein_id'])

    data_7 = pd.concat([data_1.query('protein_id in @compare_7'), data_2.query('protein_id in @compare_7')],
                       ignore_index=True)

    data_1_7 = data_7.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    # =============================================================================================================

    # contain all the intersections with common matching sequence segments

    compare_1_seq = set(df_1['Reference_sequence_segment']).intersection(
        df_2['Reference_sequence_segment'])  # abeta ==>
    # proteome + tau ==> proteome
    compare_2_seq = set(df_3['Reference_sequence_segment']).intersection(
        df_1['Reference_sequence_segment'])  # abeta ==>
    # proteome + abeta ==> neuro
    compare_3_seq = set(df_2['Reference_sequence_segment']).intersection(df_4['Reference_sequence_segment'])  # tau ==>
    # proteome + tau ==> neuro
    compare_4_seq = set(df_3['Reference_sequence_segment']).intersection(
        df_4['Reference_sequence_segment'])  # abeta ==>
    # neuro + tau ==> neuro
    compare_5_seq = set(df_1['Reference_sequence_segment']).intersection(
        df_4['Reference_sequence_segment'])  # abeta ==>
    # proteome + tau ==> neuro
    compare_6_seq = set(df_2['Reference_sequence_segment']).intersection(df_3['Reference_sequence_segment'])  # tau ==>
    # proteome + abeta ==> neuro

    print(compare_1)

    # filtering out all the details that contain intersected protein ID's
    data_1_seq = pd.concat([df_1.query('Reference_sequence_segment in @compare_1_seq'),
                            df_2.query('Reference_sequence_segment in @compare_1_seq')],
                           ignore_index=True)
    # removing duplicates in dataframe
    data_1_1_seq = data_1_seq.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    data_2_seq = pd.concat([df_3.query('Reference_sequence_segment in @compare_2_seq'),
                            df_1.query('Reference_sequence_segment in @compare_2_seq')],
                           ignore_index=True)
    data_1_2_seq = data_2_seq.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    data_3_seq = pd.concat([df_2.query('Reference_sequence_segment in @compare_3_seq'),
                            df_4.query('Reference_sequence_segment in @compare_3_seq')],
                           ignore_index=True)
    data_1_3_seq = data_3_seq.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    data_4_seq = pd.concat([df_3.query('Reference_sequence_segment in @compare_4_seq'),
                            df_4.query('Reference_sequence_segment in @compare_4_seq')],
                           ignore_index=True)
    data_1_4_seq = data_4_seq.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    data_5_seq = pd.concat([df_1.query('Reference_sequence_segment in @compare_5_seq'),
                            df_4.query('Reference_sequence_segment in @compare_5_seq')],
                           ignore_index=True)
    data_1_5_seq = data_5_seq.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    data_6_seq = pd.concat([df_2.query('Reference_sequence_segment in @compare_6_seq'),
                            df_3.query('Reference_sequence_segment in @compare_6_seq')],
                           ignore_index=True)
    data_1_6_seq = data_6_seq.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    compare_7_seq = set(data_1['Reference_sequence_segment']).intersection(data_4['Reference_sequence_segment'])

    data_7_seq = pd.concat([data_1.query('Reference_sequence_segment in @compare_7_seq'),
                            data_2.query('Reference_sequence_segment in @compare_7_seq')],
                           ignore_index=True)

    data_1_7_seq = data_7_seq.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    # =========================================================================================================

    # sending all the analysis data to the CSV files

    # path for the protein id intersection analysis
    csv_path_id = 'results/analysis/pro_id/'

    # path for the protein sequence segment intersection analysis
    csv_path_seq = 'results/analysis/pro_seq/'

    data_1_1.to_csv(csv_path_id + file_1_name + 'analysis_with' + file_2_name + '.csv', index=True)
    data_1_1_seq.to_csv(csv_path_seq + file_1_name + 'analysis_seq_with' + file_2_name + '.csv', index=True)

    data_1_2.to_csv(csv_path_id + file_1_name + 'analysis_with' + file_3_name + '.csv', index=True)
    data_1_2_seq.to_csv(csv_path_seq + file_1_name + 'analysis_seq_with' + file_3_name + '.csv', index=True)

    data_1_3.to_csv(csv_path_id + file_2_name + 'analysis_with' + file_4_name + '.csv', index=True)
    data_1_3_seq.to_csv(csv_path_seq + file_2_name + 'analysis_seq_with' + file_4_name + '.csv', index=True)

    data_1_4.to_csv(csv_path_id + file_3_name + 'analysis_with' + file_4_name + '.csv', index=True)
    data_1_4_seq.to_csv(csv_path_seq + file_3_name + 'analysis_seq_with' + file_4_name + '.csv', index=True)

    data_1_5.to_csv(csv_path_id + file_1_name + 'analysis_with' + file_4_name + '.csv', index=True)
    data_1_5_seq.to_csv(csv_path_seq + file_1_name + 'analysis_seq_with' + file_4_name + '.csv', index=True)

    data_1_6.to_csv(csv_path_id + file_2_name + 'analysis_with' + file_3_name + '.csv', index=True)
    data_1_6_seq.to_csv(csv_path_seq + file_2_name + 'analysis_seq_with' + file_3_name + '.csv', index=True)

    data_1_7.to_csv(csv_path_id + 'analysis_with_4_results.csv', index=True)
    data_1_7_seq.to_csv(csv_path_seq + 'analysis_seq_with_4_dataframes.csv', index=True)

    # ====================================================================================================

    # adding all the dataframes into one list
    frames = [df_1, df_2, df_3, df_4]

    # concat all the frames into one dataframe
    full_list = pd.concat(frames)

    # dropping duplicate records of dataframe
    unique_full_list = full_list.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    # print(unique_full_list)

    # get the all the matching sequences into one place( unique_full_list) with all the datasets and export to CSV file
    csv_path_all = 'results/analysis/'
    unique_full_list.to_csv(csv_path_all + 'all_in_4_matched_dataframes.csv', index=True)

    # =========================================================================================================
    # plotting frequency of matching sequences

    # dropping Unnamed:0 and "case detail" columns form dataframe
    table = unique_full_list.drop(columns=["Unnamed: 0", "case details"])

    # Calculating the total sum of each matching frequency of reference sequence segment
    table['total_matching'] = unique_full_list.groupby(['Reference_sequence_segment'])[
        'matching_frequency'].transform('sum')

    # dropping matching sequence column
    table = table.drop(columns=["matching_frequency"])

    # dropping duplicates considering both reference sequence segment and total matching sequence
    table = table.drop_duplicates(
        subset=["Reference_sequence_segment", "total_matching"])

    print(table.dtypes)
    table = table.drop(columns=["protein_id"])


    # defining method to return hydrophobic values of sequences to the dataframe
    def get_hydrophobicity(df_row):
        single_seq = df_row["Reference_sequence_segment"]
        peptide = peptides.Peptide(single_seq)
        hydrophobicity = peptide.hydrophobicity(scale="KyteDoolittle")
        return hydrophobicity

    # defining method to return pI values of sequences to the dataframe
    def get_isoelectric(df_row):
        single_seq = df_row["Reference_sequence_segment"]
        peptide = peptides.Peptide(single_seq)
        isoelectic_point = peptide.isoelectric_point(pKscale="EMBOSS")

        """# The pI is a variable that affects the solubility of the peptides under
            # certain conditions of pH. When the pH of the solvent is equal to the pI of the protein,
            # it tends to precipitate and lose its biological function."""

        return isoelectic_point

    # defining method to return instability index values of sequences to the dataframe
    def get_instability(df_row):
        single_seq = df_row["Reference_sequence_segment"]
        peptide = peptides.Peptide(single_seq)
        instability = peptide.instability_index()

        """If instability index <= 40 the protein is stable"""

        return instability

    # defining method to return length of sequences to the dataframe
    def get_length(df_row):
        single_seq = df_row["Reference_sequence_segment"]
        length = len(single_seq)
        return length

    # defining method to return charge of sequences to the dataframe
    def get_charge(df_row):
        single_seq = df_row["Reference_sequence_segment"]
        peptide = peptides.Peptide(single_seq)
        charge = peptide.charge(pKscale="Lehninger")
        return charge

    table['hydrophobicity'] = table.apply(get_hydrophobicity, axis=1)
    table["isoelectic_point"] = table.apply(get_isoelectric, axis=1)
    table["instability"] = table.apply(get_instability, axis=1)
    table["length"] = table.apply(get_length, axis=1)
    table["charge"] = table.apply(get_charge, axis=1)

    table.set_index('Reference_sequence_segment', inplace=True)

    max_length = table[table["length"] >= 13].sort_values("length")
    higest_matching = table[table["total_matching"] > 14].sort_values("total_matching")
    isoelectic_aggregates = table[table["isoelectic_point"].between(6.0, 8.8)].sort_values("isoelectic_point")
    hydrophobicity = table[table["hydrophobicity"] > 0].sort_values("hydrophobicity")
    instability_index = table[table["instability"] < 40].sort_values("instability")
    all_filter = table.query('instability < 40 & hydrophobicity > 0 & isoelectic_point.between(6.0, 8.8)').sort_values(
        "isoelectic_point")
    iso_fil2 = all_filter.drop(columns=["total_matching", "charge", "length"])
    iso_fil2.rename(columns={'hydrophobicity': 'Hydrophobicity', 'isoelectic_point': 'pI value',
                             'instability	': 'instability index '}, inplace=True)

    csv_path_stability = 'results/analysis/stability/'
    max_length.to_csv(csv_path_stability + 'max_length.csv', index=True)
    higest_matching.to_csv(csv_path_stability + 'highest_matching.csv', index=True)
    isoelectic_aggregates.to_csv(csv_path_stability + 'pi.csv', index=True)
    hydrophobicity.to_csv(csv_path_stability + 'hydrophobicity.csv', index=True)
    instability_index.to_csv(csv_path_stability + 'instability_index.csv', index=True)
    iso_fil2.to_csv(csv_path_stability + 'all_filters_at_once.csv', index=True)

