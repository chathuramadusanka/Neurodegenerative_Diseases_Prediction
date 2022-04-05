import pandas as pd


def analysis_common(file_1_name, file_2_name, file_3_name, file_4_name):
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

    compare_1 = set(df_1['protein_id']).intersection(df_2['protein_id'])  # abeta ==> proteome + tau ==> proteome
    compare_2 = set(df_3['protein_id']).intersection(df_1['protein_id'])  # abeta ==> proteome + abeta ==> neuro
    compare_3 = set(df_2['protein_id']).intersection(df_4['protein_id'])  # tau ==> proteome + tau ==> neuro
    compare_4 = set(df_3['protein_id']).intersection(df_4['protein_id'])  # abeta ==> neuro + tau ==> neuro
    compare_5 = set(df_1['protein_id']).intersection(df_4['protein_id'])  # abeta ==> proteome + tau ==> neuro
    compare_6 = set(df_2['protein_id']).intersection(df_3['protein_id'])  # tau ==> proteome + abeta ==> neuro

    print(compare_1)

    data_1 = pd.concat([df_1.query('protein_id in @compare_1'), df_2.query('protein_id in @compare_1')],
                       ignore_index=True)
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

    unique_dogs = data_7.drop_duplicates(
        subset=["Unnamed: 0", "Reference_sequence_segment", "case details", "protein_id", "matching_frequency"])

    print(data_1)
