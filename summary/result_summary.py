import pandas as pd


def dataframe_summarizing(open_csv_file_name, output_csv_save_name, exclude_id):
    # reading csv and creating dataframe
    df = pd.read_csv('results/' + open_csv_file_name + '.csv')

    # slicing data from the dataframe
    df_cols = df.iloc[:, 3:]
    df_initial_col = df.iloc[:, 1:3]

    # filter colums that contains even one value grater than zero

    df_remove_zero = df_cols.loc[:, [(df_cols[col] > 0).any() for col in df_cols.columns]]

    # recombine sliced dataframe row wise
    df_filter_cols = pd.concat([df_initial_col.reset_index(drop=True), df_remove_zero], axis=1)

    col_list = df_filter_cols.columns

    # create empty dataframe to get filtered row details
    result = pd.DataFrame(columns=col_list)
    print(result)

    # filter rows which contain values that don't contain tuples that have zero values in entire row
    for col_names in df_filter_cols.columns:
        if col_names != 'Reference_sequence_segment' and col_names != 'case details':
            df_filtered_row = df_filter_cols[df_filter_cols[col_names] >= 1]
            df_filtered_row = df_filtered_row[["Reference_sequence_segment", "case details", col_names]]

            # this is the final result dataframe
            result = result.append(df_filtered_row)
            print(col_names)

    print(result)
    # saving detailed summary to the location
    # result = result.sort_index(axis = 0)
    result.to_csv('results/summary/summary_detailed/' + output_csv_save_name + '.csv', index=True)

    # filtering data to get short summary
    summary = pd.DataFrame(columns=["Reference_sequence_segment", "case details", "protein_id"])
    for col_names in result.columns:
        # excepting columns that contain reference sequence segment and case details to filter nonzero records
        if col_names != 'Reference_sequence_segment' and col_names != 'case details':
            # choosing rows that contain matching frequency grater than 1
            df_filtered_row_sumery = result[result[col_names] >= 1]

            extract = df_filtered_row_sumery[["Reference_sequence_segment", "case details"]]
            # inserting protein ids and values to the filtering dataframe
            extract["protein_id"] = col_names
            extract["matching_frequency"] = df_filtered_row_sumery.loc[:, col_names]
            # appending extract to summary dataframe
            print(col_names)
            summary = summary.append(extract)

    # summary = summary.sort_index(axis=0)
    print(summary)
    # summary = summary[summary.protein_id != "sp|P05067|A4_HUMAN"]
    summary = summary[summary.protein_id != exclude_id]
    # saving short summary to the location
    summary.to_csv('results/summary/short_summary/' + output_csv_save_name + '.csv', index=True)
