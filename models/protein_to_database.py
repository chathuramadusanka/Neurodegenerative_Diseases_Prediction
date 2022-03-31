import sqlite3


def improt_sequence(fasta_seq, fasta_id, table_name):
    sequnce_dict = {}
    sequence_list = []
    for counter_id, counter_seq in zip(fasta_id, fasta_seq):
        # print(counter_id, '-->', table_name, '-->', counter_seq)
        sequnce_dict = {"protein Id": counter_id, "sequence": str(counter_seq)}

        # print(str(counter_seq))
        sequence_list.append(sequnce_dict)

    save_proteins_to_database(sequence_list, table_name)
    # print(sequence_list)
    return 0


def save_proteins_to_database(seq_contain_list, table_name):
    conn = sqlite3.connect('database.db')

    # Create a cursor to allow executing SQL commands
    cursor = conn.cursor()

    sql_command = "CREATE TABLE IF NOT EXISTS " + table_name + "( Id INTEGER PRIMARY KEY AUTOINCREMENT, Sequenceid " \
                                                               "TEXT, Sequence TEXT) "

    cursor.execute(sql_command)
    for protein in seq_contain_list:
        quary_val_1 = str(protein['protein Id'])
        quary_val_2 = str(protein['sequence'])

        quary_val_all = "VALUES ('" + quary_val_1 + "','" + quary_val_2 + "')"

        # print(quary_val_all)

        insert_data = f"INSERT INTO " + table_name + " (Sequenceid, sequence) " + quary_val_all
        print(insert_data)
        cursor.execute(insert_data)
        conn.commit()
    # Commit the changes to the database
    conn.commit()
