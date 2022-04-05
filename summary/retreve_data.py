from os import error
from re import fullmatch
import sqlite3
from sqlite3.dbapi2 import Cursor, Error


# defining starting and ending connection========================================
def start_connection():
    conn = None
    try:
        conn = sqlite3.connect("database.db")
        message = "Connection Success"
    except Error as e:
        message = e

    return conn, message


def close_connection(conn):
    conn.close()
    print("connection close")


# ===================================================================================


# defining method to get the data of a certain table================================
def data_retrive(table_name):
    conn, notification = start_connection()
    print(notification)

    conn.row_factory = sqlite3.Row
    # creating cursor object
    cursor = conn.cursor()
    sequence_list = []
    sequence_list_name = []
    full_dict = {}
    full_list = []
    # quary_search = 'SELECT * FROM "{}"'.fromat(table_name)

    try:
        # take all the meta tags with all the datatypes and sizes that defined
        # cursor.execute("PRAGMA table_info(protein_seq_1)")

        # select all the records in the database rocord where some colum quals to something using variable methods
        # cursor.execute('SELECT Sequenceid, Sequence FROM {}'.format(table_name))

        cursor.execute('SELECT * FROM {}'.format(table_name))
        records = cursor.fetchall()

        # get the all sequences from 
        for list in records:
            # for listin compare sequences data tables
            # for lists in list[2:4]:
            sequence_list_name.append(list[1])
            sequence_list.append(list[2])

        for pro_id, sequence in zip(sequence_list_name, sequence_list):
            full_dict = {"protein_id": pro_id, "sequence": str(sequence)}

            full_list.append(full_dict)

        full_dict = dict(zip(sequence_list_name, sequence_list))
        # print(full_dict)

        # print((sequence_list))
        # print(len(col_name_list))

    except Error as re_er:
        message = re_er

    finally:
        if conn:
            close_connection(conn)

    print(full_dict)
    return sequence_list, full_dict
