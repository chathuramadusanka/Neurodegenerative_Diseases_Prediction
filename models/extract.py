from Bio.PDB import *
from Bio import SeqIO
from models.protein_to_database import improt_sequence


# defining variables to load sequence of fasta files
def sequence_extract_fasta(fasta_files):
    # Defining empty list for the Fasta id and fasta sequence variables
    fasta_id = []
    fasta_seq = []

    # opening given fasta file using the file path
    with open(fasta_files, 'r') as fasta_file:
        # extracting multiple data in single fasta file using biopython
        for record in SeqIO.parse(fasta_file, 'fasta'):  # (file handle, file format)
            # print("\n",'>' + record.id)
            # print(record.seq)

            # appending extracted fasta data to empty lists variables
            fasta_seq.append(record.seq)
            fasta_id.append(record.id)
            fasta_name = record.name

    # returning fasta_id and fasta sequence to both call_compare_fasta and call_reference_fasta
    return fasta_id, fasta_seq


# defining method for call reference fasta file
def call_reference_fasta(reference_path):
    reference_fasta_id, reference_fasta_sequence = sequence_extract_fasta(reference_path)

    # returning reference fasta id and fasta sequence to extract_fasta_data() method
    return reference_fasta_id, reference_fasta_sequence


# defining method for call compare fasta file

def call_compare_fasta(compare_file_path):
    compare_fasta_id, compare_fasta_seq = sequence_extract_fasta(compare_file_path)
    # returning compare fasta id and fasta sequence to extract_fasta_data() method
    return compare_fasta_id, compare_fasta_seq


# this method containsing the commands to give to save fasta data to database
def extract_fasta_data():
    # ===========================================================================================================
    # giving path to the reference fasta file to extract data

    # reference_path = 'data_set/rcsb_pdb_1IYT.fasta'  # for amyloid beta to reference table
    reference_path = 'data_set/human_tau.fasta'    # for human tau to reference table

    # calling call_reference_fasta method to extract reference fasta data
    fasta_id_ref, fasta_seq_ref = call_reference_fasta(reference_path)

    # table name variable for insert to database as reference sequences
    # table_name = "abeta"  # for amyloid beta
    table_name = "human_tau"              # for human tau

    # ============================================================================================================
    # giving path to the compare fasta file to extract data
    # sequence and ids are output

    # compare_fasta = 'data_set/uniprot-proteome_UP000005640+reviewed_yes.fasta'             # for human proteome 21000
    compare_fasta = 'data_set/uniprot-neurodegenerative+disease-filtered-reviewed_yes.fasta'  # for
    # neurodegenerative diseses data neurodegenarative diseases compare_fasta = 'data_set/test_compare.fasta' for
    # testing fasta sequences

    fasta_id_comp, fasta_seq_comp = call_compare_fasta(compare_fasta)
    # print("compare >", (fasta_seq_comp[0]))

    # table name variable for insert to database as compare sequences
    # table_name = "humaproteome"                    #with respect to human proteome
    # table_name = "neurodegenerative"                 # with respect to neurodegenerative disease dataset

    # ===============================================================================================================

    # for saving reference sequences for the database
    table_details_all = improt_sequence(fasta_seq_ref, fasta_id_ref, table_name)

    # for saving compare sequences for the database
    # table_details_all = improt_sequence(fasta_seq_comp, fasta_id_comp, table_name)

    # print(table_details_all)
