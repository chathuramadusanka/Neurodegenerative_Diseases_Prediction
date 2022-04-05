from models.extract import extract_fasta_data
from summary.summary import generate_summary
from summary.result_summary import dataframe_summarizing
from analysis.analysis import analysis_common

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # extract_fasta_data()

    # keywords      for human proteome database table ==> "humanproteome"
    #               for amoloid beta database table ==> "abeta"
    #               for neurodegenerative disease table  ==> "neurodegenerative"
    # generate_summary(5, 0, "humaproteome", "abeta", "human_proteins_with_abeta")
    # generate_summary(5, 743, "humaproteome", "human_tau", "human_proteins_with_humantau")
    # generate_summary(5, 0, "neurodegenerative", "abeta", "neuro_with_abeta")
    # generate_summary(5, 743, "neurodegenerative", "human_tau", "human_proteins_with_humantau")


    # exclude_id sp|P05067|A4_HUMAN for human tau
    # exclude_id sp|P10636|TAU_HUMAN for amyloid beta
    # dataframe_summarizing('human_proteins_with_abeta', 'human_proteins_with_abeta_sumary', "sp|P05067|A4_HUMAN")
    # dataframe_summarizing('human_proteins_with_humantau', 'human_proteins_with_humantau_sumary', "sp|P10636|TAU_HUMAN")
    # dataframe_summarizing('neuro_with_abeta', 'neuro_with_abeta_sumary', "sp|P05067|A4_HUMAN")
    # dataframe_summarizing('neuro_with_humantau', 'neuro_with_humantau_sumary', "sp|P10636|TAU_HUMAN")

    analysis_common("human_proteins_with_abeta_sumary", "human_proteins_with_humantau_sumary", "neuro_with_abeta_sumary", "neuro_with_humantau_sumary")
