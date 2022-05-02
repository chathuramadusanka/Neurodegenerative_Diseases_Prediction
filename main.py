from models.extract import extract_fasta_data
from summary.summary import generate_summary
from summary.result_summary import dataframe_summarizing
from analysis.analysis import analysis_common

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # ---------------------------------------------------------------------------------------------------------
    # This method is used to extract sequences from fasta files
    # extract_fasta_data()

    # ----------------------------------------------------------------------------------------------------------
    """This section is used to find the proteins of human proteome or common nurodegenerative diseses
        that are similar to sequence segments of a reference protein"""

    # keywords      for human proteome database table ==> "humanproteome"
    #               for amyloid beta database table ==> "abeta"
    #               for neurodegenerative disease table  ==> "neurodegenerative"

    # generate_summary(5, 0, "humaproteome", "abeta", "human_proteins_with_abeta")
    # generate_summary(5, 743, "humaproteome", "human_tau", "human_proteins_with_humantau")
    # generate_summary(5, 0, "neurodegenerative", "abeta", "neuro_with_abeta")
    # generate_summary(5, 743, "neurodegenerative", "human_tau", "human_proteins_with_humantau")

    # ------------------------------------------------------------------------------------------------------------
    # in this section the results of generated_summary method is used and summarize the data removing unwanted details
    # from dataframes
    """Since we used the amyloid beta and human tau sequences as the reference we have to exclude their protein
        Id from the results"""
    # exclude_id sp|P05067|A4_HUMAN for human tau
    # exclude_id sp|P10636|TAU_HUMAN for amyloid beta

    # dataframe_summarizing('path of the raw data file', 'path to save new summary file', 'exlude_id') should be given
    # dataframe_summarizing('human_proteins_with_abeta', 'human_proteins_with_abeta_sumary', "sp|P05067|A4_HUMAN")
    # dataframe_summarizing('human_proteins_with_humantau', 'human_proteins_with_humantau_sumary', "sp|P10636|TAU_HUMAN")
    # dataframe_summarizing('neuro_with_abeta', 'neuro_with_abeta_sumary', "sp|P05067|A4_HUMAN")
    # dataframe_summarizing('neuro_with_humantau', 'neuro_with_humantau_sumary', "sp|P10636|TAU_HUMAN")

    # -------------------------------------------------------------------------------------------------------------
    # this section contain calls the analysis_common function of summarized data generated in previous section

    # analysis_common('1st summary file name', '2nd summary file name' ,'3rd file summary name, 4th file name)
    analysis_common("human_proteins_with_abeta_sumary", "human_proteins_with_humantau_sumary",
                    "neuro_with_abeta_sumary", "neuro_with_humantau_sumary")
