# Author: Lalitha Viswanathan
# Package returns position along bunch of supporting reads, upstream trimers, downstream trimers
# Naive ECGR Mutation detector (May 2021)
import pandas


#######################################################################################################
def findupstreamtrimers(supporting_reads_tofindpointmutations: pandas.DataFrame.astype,
                        counter: int) -> pandas.DataFrame.astype:
    """

    :type supporting_reads_tofindpointmutations: pandas.DataFrame.astype
    :rtype: pandas.DataFrame.astype
    :param supporting_reads_tofindpointmutations:
    :param counter:
    :return: pandas.DataFrame.astype
    """
    supporting_reads_tofindpointmutations["3bpdownstream"], supporting_reads_tofindpointmutations[
        "downstreamtrimerindex"] = [
        [supportingread["sequence"][supportingread["variant_index"]:supportingread["variant_index"] + 3]
         for rownum, supportingread in supporting_reads_tofindpointmutations.iterrows()], counter]
    assert isinstance(supporting_reads_tofindpointmutations, pandas.DataFrame)
    return supporting_reads_tofindpointmutations["3bpdownstream"], supporting_reads_tofindpointmutations["downstreamtrimerindex"]


#######################################################################################################
def finddownstreamtrimers(supporting_reads_tofindpointmutations: pandas.DataFrame,
                          counter: int) -> tuple:
    """

    :type supporting_reads_tofindpointmutations: pandas.DataFrame.astype
    :param supporting_reads_tofindpointmutations:
    :param counter:
    :return:
    """
    supporting_reads_tofindpointmutations["3bpupstream"], supporting_reads_tofindpointmutations[
        "upstreamtrimerindex"] = [[supportingread["sequence"][
                                   supportingread["variant_index"] - 2:supportingread[
                                                                           "variant_index"] + 1] for
                                   rownum, supportingread in
                                   supporting_reads_tofindpointmutations.iterrows()], counter]
    assert isinstance(supporting_reads_tofindpointmutations, pandas.DataFrame)
    return supporting_reads_tofindpointmutations["3bpupstream"], supporting_reads_tofindpointmutations["upstreamtrimerindex"]
