import pandas


#######################################################################################################
def findupstreamtrimers(supporting_reads_pointmutations: pandas.DataFrame.astype,
                        counter: int) -> pandas.DataFrame.astype:
    """

    :type supporting_reads_pointmutations: pandas.DataFrame.astype
    :rtype: pandas.DataFrame.astype
    :param supporting_reads_pointmutations:
    :param counter:
    :return: pandas.DataFrame.astype
    """
    supporting_reads_pointmutations["3bpdownstream"], supporting_reads_pointmutations[
        "downstreamtrimerindex"] = [
        [supportingread["sequence"][supportingread["variant_index"]:supportingread["variant_index"] + 3]
         for rownum, supportingread in supporting_reads_pointmutations.iterrows()], counter]
    assert isinstance(supporting_reads_pointmutations, pandas.DataFrame)
    return supporting_reads_pointmutations["3bpdownstream"], supporting_reads_pointmutations["downstreamtrimerindex"]


#######################################################################################################
def finddownstreamtrimers(supporting_reads_pointmutations: pandas.DataFrame,
                          counter: int) -> tuple:
    """

    :type supporting_reads_pointmutations: pandas.DataFrame.astype
    :param supporting_reads_pointmutations:
    :param counter:
    :return:
    """
    supporting_reads_pointmutations["3bpupstream"], supporting_reads_pointmutations[
        "upstreamtrimerindex"] = [[supportingread["sequence"][
                                   supportingread["variant_index"] - 2:supportingread[
                                                                           "variant_index"] + 1] for
                                   rownum, supportingread in
                                   supporting_reads_pointmutations.iterrows()], counter]
    assert isinstance(supporting_reads_pointmutations, pandas.DataFrame)
    return supporting_reads_pointmutations["3bpupstream"], supporting_reads_pointmutations["upstreamtrimerindex"]
