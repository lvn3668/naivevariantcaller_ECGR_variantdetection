# Author: Lalitha Viswanathan
# Naive ECGR mutation detector (May 2021)
# Package builds all dimers upstream of a given position along bunch of supporting reads
# Package builds all dimers downstream of a given position along bunch of supporting reads
# Returns 4 artefacts (the position along the read, all dimers upstream, all dimers downstream)
import pandas


#######################################################################################################
def findupstreamdimers(supporting_reads_tofindpointmutations: pandas.DataFrame.astype,
                       counter: int) -> tuple:
    """

    :type supporting_reads_tofindpointmutations: pandas.DataFrame.astype
    :param supporting_reads_tofindpointmutations:
    :param counter:
    :return: pandas.DataFrame.astype
    :rtype: pandas.DataFrame.astype
    """
    supporting_reads_tofindpointmutations['2bpupstream'], supporting_reads_tofindpointmutations[
        'upstreamdimerindex'] = [[supportingread['sequence'][
                                  supportingread['variant_index'] - 1:supportingread[
                                                                          'variant_index'] + 1] for
                                  rownum, supportingread in supporting_reads_tofindpointmutations.iterrows()],
                                 counter]
    assert isinstance(supporting_reads_tofindpointmutations, pandas.DataFrame)
    return supporting_reads_tofindpointmutations['2bpupstream'], supporting_reads_tofindpointmutations['upstreamdimerindex']


#######################################################################################################
def finddownstreamdimers(supporting_reads_pointmutations: pandas.DataFrame,
                         counter: int) -> tuple:
    """

    :type counter: int
    :type supporting_reads_pointmutations: pandas.DataFrame
    :rtype: pandas.Dataframe
    :param supporting_reads_pointmutations:
    :param counter:
    :return:
    """
    rownum: int
    supportingread: tuple
    supporting_reads_pointmutations['2bpdownstream'], supporting_reads_pointmutations[
        'downstreamdimerindex'] = [[supportingread["sequence"][
                                    supportingread["variant_index"]:supportingread[
                                                                        "variant_index"] + 2] for
                                    rownum, supportingread in supporting_reads_pointmutations.iterrows()],
                                   counter]
    assert isinstance(supporting_reads_pointmutations, pandas.DataFrame)
    return supporting_reads_pointmutations['2bpdownstream'], supporting_reads_pointmutations['downstreamdimerindex']
