import pandas


#######################################################################################################
# function that finds all supporting reads at that position
def findsupportingreads(reads: pandas.DataFrame.astype, readpos: int) -> pandas.DataFrame.astype:
    """

    :param reads:
    :param readpos:
    :return: pandas.DataFrame.astype
    :rtype: pandas.DataFrame.astype
    """
    # for each read
    # find all supporting reads spanning the start of the read (Readpos)

    # supporting reads is defined as reads going from reads["position"] till the length of the longest read
    # it is a panda DataFrame containing all the reads
    supporting_reads_pointmutations: pandas.DataFrame.astype = reads[
        (readpos >= reads["position"]) &
        (readpos < (reads["position"] + reads["sequence"].str.len()))].copy()
    assert isinstance(supporting_reads_pointmutations, pandas.DataFrame)
    return supporting_reads_pointmutations
