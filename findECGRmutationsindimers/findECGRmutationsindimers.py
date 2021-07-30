# Author: Lalitha Viswanathan
# Naive ECGR Mutation detector
# May 2021
# Python package checks if reference dimer differs from dimer in supporting reads at a given position
# For both upstream and downstream dimers
from typing import Any, Union

import pandas
import csv

from pandas import Series

from findreadinreference.findreadinreference import searchtrimerinreference


#######################################################################################################
def findECGRmutationsinupstreamdimers(supporting_reads_pointmutations: pandas.DataFrame.astype, counter: int,
                                      readpos: int, reference: str, reference_dict: dict,
                                      variantpositions: dict, writer: csv) -> (dict, csv):
    """

    :rtype: object
    :param writer:
    :type variantpositions: object
    :type reference_dict: object
    :param reference_dict:
    :param supporting_reads_pointmutations:
    :param counter:
    :param readpos:
    :param reference:
    """
    if counter >= readpos + 2:
        supporting_reads_pointmutations["upstreamdimerexists"] = supporting_reads_pointmutations[
            ["2bpupstream", "upstreamdimerindex"]].apply(
            lambda row: searchtrimerinreference(row["2bpupstream"], row["upstreamdimerindex"], reference_dict), axis=1)
        dimercounts = supporting_reads_pointmutations["2bpupstream"].value_counts()
        uniquedimers = supporting_reads_pointmutations["2bpupstream"].unique()
        referencedimer = reference[counter - 1:counter + 1]
        for uniq in uniquedimers:
            if (
                    ((referencedimer.strip()) != uniq.strip())
                    &
                    (referencedimer[0:1].strip() != uniq[0:1].strip())
                    &
                    (referencedimer[1:2].strip() != uniq[1:2].strip())
            ):
                if dimercounts[uniq] >= 5:
                    if counter not in variantpositions.keys():
                        variantpositions[counter] = dimercounts[uniq]
                        writer.writerow(
                            [counter, reference[counter - 1:counter + 1], uniq, dimercounts[uniq]])
                        print([counter, reference[counter - 1:counter + 1], uniq, dimercounts[uniq]])

    return variantpositions, writer


#######################################################################################################
def findECGRmutationsindownstreamdimers(supporting_reads_pointmutations: pandas.DataFrame.astype, counter: int,
                                        readpos: int,
                                        lensmallestread: int,
                                        reference: str, reference_dict: dict, variantpositions: dict,
                                        writer: csv) -> (dict, csv):
    """

    :type counter: int
    :param reference:
    :type lensmallestread: int
    :param lensmallestread:
    :param counter:
    :type readpos: int
    :param readpos:
    :param supporting_reads_pointmutations:
    :param writer:
    :type variantpositions: object
    :param variantpositions:
    :type reference_dict: object
    """

    if counter <= ((readpos + lensmallestread) + 2):
        # print 2 columns from the downstream trimer df
        # To find if the reads supporting the trimer is greater than 5
        assert isinstance(supporting_reads_pointmutations, pandas.DataFrame)
        supporting_reads_pointmutations["downstreamdimerexists"] = supporting_reads_pointmutations[
            ["2bpdownstream", "downstreamdimerindex"]].apply(
            lambda row: searchtrimerinreference(row["2bpdownstream"], row["downstreamdimerindex"], reference_dict),
            axis=1)
        dimercounts: Union[Series, Any] = supporting_reads_pointmutations["2bpdownstream"].value_counts()
        uniquedimers = supporting_reads_pointmutations["2bpdownstream"].unique()
        referencedimer: str = reference[counter:counter + 2]
        for uniq in uniquedimers:
            if (
                    (referencedimer.strip() != uniq.strip())
                    &
                    (referencedimer[0:1].strip() != uniq[0:1].strip())
                    &
                    (referencedimer[1:2].strip() != uniq[1:2].strip())

            ):
                if dimercounts[uniq] >= 5:
                    if counter not in variantpositions.keys():
                        variantpositions[counter] = dimercounts[uniq]
                        writer.writerow([counter, reference[counter:counter + 2], uniq, dimercounts[uniq]])
                        print([counter, reference[counter:counter + 2], uniq, dimercounts[uniq]])

        return variantpositions, writer

#######################################################################################################
