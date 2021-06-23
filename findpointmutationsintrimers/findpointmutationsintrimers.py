import pandas
import csv
import findreadinreference.findreadinreference as frr


#######################################################################################################
def findpointmutationsinupstreamtrimers(supporting_reads_pointmutations: pandas.DataFrame.astype,
                                        counter: int,
                                        readpos: int, reference: str, reference_dict: dict,
                                        variantpositions: dict, writer: csv) -> (dict, csv):
    """

    :type readpos: int
    :type counter: int
    :type supporting_reads_pointmutations: pandas.Dataframe
    :param writer:
    :type variantpositions: dict
    :param variantpositions:
    :param reference_dict: dict
    :param supporting_reads_pointmutations:
    :param counter:
    :param readpos:
    :param reference:
    :rtype: object
    """
    # Check ref and alt alleles for trimers
    if counter >= readpos + 3:
        supporting_reads_pointmutations["upstreamtrimerexists"] = supporting_reads_pointmutations[
            ["3bpupstream", "upstreamtrimerindex"]].apply(
            lambda row: frr.searchtrimerinreference(row["3bpupstream"], row["upstreamtrimerindex"], reference_dict),
            axis=1)
        trimercounts = supporting_reads_pointmutations["3bpupstream"].value_counts()
        uniquetrimers = supporting_reads_pointmutations["3bpupstream"].unique()

        # Get the trimer from counter-2 up to counter
        referencetrimer: str = reference[counter - 2:counter + 1]
        for uniq in uniquetrimers:
            if (
                    (referencetrimer.strip() != uniq.strip())
                    &
                    (referencetrimer[0:1].strip() != uniq[0:1].strip())
                    &
                    (referencetrimer[1:2].strip() != uniq[1:2].strip())
                    &
                    (referencetrimer[2:3].strip() != uniq[2:3].strip())
            ):

                # Noise is defined as any point mutation with less than 5 supporting reads
                if trimercounts[uniq] >= 5:
                    if counter not in variantpositions.keys():
                        variantpositions[counter] = trimercounts[uniq]
                        writer.writerow(
                            [counter, reference[counter - 2:counter + 1], uniq, trimercounts[uniq]])
                        print([counter, reference[counter - 2:counter + 1], uniq, trimercounts[uniq]])

    return variantpositions, writer


#######################################################################################################
def findpointmutationsindownstreamtrimers(supporting_reads_pointmutations: pandas.DataFrame, counter: int,
                                          readpos: int, reference: str, lensmallestread: int,
                                          reference_dict: dict,
                                          variantpositions: dict, writer: csv) -> (dict, csv):
    """

    :type lensmallestread: int
    :type counter: int
    :type supporting_reads_pointmutations: pandas.Dataframe
    :type writer: object
    :param writer:
    :type variantpositions: object
    :param variantpositions:
    :type reference_dict: object
    :param reference_dict:
    :rtype: object
    :param supporting_reads_pointmutations:
    :param counter:
    :param readpos:
    :param reference:
    :param lensmallestread:
    """
    if counter <= ((readpos + lensmallestread) + 3):
        # print 2 columns from the downstream trimer df
        # To find if the reads supporting the trimer is greater than 5
        supporting_reads_pointmutations["downstreamtrimerexists"] = supporting_reads_pointmutations[
            ["3bpdownstream", "downstreamtrimerindex"]].apply(
            lambda row: frr.searchtrimerinreference(row["3bpdownstream"],
                                                    row["downstreamtrimerindex"], reference_dict),
            axis=1)

        trimercounts = supporting_reads_pointmutations["3bpdownstream"].value_counts()
        uniquetrimers = supporting_reads_pointmutations["3bpdownstream"].unique()
        referencetrimer = reference[counter:counter + 3]
        for uniq in uniquetrimers:
            if (
                    (referencetrimer.strip() != uniq.strip())
                    &
                    (referencetrimer[0:1].strip() != uniq[0:1].strip())
                    &
                    (referencetrimer[1:2].strip() != uniq[1:2].strip())
                    &
                    (referencetrimer[2:3].strip() != uniq[2:3].strip())
            ):
                if trimercounts[uniq] >= 5:
                    variantpositions[counter] = trimercounts[uniq]
                    writer.writerow([counter, reference[counter:counter + 3], uniq, trimercounts[uniq]])
                    print([counter, reference[counter:counter + 3], uniq, trimercounts[uniq]])

        return variantpositions, writer
