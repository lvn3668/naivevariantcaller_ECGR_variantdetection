# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import sys
import time
from typing import TextIO

import referenceutilities.indexreference as refindex
import referenceutilities.readreference as refutilites
import csv
import pandas
import findsupportingreads.findsupportingreads as fsr
import findalldimers.findalldimers as dimers
import findalltrimers.findalltrimers as trimers
from findpointmutationsintrimers.findpointmutationsintrimers import findpointmutationsinupstreamtrimers, \
    findpointmutationsindownstreamtrimers
import findpointmutationsindimers.findpointmutationsindimers as fpmtindmrs


def main():
    if len(sys.argv) != 4:
        sys.exit(
            "The program takes 2 arguments: name of the reference file and name of the reads file name of the outputfile")
    reference: str
    reads: pandas.DataFrame
    length_of_reference: int
    reference, reads, length_of_reference = refutilites.readreference(sys.argv[1], sys.argv[2])
    if length_of_reference == 0:
        sys.exit("The reference file is empty")
    if len(reads) == 0:
        sys.exit("There are no supporting reads in the file")
    reference_dict: dict = dict()

    # More error checks can be added for read quality
    # interface to fastqc etc.
    # interface to GATK / BioPython
    # Index the reference genome to find locations of monomers, dimers and trimers
    reference_dict = refindex.indexreferencetofinddimers(reference, length_of_reference, reference_dict)
    reference_dict = refindex.indexreferencetofindtrimers(reference, length_of_reference, reference_dict)
    reference_dict = refindex.indexreferencetofindmonomers(reference, length_of_reference, reference_dict)

    # After indexing reference genome
    index_of_read_most_upstream: int = int(reads["position"].min())

    # No reads before this location along the reference genome; Hence cannot find point mutations
    assert isinstance(index_of_read_most_upstream, int)
    readpos: int = int(index_of_read_most_upstream)

    variantcalls: TextIO = open(sys.argv[3], "w")
    writer = csv.writer(variantcalls)
    writer.writerow("Position,Ref,Alt,Count")
    # start with the position of the read most upstream of the reference
    # aggregate all supporting reads for that position
    while readpos < length_of_reference:

        print("Checking for variants in position ", readpos)
        length_of_longest_read = 0
        supporting_reads_pointmutations: pandas.DataFrame.astype = fsr.findsupportingreads(reads, readpos)

        # if there are no reads at that position, advance the counter
        if supporting_reads_pointmutations.empty:
            # advance by one
            length_of_longest_read = 1
        else:
            # find the length of the longest and shortest reads
            length_of_longest_read: int = supporting_reads_pointmutations["sequence"].str.len().max()
            lensmallestread: int = int(supporting_reads_pointmutations["sequence"].str.len().min())
            supporting_reads_pointmutations["variant_index"] = ""
            supporting_reads_pointmutations["variant"] = ""
            supporting_reads_pointmutations["3bpdownstream"] = ""
            supporting_reads_pointmutations["2bpdownstream"] = ""
            supporting_reads_pointmutations["3bpupstream"] = ""
            supporting_reads_pointmutations["2bpupstream"] = ""
            supporting_reads_pointmutations["upstreamtrimerindex"] = ""
            supporting_reads_pointmutations["upstreamdimerindex"] = ""
            supporting_reads_pointmutations["downstreamtrimerindex"] = ""
            supporting_reads_pointmutations["downstreamdimerindex"] = ""

            # for each position in the set of supporting reads, calculate variant index, variant,
            # 3 bp upstream, 3 bp downstream, 2 bp upstream, 2 bp downstream
            for counter in range(readpos, readpos + lensmallestread):
                supporting_reads_pointmutations["variant_index"] = counter - supporting_reads_pointmutations["position"]
                # for each column in the set of reads,
                # iterate through all reads and
                # get "Variant" as the base at that position in the read
                supporting_reads_pointmutations["variant"] = [
                    supportingread["sequence"][supportingread["variant_index"]] for rownum, supportingread in
                    supporting_reads_pointmutations.iterrows()]

                #######################################################################################################
                # Looking for point mutations 1, 2, 3 bases long along the read length
                if counter >= readpos + 2:
                    supporting_reads_pointmutations["2bpupstream"], supporting_reads_pointmutations[
                        "upstreamdimerindex"] = dimers.findupstreamdimers(supporting_reads_pointmutations, counter)
                if counter <= ((readpos + lensmallestread) - 2):
                    supporting_reads_pointmutations["2bpdownstream"], supporting_reads_pointmutations[
                        "downstreamdimerindex"] = dimers.finddownstreamdimers(supporting_reads_pointmutations, counter)

                #######################################################################################################
                if counter >= readpos + 3:
                    supporting_reads_pointmutations["3bpupstream"], supporting_reads_pointmutations[
                        "upstreamtrimerindex"] = trimers.finddownstreamtrimers(supporting_reads_pointmutations, counter)

                if counter <= ((readpos + lensmallestread) - 3):
                    supporting_reads_pointmutations["3bpdownstream"], supporting_reads_pointmutations[
                        "downstreamtrimerindex"] = trimers.findupstreamtrimers(supporting_reads_pointmutations, counter)

                #######################################################################################################
                variantpositions = dict()
                # check for trimers as below
                # loop through upstream and downstream trimers
                # and look for positions in reference genome"s dict

                assert isinstance(supporting_reads_pointmutations, pandas.DataFrame)
                assert isinstance(counter, int)
                assert isinstance(int(readpos), int)
                assert isinstance(reference, str)
                assert isinstance(lensmallestread, int)
                (variantpositions, writer) = findpointmutationsinupstreamtrimers(
                    supporting_reads_pointmutations,
                    counter,
                    readpos, reference,
                    reference_dict, variantpositions,
                    writer)
                (variantpositions, writer) = findpointmutationsindownstreamtrimers(
                    supporting_reads_pointmutations,
                    counter,
                    readpos, reference,
                    lensmallestread, reference_dict,
                    variantpositions, writer)
                (variantpositions, writer) = fpmtindmrs.findpointmutationsinupstreamdimers(
                    supporting_reads_pointmutations,
                    counter, readpos,
                    reference,
                    reference_dict, variantpositions,
                    writer)
                (variantpositions, writer) = fpmtindmrs.findpointmutationsindownstreamdimers(
                    supporting_reads_pointmutations,
                    counter,
                    readpos, lensmallestread, reference,
                    reference_dict, variantpositions,
                    writer)

                # At this stage point mutations,  and
                # mutations 3bp long, both upstream and downstream
                # count reference and alt alleles for point mutations
                variantcounts = supporting_reads_pointmutations["variant"].value_counts()
                uniquevariants = supporting_reads_pointmutations["variant"].unique()
            
                for uniq in uniquevariants:
                    if (reference[counter].strip()) != uniq.strip():
                        if variantcounts[uniq] >= 5:
                            if counter not in variantpositions.keys():
                                variantpositions[counter] = variantcounts[uniq]
                                writer.writerow([counter, reference[counter], uniq, variantcounts[uniq]])
                                print([counter, reference[counter], uniq, variantcounts[uniq]])

        readpos = readpos + length_of_longest_read
    variantcalls.close()


#######################################################################################################


# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    t = time.process_time()
    main()
    elapsed_time = time.process_time() - t
    print(elapsed_time, " in seconds to process ", sys.argv[0], " with reads in ", sys.argv[1], " with variants written to ", sys.argv[2])

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
