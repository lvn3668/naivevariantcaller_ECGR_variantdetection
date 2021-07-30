# Author: Lalitha Viswanathan
# Python Package to index the reference genome
# Builds index of all nmers (n=1-3) as point mutations is defined as
# 1-3 bp mutations
# 3 consecutive point mutations to be called as 1 variant
# This code reads in reference genome into dataframe
# Reads all read fragments
# Returns reference, reads and length of reference
from typing import TextIO

import pandas
import os


def readreference(referencefilename: str, readsfilename: str) -> (str, pandas.DataFrame, int):
    """

    :rtype: object
    :type readsfilename: str
    :type referencefilename: str
    """
    file: TextIO
    print("Reference file name is ", referencefilename)
    with open(os.path.join(os.environ['DATA'], referencefilename), "r") as file:
        reference: str = file.read()  # reference sequence as a string
    file.close()

    print("Reads file name is ", readsfilename)
    reads = pandas.read_csv(os.path.join(os.environ['DATA'], readsfilename),
                            index_col="read_id")
    length_of_reference: int = len(reference)
    return reference, reads, length_of_reference
