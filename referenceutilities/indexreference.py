

################ Index the reference assembly (in this case ECGR) #################################
# function that indexes the reference reads as 1-3 bp
def indexreferencetofindtrimers(reference: str, length_of_reference: int, reference_dict: dict) -> dict:
    """

    :type length_of_reference: int
    :type reference_dict: dict
    :rtype: bool
    :param reference:
    :param length_of_reference:
    """
    # break up reference into trimers
    counter: int
    for counter in range(1, length_of_reference - 3):
        # move along the length of the genome, and find all the trimers, and add to hash
        # hash creates a key of trimer, and the values are the positions where the trimer occurs
        if reference[counter:counter + 3] in reference_dict.keys():
            reference_dict[str(reference[counter:counter + 3])].append(counter)
        else:
            reference_dict[reference[counter:counter + 3]] = [1]
    return reference_dict


#######################################################################################################
def indexreferencetofindmonomers(reference: str, length_of_reference: int, reference_dict: dict) -> dict:
    """

    :type length_of_reference: int
    :param length_of_reference:
    :type reference: str
    :param reference:
    :rtype: object
    :type reference_dict: dict
    """
    # break up the reference into single bases and
    # build a hash of the single base, and all the locations where it occurs
    counter: int
    for counter in range(1, length_of_reference):
        if reference[counter] in reference_dict.keys():
            reference_dict[str(reference[counter])].append(counter)
        else:
            reference_dict[reference[counter]] = [1]
    return reference_dict


#######################################################################################################
# function that indexes the reference reads as 1-3 bp
def indexreferencetofinddimers(reference: str, length_of_reference: int, reference_dict: dict) -> dict:
    """

    :rtype: dict
    :param reference_dict:
    :type reference: str
    :param reference:
    :param length_of_reference: int
    """
    # break up reference into dimers
    # build a hash of dimer and its location
    counter: int
    for counter in range(1, length_of_reference - 2):
        if reference[counter:counter + 2] in reference_dict.keys():
            reference_dict[str(reference[counter:counter + 2])].append(counter)
        else:
            reference_dict[reference[counter:counter + 2]] = [1]
    return reference_dict
