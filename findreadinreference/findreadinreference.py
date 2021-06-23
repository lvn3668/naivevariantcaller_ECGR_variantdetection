
#######################################################################################################
# function that searches for trimer in reference genome that was indexed


def searchtrimerinreference(trimer: str, variant_index: int, reference_dict: dict()) -> bool:
    """

    :rtype: bool
    :param trimer:
    :param variant_index:
    :return:
    """
    if trimer in reference_dict.keys():
        trimerlocs = reference_dict.get(trimer)
        if variant_index in trimerlocs:
            return True
        else:
            return False

    return False
