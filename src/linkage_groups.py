from enum import Enum


class LinkageGroup(str, Enum):
    """
    Assign all chromosome names to their linkage groups in bruchini
    Here listed are the Bsil and Dcar chromosome names respectively
    """
    BLG1 = "1",
    BLG2 = "2",
    BLG3 = "3",
    BLG4 = "4",
    BLG5 = "5",
    BLG6 = "6",
    BLG7 = "7",
    BLG8 = "8",
    BLG9 = "9",
    BLGX = "X",
    BLGY = "Y",
    DLG1 = "NC_079476.1",
    DLG2 = "NC_079478.1",
    DLG3 = "NC_079479.1",
    DLG4 = "NC_079480.1",
    DLG5 = "NC_079481.1",
    DLG6 = "NC_079477.1",
    DLG7 = "NC_079484.1",
    DLGX = "NC_079460.1",
    DLGY = "NC_079473.1",


contigname_to_linkagegroup = {
    # B. siliquastri
    "1" : LinkageGroup.BLG1,
    "2" : LinkageGroup.BLG2,
    "3" : LinkageGroup.BLG3,
    "4" : LinkageGroup.BLG4,
    "5" : LinkageGroup.BLG5,
    "6" : LinkageGroup.BLG6,
    "7" : LinkageGroup.BLG7,
    "8" : LinkageGroup.BLG8,
    "9" : LinkageGroup.BLG9,
    "X" : LinkageGroup.BLGX,
    "Y" : LinkageGroup.BLGY,
    # B. varius
    "OZ123444.1" : LinkageGroup.BLG1,
    "OZ123447.1" : LinkageGroup.BLG2,
    "OZ123445.1" : LinkageGroup.BLG3,
    "OZ123446.1" : LinkageGroup.BLG4,
    "OZ123448.1" : LinkageGroup.BLG5,
    "OZ123450.1" : LinkageGroup.BLG6,
    "OZ123449.1" : LinkageGroup.BLG7,
    "OZ123443.1" : LinkageGroup.BLG8,
    "OZ123451.1" : LinkageGroup.BLGX,
    "OZ123452.1" : LinkageGroup.BLGY,
    # "OZ123443.1" : LinkageGroup.BLG9, # I will just pick the fusion for LG8
    # A. obtectus 
    "CAVLJG010000005.1" : LinkageGroup.BLG1,
    "CAVLJG010000003.1" : LinkageGroup.BLG2,
    "CAVLJG010000008.1" : LinkageGroup.BLG3,
    "CAVLJG010000004.1" : LinkageGroup.BLG4,
    "CAVLJG010000007.1" : LinkageGroup.BLG5,
    "CAVLJG010000001.1" : LinkageGroup.BLG6,
    "CAVLJG010000006.1" : LinkageGroup.BLG7,
    "CAVLJG010000009.1" : LinkageGroup.BLG8,
    "CAVLJG010000010.1" : LinkageGroup.BLG9,
    "CAVLJG010000002.1" : LinkageGroup.BLGX,
    "CAVLJG010003236.1" : LinkageGroup.BLGX,
    "CAVLJG010003544.1" : LinkageGroup.BLGX,
    "CAVLJG010000099.1" : LinkageGroup.BLGX,
    "CAVLJG010000155.1" : LinkageGroup.BLGX,
    "CAVLJG010000244.1" : LinkageGroup.BLGX,
    "CAVLJG010000377.1" : LinkageGroup.BLGX,
    "CAVLJG010000488.1" : LinkageGroup.BLGX,
    "CAVLJG010000343.1" : LinkageGroup.BLGY,
    "CAVLJG010002896.1" : LinkageGroup.BLGY,
    "CAVLJG010000233.1" : LinkageGroup.BLGY,
    "CAVLJG010000566.1" : LinkageGroup.BLGY,
    "CAVLJG010000588.1" : LinkageGroup.BLGY,
    # C. maculatus
    "scaffold_1" : LinkageGroup.BLG1,
    "scaffold_6" : LinkageGroup.BLG2,
    "scaffold_4" : LinkageGroup.BLG3,
    "scaffold_13" : LinkageGroup.BLG4,
    "scaffold_7" : LinkageGroup.BLG4,
    "scaffold_15" : LinkageGroup.BLG4,
    "scaffold_3" : LinkageGroup.BLG5,
    "scaffold_5" : LinkageGroup.BLG6,
    "scaffold_2" : LinkageGroup.BLG7,
    "scaffold_8" : LinkageGroup.BLG8,
    "scaffold_12" : LinkageGroup.BLG8,
    "scaffold_11" : LinkageGroup.BLG9,
    "scaffold_9" : LinkageGroup.BLG9,
    # D. carinulata
    "NC_079472.1" :	LinkageGroup.DLG1,
    "NC_079461.1" :	LinkageGroup.DLG1,
    "NC_079467.1" :	LinkageGroup.DLG1,
    "NC_079462.1" :	LinkageGroup.DLG2,
    "NC_079463.1" :	LinkageGroup.DLG3,
    "NC_079464.1" :	LinkageGroup.DLG4,
    "NC_079465.1" :	LinkageGroup.DLG5,
    "NC_079466.1" :	LinkageGroup.DLG6,
    "NC_079469.1" :	LinkageGroup.DLG6,
    "NC_079470.1" :	LinkageGroup.DLG6,
    "NC_079471.1" :	LinkageGroup.DLG7,
    "NC_079460.1" : LinkageGroup.DLGX,
    "NC_079473.1" : LinkageGroup.DLGY,
    # D. sublineata
    "NC_079474.1" : LinkageGroup.DLG1,
    "NC_079475.1" : LinkageGroup.DLG1,
    "NC_079476.1" : LinkageGroup.DLG1,
    "NC_079478.1" : LinkageGroup.DLG2,
    "NC_079479.1" : LinkageGroup.DLG3,
    "NC_079480.1" : LinkageGroup.DLG4,
    "NC_079481.1" : LinkageGroup.DLG5,
    "NC_079477.1" : LinkageGroup.DLG6,
    "NC_079483.1" : LinkageGroup.DLG6,
    "NC_079484.1" : LinkageGroup.DLG7,
    "NC_079485.1" : LinkageGroup.DLGX,
    "NC_079486.1" : LinkageGroup.DLGY,
}


def assign_linkagegoup(s):
    """
    sort the contig/scaffold name into one of the existing linkage groups
    """
    return contigname_to_linkagegroup.get(s)
