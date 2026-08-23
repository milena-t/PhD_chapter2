def get_contig_names():
    out_dict = {
        "A_obtectus" : {
            "X" : ["CAVLJG010000002.1","CAVLJG010003236.1","CAVLJG010003544.1","CAVLJG010000099.1","CAVLJG010000155.1","CAVLJG010000244.1","CAVLJG010000377.1","CAVLJG010000488.1",],
            "Y" : ["CAVLJG010000343.1","CAVLJG010002896.1","CAVLJG010000233.1","CAVLJG010000566.1","CAVLJG010000588.1",]
        },
        "B_siliquastri" : {
            "X": ["X"],
            "Y": ["Y"],
        },
        "B_varius" : {
            "X": ["OZ123451.1"],
            "Y": ["OZ123452.1"],
        },
        "C_chinensis" : {
            "X": ["125_quiver","151_quiver","161_quiver","182_quiver","252_quiver","274_quiver","277_quiver","310_quiver","313_quiver","325_quiver","326_quiver","342_quiver","347_quiver","353_quiver","360_quiver","358_quiver","370_quiver","376_quiver","411_quiver","413_quiver","414_quiver","419_quiver","462_quiver","476_quiver","474_quiver","500_quiver","505_quiver","509_quiver","518_quiver","525_quiver","537_quiver","593_quiver","613_quiver","615_quiver","619_quiver","643_quiver","658_quiver","693_quiver","682_quiver","703_quiver","700_quiver","718_quiver","739_quiver","767_quiver","777_quiver","769_quiver","791_quiver","799_quiver","805_quiver","824_quiver","839_quiver","849_quiver","854_quiver","868_quiver","882_quiver","910_quiver","919_quiver","941_quiver","955_quiver","851_quiver","959_quiver","968_quiver","971_quiver","977_quiver","988_quiver","993_quiver","1012_quiver","1024_quiver","1025_quiver","1042_quiver","1000_quiver","1054_quiver"],
            "Y": ["850_quiver","895_quiver","949_quiver","1088_quiver"],
        },
        "C_maculatus" : { 
            "X" : ['scaffold_10','scaffold_14','scaffold_23','scaffold_31','scaffold_34','scaffold_83'],
            "Y" : ['scaffold_26','scaffold_48','scaffold_103','scaffold_112','scaffold_164']
        },
        "D_carinulata" : {
            "X": ["NC_079460.1"], # this was originally falsely annotated as chromosome 1, but it is the X. The annotation on the NCBI has since received a revision
            "Y": ["NC_079473.1"],
        },
        "D_sublineata" : {
            "X": ["NC_079485.1"],
            "Y": ["NC_079486.1"],
        },
    }
    return out_dict


def Cmac_S_L_nonscaffolded_contig_names():
    sex_chr_contigs = { 
        "Y_s" : {
            "X" : ['utg000057l_1','utg000114l_1','utg000139l_1','utg000191l_1','utg000326l_1','utg000359l_1','utg000532l_1','utg000602l_1'],
            "Y" : ['utg000322l_1','utg000312c_1','utg000610l_1','utg001235l_1']
        },
        "Y_l" : {
            "X" : ['utg000006l_1','utg000025l_1','utg000027l_1','utg000128l_1','utg000151l_1','utg000238l_1','utg000327l_1','utg000342l_1','utg000486l_1','utg001894l_1'],
            "Y" : ['utg000049l_1','utg000385c_1','utg001455l_1','utg001921l_1','utg000152l_1']
        }
    }
    return sex_chr_contigs