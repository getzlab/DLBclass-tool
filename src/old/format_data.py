import pandas as pd


def gsm_checker(data):
    # to-do. will check GSM for various things like correct features, formatting, encoding, etc.
    pass


def construct_reduced_winning_version(data):
    gsm_checker(data)

    # This function expects samples in rows, features in columns
    if 'MYD88' in data.index:
        data = data.T

    data = data.astype(float).astype(int)

    # These 5 features removed in winning model.
    # X1Q.AMP
    # X5Q.AMP
    # X4Q35.1.DEL
    # X1Q23.3.AMP
    # X9Q21.13.DEL

    # Don't change any of the recipe code, just add these regions with full zeros.

    data['X1Q.AMP'] = 0
    data['X5Q.AMP'] = 0
    data['X4Q35.1.DEL'] = 0
    data['X1Q23.3.AMP'] = 0
    data['X9Q21.13.DEL'] = 0

    # C1 ##############

    list_bcl6 = ["SV.BCL6", "BCL6"]
    BCL6_ALT = data[list_bcl6].sum(axis=1)

    list_notch2_vec = ["NOTCH2", "SPEN", "DTX1"]
    NOTCH2_vec = data[list_notch2_vec].sum(axis=1)

    list_M88O_vec = ["MYD88.OTHER", "TNFAIP3", "TNIP1", "BCL10", "NFKBIE"]
    M88O_vec = data[list_M88O_vec].sum(axis=1)

    list_C1_vec4 = ["UBE2A", "TMEM30A", "ZEB2", "GNAI2", "X5P.AMP", "X5Q.AMP",
                    "POU2F2", "IKZF3", "X3Q28.DEL", "EBF1", "LYN", "HIST1H2BC",
                    "BCL7A", "CXCR4", "CCDC27", "TUBGCP5", "SMG7", "RHOA",
                    "BTG2", "HLA.B", "ETS1"]
    C1_vec4 = data[list_C1_vec4].sum(axis=1)

    list_CD70_vec = ["CD70", "FAS", "CD58", "B2M", "FADD"]
    CD70_vec = data[list_CD70_vec].sum(axis=1)

    # C2 ##############

    list_tp53 = ["TP53", "X17P.DEL"]
    TP53_biallelic = data[list_tp53].sum(axis=1)

    X21Q_AMP = data["X21Q.AMP"]

    list_C2_sum_arm = ["X17P.DEL", "X21Q.AMP", "X11Q.AMP", "X6P.AMP", "X11P.AMP", "X6Q.DEL", "X7P.AMP", "X13Q.AMP", "X7Q.AMP",
                       "X3Q.AMP", "X5P.AMP", "X5Q.AMP", "X18P.AMP", "X3P.AMP", "X19Q.AMP", "X9Q.AMP", "X1Q.AMP", "X12P.AMP"]
    Sum_C2_ARM = data[list_C2_sum_arm].sum(axis=1)

    list_C2_sum_focal = ["X1P36.11.DEL", "X1P31.1.DEL", "X1P13.1.DEL", "X2Q22.2.DEL",
                         "X16Q12.1.DEL", "X14Q32.31.DEL", "X1P36.32.DEL", "X4Q35.1.DEL",
                         "X9Q21.13.DEL", "X15Q15.3.DEL", "X4Q21.22.DEL", "X9P21.3.DEL",
                         "X8Q24.22.AMP", "X12P13.2.DEL", "X2P16.1.AMP", "X8Q12.1.DEL",
                         "X19P13.2.DEL", "X17Q25.1.DEL", "X1Q42.12.DEL", "X3P21.31.DEL",
                         "X18Q23.DEL", "X19P13.3.DEL", "X13Q34.DEL", "X7Q22.1.AMP",
                         "X10Q23.31.DEL", "X9P24.1.AMP", "X1Q23.3.AMP", "X3Q28.AMP",
                         "X11Q23.3.AMP", "X17Q24.3.AMP", "X3Q28.DEL", "X13Q14.2.DEL",
                         "X18Q21.32.AMP", "X19Q13.32.1.DEL", "X6P21.1.AMP", "X18Q22.2.AMP",
                         "EP300", "CD274", "ZNF423"]
    Sum_C2_FOCAL = data[list_C2_sum_focal].sum(axis=1)

    # C3 ##############

    list_bcl2 = ["BCL2", "SV.BCL2"]
    BCL2_combined = data[list_bcl2].sum(axis=1)

    list_CREBBP = ["CREBBP", "EZH2", "KMT2D", "EP300"]
    CREBBP_vec = data[list_CREBBP].sum(axis=1)

    list_GNA13 = ["GNA13", "TNFRSF14", "MAP2K1", "MEF2B", "IRF8", "HVCN1", "GNAI2",
                  "MEF2C", "POU2AF1", "RAC2", "X12P.AMP", "X12Q.AMP", "X6Q14.1.DEL"]
    GNA13_vec = data[list_GNA13].sum(axis=1)

    PTEN = data["PTEN"] + data["X10Q23.31.DEL"] + data["X13Q14.2.DEL"]

    SV_MYC = data["SV.MYC"]

    # C4 ##############

    list_Hist_comp = ["HIST1H2AC", "HIST1H1E", "HIST1H1B", "HIST1H2AM",
                      "HIST1H1C", "HIST1H1D", "HIST1H2BC", "HIST1H2BD"]
    Hist_comp = data[list_Hist_comp].sum(axis=1)

    list_SGK1_vec = ["SGK1", "TET2", "NFKBIA", "STAT3", "PTPN6", "BRAF", "KRAS", "CD83",
                     "SF3B1", "CD274", "MEF2C", "KLHL6", "CXCR4", "PTEN", "RAC2", "SESN3"]
    SGK1_vec = data[list_SGK1_vec].sum(axis=1)

    list_DUSP2_vec = ["DUSP2", "SOCS1", "ZFP36L1", "CRIP1", "ACTB", "LTB",
                      "YY1", "ZNF608", "PABPC1", "EEF1A1"]
    DUSP2_vec = data[list_DUSP2_vec].sum(axis=1)

    # C5 ##############

    list_TBL1XR1_vec = ["TBL1XR1", "PIM1", "PRDM1", "ETV6", "ZC3H12A", "BTG1", "BTG2",
                        "IGLL5", "TMSB4X", "GRHPR", "HLA.C", "MYD88", "TOX", "LYN", "POU2F2",
                        "IKZF3", "HLA.A", "ZFP36L1", "CARD11", "SF3B1", "HLA.B",
                        "IRF2BP2", "OSBPL10", "VMP1", "ATP2A2", "CCDC27", "ETS1"]
    TBL1XR1_vec = data[list_TBL1XR1_vec].sum(axis=1)

    MYD88_L265P_CD79B = data["MYD88.L265P"] + data["CD79B"]

    list_Sum_C5_sig = ["X18Q.AMP", "X3Q.AMP", "X3P.AMP", "X19Q13.42.AMP", "X6Q21.DEL",
                       "X18P.AMP", "X19Q.AMP", "X8Q12.1.DEL", "X6Q14.1.DEL", "X19P13.2.DEL",
                       "X9P21.3.DEL", "X18Q21.32.AMP", "X18Q22.2.AMP", "X13Q.AMP",
                       "X1Q42.12.DEL", "X9Q.AMP", "X1Q32.1.AMP", "X6P21.33.DEL"]
    Sum_C5_CNA = data[list_Sum_C5_sig].sum(axis=1)

    # MISC ############

    CN_2P16_1_AMP = data["X2P16.1.AMP"]

    reduced_feature_dict = \
        {'BCL6_ALT': BCL6_ALT,
         'NOTCH2_vec': NOTCH2_vec,
         'M88O_vec': M88O_vec,
         'C1_vec4': C1_vec4,
         'CD70_vec': CD70_vec,
         'TP53_biallelic': TP53_biallelic,
         'X21Q_AMP': X21Q_AMP,
         'Sum_C2_ARM': Sum_C2_ARM,
         'Sum_C2_FOCAL': Sum_C2_FOCAL,
         'BCL2_combined': BCL2_combined,
         'CREBBP_vec': CREBBP_vec,
         'GNA13_vec': GNA13_vec,
         'PTEN': PTEN,
         'SV_MYC': SV_MYC,
         'Hist_comp': Hist_comp,
         'SGK1_vec': SGK1_vec,
         'DUSP2_vec': DUSP2_vec,
         'CN_2P16_1_AMP': CN_2P16_1_AMP,
         'TBL1XR1_vec': TBL1XR1_vec,
         'MYD88_L265P_CD79B': MYD88_L265P_CD79B,
         'Sum_C5_CNA': Sum_C5_CNA}

    data = pd.DataFrame.from_dict(reduced_feature_dict)

    return data
