import math


def model_dominant(row_list, effect_size_type):
    effect_size = 0
    var = 0
    aa1 = float(row_list[0])
    ab1 = float(row_list[1])
    bb1 = float(row_list[2])
    aa0 = float(row_list[3])
    ab0 = float(row_list[4])
    bb0 = float(row_list[5])
    R = aa1 + ab1 + bb1
    S = aa0 + ab0 + bb0
    N = R + S
    n2 = bb0 + bb1
    n1 = ab1 + ab0
    if (effect_size_type == 'OR'):
        effect_size = math.log(((aa1 + ab1) / (aa0 + ab0)) / (bb1 / bb0))
        var = 1 / (aa1 + ab1) + 1 / (aa0 + ab0) + 1 / bb1 + 1 / bb0

    # if (effect_size_type == 'RR'):
    #     effect_size = ((aa1 + ab1) * (bb1 + bb0)) / (bb1 * (aa1 + ab1 + aa0 + ab0))
    #
    #     var = 1 / (aa1 + ab1) - 1 / (aa1 + ab1 + aa0 + ab0) + 1 / bb1 - 1 / bb1 / bb0
    #
    # if (effect_size_type == 'RISK_DIFF'):
    #     effect_size = (aa1 + ab1) / (aa1 + ab1 + aa0 + ab0) - bb1 / (bb1 + bb0)
    #
    #     var = ((aa1 + ab1) * (aa0 + ab0)) / (aa1 + ab1 + aa0 + ab0) ** 3 - bb1 / (bb1 + bb0)

    if (effect_size_type == 'CATT'):
        effect_size = 1 / N * ((S * (ab1 + bb1)) - (R * (ab0 + bb0)))

        var = (R * S / N) * (((n1) + n2) / N - (((n1 + n2) / N) * ((n1 + n2) / N)))

    return effect_size, var


def model_additive(row_list, effect_size_type):  # allelic model
    effect_size = 0
    var = 0
    aa1 = float(row_list[0])
    ab1 = float(row_list[1])
    bb1 = float(row_list[2])
    aa0 = float(row_list[3])
    ab0 = float(row_list[4])
    bb0 = float(row_list[5])
    R = aa1 + ab1 + bb1
    S = aa0 + ab0 + bb0
    N = R + S
    n2 = bb0 + bb1
    n1 = ab1 + ab0

    if (effect_size_type == 'OR'):
        effect_size = math.log(((2 * bb1 + ab1) / (2 * bb0 + ab0)) / ((2 * aa1 + ab1) / (2 * aa0 + ab0)))
        var = 1 / (2 * bb1 + ab1) + 1 / (2 * bb0 + ab0) + 1 / (2 * aa1 + ab1) + 1 / (2 * aa0 + ab0)

    # if (effect_size_type == 'RR'):
    #     effect_size = math.log(((2 * bb1 + ab1) * ((2 * aa1 + ab1) + (2 * aa0 + ab0))) / (
    #             (2 * aa1 + ab1) * (2 * bb1 + ab1 + 2 * bb0 + ab0)))
    #
    #     var = 1 / (2 * bb1 + ab1) - 1 / ((2 * bb1 + ab1) + (2 * bb0 + ab0)) + 1 / (2 * aa1 + ab1) - 1 / (
    #             2 * aa1 + ab1 + 2 * aa0 + ab0)

    # if (effect_size_type == 'RISK_DIFF'):
    #     effect_size = (2 * bb1 + ab1) / (2 * bb1 + ab1 + 2 * bb0 + ab0) - (2 * aa1 + ab1) / (
    #             2 * aa1 + ab1 + 2 * aa0 + ab0)
    #
    #     var = ((2 * bb1 + ab1) * (2 * bb0 + ab0)) / (2 * bb1 + ab1 + 2 * bb0 + ab0) ** 3 + (
    #             (2 * aa1 + ab1) * (2 * aa0 + ab0)) / (2 * aa1 + ab1 + 2 * aa0 + ab0) ** 3

    if (effect_size_type == 'CATT'):
        effect_size = 1 / N * ((S * (0.5 * ab1 + bb1)) - (R * (0.5 * ab0 + bb0)))

        var = (R * S / N) * (((0.5 * 0.5 * n1) + n2) / N - (((0.5 * n1 + n2) / N) * ((0.5 * n1 + n2) / N)))

    return effect_size, var


def model_recessive(row_list, effect_size_type):
    effect_size = 0
    var = 0
    aa1 = float(row_list[0])
    ab1 = float(row_list[1])
    bb1 = float(row_list[2])
    aa0 = float(row_list[3])
    ab0 = float(row_list[4])
    bb0 = float(row_list[5])
    R = aa1 + ab1 + bb1
    S = aa0 + ab0 + bb0
    N = R + S
    n2 = bb0 + bb1
    n1 = ab1 + ab0

    # print(aa1,ab1,bb1,aa0,ab0,bb0)
    if (effect_size_type == 'OR'):
        effect_size = math.log(((bb1 + ab1) /
                                (bb0 + ab0)) /
                               (aa1 / aa0))

        var = 1 / (bb1 + ab1) + 1 / (bb0 + ab0) + 1 / aa1 + 1 / aa0
    #
    # if (effect_size_type == 'RR'):
    #     effect_size = ((aa1 + ab1) * (bb1 + bb0)) / (bb1 * (aa1 + ab1 + ab0 + ab0))
    #     var = 1 / (aa1 + ab1) - 1 / (aa1 + ab1 + aa0 + ab0) + 1 / bb1 - 1 / (bb1 + bb0)
    #
    # if (effect_size_type == 'RISK_DIFF'):
    #     effect_size = (aa1 + ab1) / (aa1 + ab1 + aa0 + ab0) - bb1 / (bb1 + bb0)
    #
    #     var = ((aa1 + ab1) * (aa0 + ab0)) / (aa1 + ab1 + aa0 + ab0) ** 3 + (bb1 * bb0) / (bb1 + bb0) ** 3

    if (effect_size_type == 'CATT'):
        effect_size = 1 / N * (S * bb1 - R * bb0)

        var = (R * S / N) * ((n2 / N) - ((n2 / N) * (n2 / N)))
    return effect_size, var
