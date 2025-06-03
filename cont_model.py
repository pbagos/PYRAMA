import math



def calcVariance(sd1,n1,sd2,n2):
    return(sd1/n1 + sd2/n2)

def modelWeightedPar(n1,x1,n2,x2):
    return(((n1 * x1) + (n2 * x2)) / (n1 + n2))


def calcCombinedSD(n1,n2,sd1,sd2):
    return(math.sqrt (((n1 - 1) * sd1 * sd1 + (n2 - 1) * sd2 * sd2) / (n1 + n2)))



def cont_model_dominant(row_list):
    xaa = row_list[0]
    sdaa = row_list[1]
    naa = row_list[2]
    xab = row_list[3]
    sdab = row_list[4]
    nab = row_list[5]
    xbb = row_list[6]
    sdbb = row_list[7]
    nbb = row_list[8]

    rmd_variance = calcVariance(calcCombinedSD(nbb,nab,sdbb,sdab),
                                nbb +nab -2,
                                sdaa*sdaa,
                                naa)
    rmd = modelWeightedPar(nbb,xbb,nab,xab) - xaa
    return rmd, rmd_variance

def cont_model_additive(row_list):
    xaa = row_list[0]
    sdaa = row_list[1]
    naa = row_list[2]
    xab = row_list[3]
    sdab = row_list[4]
    nab = row_list[5]
    xbb = row_list[6]
    sdbb = row_list[7]
    nbb = row_list[8]

    rmd_variance = calcVariance(calcCombinedSD(2*nbb,nab,sdbb,sdab),
                                2*nbb +nab -2,
                                calcCombinedSD(2 * naa,nab,  sdaa, sdab),
                                2 * naa + nab- 2)
    rmd = modelWeightedPar(2*nbb,xbb,nab,xab) - modelWeightedPar(2 * naa, xaa, nab, xab)
    return rmd, rmd_variance

def cont_model_recessive(row_list):
    xaa = row_list[0]
    sdaa = row_list[1]
    naa = row_list[2]
    xab = row_list[3]
    sdab = row_list[4]
    nab = row_list[5]
    xbb = row_list[6]
    sdbb = row_list[7]
    nbb = row_list[8]

    rmd = xbb - modelWeightedPar(naa, xaa,nab, xab)

    rmd_variance = calcVariance((sdbb * sdbb) , nbb,
                calcCombinedSD(naa,nab,  sdaa, sdab) ,
                (naa + nab - 2))


    return rmd, rmd_variance

