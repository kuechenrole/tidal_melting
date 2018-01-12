import numpy as np
from netCDF4 import Dataset

# This code is adapted from the matlab code
# "LP Bathymetry" by Mathieu Dutour Sikiric
# http://drobilica.irb.hr/~mathieu/Bathymetry/index.html
# For a description of the method, see
# M. Dutour Sikiric, I. Janekovic, M. Kuzmic, A new approach to
# bathymetry smoothing in sigma-coordinate ocean models, Ocean
# Modelling 29 (2009) 128--136.
def RoughnessMatrix(DEP, MSK):
    """
    RoughMat=GRID_RoughnessMatrix(DEP, MSK)

    ---DEP is the bathymetry of the gridcd
    ---MSK is the mask of the grid
    """

    eta_rho, xi_rho = DEP.shape

    Umat = np.array([[0, 1],
                                            [1, 0],
                                            [0, -1],
                                            [-1, 0]])

    RoughMat = np.zeros(DEP.shape)

    for iEta in range(1,eta_rho-1):
        for iXi in range(1,xi_rho-1):
            if (MSK[iEta,iXi] == 1):
                rough = 0
                for i in range(4):
                    iEtaB = iEta + Umat[i,0]
                    iXiB = iXi + Umat[i,1]
                    if (MSK[iEtaB,iXiB] == 1):
                        dep1 = DEP[iEta,iXi]
                        dep2 = DEP[iEtaB,iXiB]
                        delta = abs((dep1 - dep2) / (dep1 + dep2))
                        rough = np.maximum(rough, delta)

                RoughMat[iEta,iXi] = rough


    return RoughMat

# %load bathy_smoothing
import numpy as np

# This code is adapted from the matlab code
# "LP Bathymetry" by Mathieu Dutour Sikiric
# http://drobilica.irb.hr/~mathieu/Bathymetry/index.html
# For a description of the method, see
# M. Dutour Sikiric, I. Janekovic, M. Kuzmic, A new approach to
# bathymetry smoothing in sigma-coordinate ocean models, Ocean
# Modelling 29 (2009) 128--136.

def smoothing_Positive_rx0(MSK, Hobs, rx0max):
    """
    This program use the direct iterative method from Martinho and Batteen (2006)
    The bathymetry is optimized for a given rx0 factor by increasing it.

    Usage:
    RetBathy = smoothing_Positive_rx0(MSK, Hobs, rx0max)

    ---MSK(eta_rho,xi_rho) is the mask of the grid
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    """

    eta_rho, xi_rho = Hobs.shape

    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    RetBathy = Hobs.copy()

    nbModif = 0
    tol = 0.000001

    while(True):
        IsFinished = 1
        for iEta in range(eta_rho):
            for iXi in range(xi_rho):
                if (MSK[iEta,iXi] == 1):
                    for ineigh in range(4):
                        iEtaN = iEta + ListNeigh[ineigh,0]
                        iXiN = iXi + ListNeigh[ineigh,1]
                        if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                                and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                            LowerBound = RetBathy[iEtaN,iXiN] * (1-rx0max)/(1+rx0max)
                            if ((RetBathy[iEta,iXi] - LowerBound) < -tol):
                                IsFinished = 0
                                RetBathy[iEta,iXi] = LowerBound
                                nbModif = nbModif + 1

        if (IsFinished == 1):
            break

    print('     nbModif=', nbModif)

    return RetBathy



def smoothing_Negative_rx0(MSK, Hobs, rx0max):
    """
    This program use an opposite methode to the direct iterative method from 
    Martinho and Batteen (2006). This program optimizes the bathymetry for 
    a given rx0 factor by decreasing it.

    Usage:
    RetBathy = smoothing_Negative_rx0(MSK, Hobs, rx0max)

    ---MSK(eta_rho,xi_rho) is the mask of the grid
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    """

    eta_rho, xi_rho = Hobs.shape

    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    RetBathy = Hobs.copy()

    nbModif = 0
    tol = 0.000001

    while(True):
        IsFinished = 1
        for iEta in range(eta_rho):
            for iXi in range(xi_rho):
                if (MSK[iEta, iXi] == 1):
                    for ineigh in range(4):
                        iEtaN = iEta + ListNeigh[ineigh,0]
                        iXiN = iXi + ListNeigh[ineigh,1]
                        if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                                and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                            UpperBound = RetBathy[iEtaN, iXiN] * (1+rx0max)/(1-rx0max)
                            if (RetBathy[iEta,iXi] > (UpperBound + tol)):
                                IsFinished = 0
                                RetBathy[iEta, iXi] = UpperBound
                                nbModif = nbModif + 1

        if (IsFinished == 1):
            break

    print('     nbModif=', nbModif)
    return RetBathy



def smoothing_PositiveVolume_rx0(MSK, Hobs, rx0max, AreaMatrix):
    """
    This program use the direct iterative method from Martinho and Batteen (2006)
    The bathymetry is optimized for a given rx0 factor by increasing it. All depth 
    are then multiplied by the coeficient K = Vol_init/Vol_final in order to 
    insure volume conservation.

    Usage:
    RetBathy = smoothing_Positive_rx0(MSK, Hobs, rx0max, AreaMatrix)

    ---MSK(eta_rho,xi_rho) is the mask of the grid
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    ---AreaMatrix(eta_rho,xi_rho) is the matrix of areas at
       rho point
    """

    eta_rho, xi_rho = Hobs.shape

    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    WorkBathy = Hobs.copy()

    nbModif = 0
    tol = 0.000001

    while(True):
        IsFinished = 1
        for iEta in range(eta_rho):
            for iXi in range(xi_rho):
                if (MSK[iEta,iXi] == 1):
                    for ineigh in range(4):
                        iEtaN = iEta + ListNeigh[ineigh,0]
                        iXiN = iXi + ListNeigh[ineigh,1]
                        if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                                and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                            LowerBound = WorkBathy[iEtaN,iXiN] * (1-rx0max)/(1+rx0max)
                            if ((WorkBathy[iEta,iXi] - LowerBound) < -tol):
                                IsFinished = 0
                                WorkBathy[iEta,iXi] = LowerBound
                                nbModif = nbModif + 1

        if (IsFinished == 1):
            break

    print('     nbModif=', nbModif)

    VolOrig=0
    VolWork=0
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            if (MSK[iEta, iXi] == 1):
                VolOrig = VolOrig + AreaMatrix[iEta,iXi] * Hobs[iEta,iXi]
                VolWork = VolWork + AreaMatrix[iEta,iXi] * WorkBathy[iEta,iXi]

    print('volorig/volwork = ',((VolOrig / VolWork)))

    RetBathy = WorkBathy * (VolOrig / VolWork)
    DeltaBathy = (WorkBathy  - WorkBathy * (VolOrig / VolWork)) * MSK
    RetBathy = WorkBathy - DeltaBathy

    return RetBathy



def smoothing_NegativeVolume_rx0(MSK, Hobs, rx0max, AreaMatrix):
    """
    This program use an opposite methode to the direct iterative method from 
    Martinho and Batteen (2006). This program optimizes the bathymetry for 
    a given rx0 factor by decreasing it. All depth are then multiplied by 
    the coeficient K = Vol_init/Vol_final in order to insure volume conservation.

    Usage:
    RetBathy = smoothing_Negative_rx0(MSK, Hobs, rx0max, AreaMatrix)

    ---MSK(eta_rho,xi_rho) is the mask of the grid
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    ---AreaMatrix(eta_rho,xi_rho) is the matrix of areas at
       rho point
    """

    eta_rho, xi_rho = Hobs.shape

    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    WorkBathy = Hobs.copy()

    nbModif = 0
    tol = 0.000001

    while(True):
        #IsFinished = 1
        for iEta in range(eta_rho):
            #print("processing j: ",iEta)
            for iXi in range(xi_rho):
                if (MSK[iEta, iXi] == 1):
                    for ineigh in range(4):
                        iEtaN = iEta + ListNeigh[ineigh,0]
                        iXiN = iXi + ListNeigh[ineigh,1]
                        if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                                and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                            UpperBound = WorkBathy[iEtaN, iXiN] * (1+rx0max)/(1-rx0max)
                            if (WorkBathy[iEta,iXi] > (UpperBound + tol)):
                                IsFinished = 0
                                WorkBathy[iEta, iXi] = UpperBound
                                nbModif = nbModif + 1

        #if (IsFinished == 1):
        break

    print('     nbModif=', nbModif)

    VolOrig=0
    VolWork=0
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            if (MSK[iEta, iXi] == 1):
                VolOrig = VolOrig + AreaMatrix[iEta,iXi] * Hobs[iEta,iXi]
                VolWork = VolWork + AreaMatrix[iEta,iXi] * WorkBathy[iEta,iXi]

    print('     percent lowered to conserve Volume: ',(VolOrig / VolWork))

    DeltaBathy = (WorkBathy  - WorkBathy * (VolOrig / VolWork)) * MSK
    RetBathy = WorkBathy - DeltaBathy

    return RetBathy,nbModif
    #return WorkBathy


def smoothing_PlusMinus_rx0(MSK, Hobs, rx0max, AreaMatrix):
    """
    This program use the Mellor-Ezer-Oey method (Mellor et al., 1994).
    The bathymetry is optimized for a given rx0 factor by doing a sequence
    of increase/decrease at adjacent cells.

    Usage:
    RetBathy, HmodifVal, ValueFct = smoothing_PlusMinus_rx0(MSK, Hobs, rx0max, AreaMatrix)

    ---MSK(eta_rho,xi_rho) is the mask of the grid
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    ---AreaMatrix(eta_rho,xi_rho) is the matrix of areas at
       rho-points.
    """

    eta_rho, xi_rho = Hobs.shape

    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    RetBathy = Hobs.copy()

    HmodifVal = 0
    TheMultiplier = (1 - rx0max) / (1 + rx0max)
    tol = 0.000001
    ValueFct = 0
    nbModif_total=0

    while(True):
        IsFinished = 1
        nbModif=0
        for iEta in range(eta_rho):
            for iXi in range(xi_rho):
                if (MSK[iEta, iXi] == 1):
                    Area = AreaMatrix[iEta, iXi]
                    for ineigh in range(4):
                        iEtaN = iEta + ListNeigh[ineigh,0]
                        iXiN = iXi + ListNeigh[ineigh,1]
                        if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                            and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                            AreaN = AreaMatrix[iEtaN,iXiN]
                            LowerBound = RetBathy[iEtaN,iXiN] * TheMultiplier
                            if ((RetBathy[iEta,iXi] - LowerBound) < -tol):
                                IsFinished = 0
                                h = (TheMultiplier * RetBathy[iEtaN,iXiN] - RetBathy[iEta,iXi]) \
                                         / (AreaN + TheMultiplier * Area)
                                RetBathy[iEta,iXi] = RetBathy[iEta,iXi] + AreaN * h
                                RetBathy[iEtaN,iXiN] = RetBathy[iEtaN,iXiN] - Area * h
                                HmodifVal = HmodifVal + abs(h)
                                ValueFct = ValueFct + abs(h) * (Area + AreaN)
                                nbModif+=1
                                nbModif_total+=1
        print('     nbModif=', nbModif)

        if (IsFinished == 1):
            break

    H = AreaMatrix * Hobs * MSK
    TheBathymetry1 = H.sum()
    H = AreaMatrix * RetBathy * MSK
    TheBathymetry2 = H.sum()
    DeltaBathymetry = TheBathymetry1 - TheBathymetry2
    print('DeltaBathymetry = ', DeltaBathymetry)

    return RetBathy, HmodifVal, ValueFct


def smoothing_Laplacian_rx0(MSK, Hobs, rx0max):
    """
    This program use Laplacian filter.
    The bathymetry is optimized for a given rx0 factor by doing an iterated
    sequence of Laplacian filterings.

    Usage:
    RetBathy = smoothing_Laplacian_rx0(MSK, Hobs, rx0max)

    ---MSK(eta_rho,xi_rho) is the mask of the grid
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    """

    eta_rho, xi_rho = Hobs.shape

    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    RetBathy = Hobs.copy()

    tol = 0.00001
    WeightMatrix = np.zeros((eta_rho, xi_rho))
    for iEta in range(eta_rho):
        for iXi in range(xi_rho):
            WeightSum = 0
            for ineigh in range(4):
                iEtaN = iEta + ListNeigh[ineigh,0]
                iXiN = iXi + ListNeigh[ineigh,1]
                if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                      and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                    WeightSum = WeightSum + 1

            WeightMatrix[iEta,iXi] = WeightSum

    Iter = 1
    NumberDones = np.zeros((eta_rho, xi_rho))
    while(True):
        RoughMat = bathy_tools.RoughnessMatrix(RetBathy, MSK)
        Kbefore = np.where(RoughMat > rx0max)
        nbPtBefore = np.size(Kbefore, 1)
        realR = RoughMat.max()
        TheCorrect = np.zeros((eta_rho,xi_rho))
        IsFinished = 1
        nbPointMod = 0
        AdditionalDone = np.zeros((eta_rho, xi_rho))
        for iEta in range(eta_rho):
            for iXi in range(xi_rho):  
                Weight = 0
                WeightSum = 0
                for ineigh in range(4):
                    iEtaN = iEta + ListNeigh[ineigh,0]
                    iXiN = iXi + ListNeigh[ineigh,1]
                    if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                          and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                        Weight = Weight + RetBathy[iEtaN,iXiN]
                        AdditionalDone[iEtaN,iXiN] = AdditionalDone[iEtaN,iXiN] + NumberDones[iEta,iXi]

                TheWeight = WeightMatrix[iEta,iXi]
                WeDo = 0
                if TheWeight > tol:
                    if RoughMat[iEta,iXi] > rx0max:
                        WeDo = 1
                    if NumberDones[iEta,iXi] > 0:
                        WeDo = 1

                if WeDo == 1:
                    IsFinished = 0
                    TheDelta = (Weight - TheWeight * RetBathy[iEta,iXi]) / (2 * TheWeight)
                    TheCorrect[iEta,iXi] = TheCorrect[iEta,iXi] + TheDelta
                    nbPointMod = nbPointMod + 1
                    NumberDones[iEta,iXi] = 1

        NumberDones = NumberDones + AdditionalDone
        RetBathy = RetBathy + TheCorrect
        NewRoughMat = bathy_tools.RoughnessMatrix(RetBathy, MSK)
        Kafter = np.where(NewRoughMat > rx0max)
        nbPtAfter = np.size(Kafter, 1)
        TheProd = (RoughMat > rx0max) * (NewRoughMat > rx0max)
        nbPtInt = TheProd.sum()
        if (nbPtInt == nbPtAfter and nbPtBefore == nbPtAfter):
            eStr=' no erase'
        else:
            eStr='';
            NumberDones = np.zeros((eta_rho, xi_rho))

        print('Iteration #', Iter)
        print('current r=', realR, '  nbPointMod=', nbPointMod, eStr)
        print(' ')

        Iter = Iter + 1

        if (IsFinished == 1):
            break
 
    return RetBathy# main routine
def smooth_bathy(grid_file,rx0max,method):
    
    id = Dataset(grid_file,'a')
    h = id.variables['h'][:,:]
    mask_rho = id.variables['mask_rho'][:,:]
    mask_deep = np.where(h>=2000.0,1,0)
    mask_box = np.zeros(np.shape(h))
    mask_box[:21,:]=1.0
    mask_box[509:,:]=1.0
    mask_box[:,:21]=1.0
    mask_box[:,609:]=1.0
   # mask_box[11:16,:]=1.0
   # mask_box[514:519,:]=1.0
   # mask_box[:,11:16]=1.0
   # mask_box[:,614:619]=1.0
    pn = id.variables['pn'][:,:]
    pm = id.variables['pm'][:,:]
    
    areaMatrix = 1.0/(pn*pm)
    
    if method=='positive':
        h_smooth = smoothing_Positive_rx0(mask_box,h,rx0max)
    elif method=='negative':
        h_smooth = smoothing_Negative_rx0(mask_rho,h,rx0max)
    elif method=='positiveVolume':
        h_smooth = smoothing_PositiveVolume_rx0(mask_deep,h,rx0max,areaMatrix)
    elif method=='negativeVolume':
        h_smooth = smoothing_NegativeVolume_rx0(mask_off,h,rx0max,areaMatrix)
    elif method=='plusMinus':
        h_smooth = smoothing_PlusMinus_rx0(mask_box,h,rx0max,areaMatrix)
    elif method=='laplace':
        h_smooth = smoothing_Laplacian_rx0(maks_rho,h,rx0max)
        
    id.variables['h'][:,:]= h_smooth
    id.close()# Command-line interface
if __name__ == "__main__":

    grid_file='/home/ubuntu/bigStick/waom10Grids/sledge.nc'
    method ='plusMinus'
    rx0max =0.02
    
    smooth_bathy(grid_file,rx0max,method)
