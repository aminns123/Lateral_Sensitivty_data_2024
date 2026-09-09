"""Selected author helper functions. See config/function_provenance.json."""

import numpy as np

import os, itertools, math, statistics

from itertools import combinations, permutations, product

from scipy.optimize import curve_fit

from scipy.stats import norm, skewnorm, binom, gaussian_kde, mode

def indexInListRemove(xAxis,yAxis,indexList):
    
    newX, newY = [], []
    for j in range(len(yAxis)):
        if j in indexList:
            pass
        elif j not in indexList: 
            newX.append(xAxis[j])
            newY.append(yAxis[j])
        
    return newX, newY

def findCorrespondingIndices(valueList, findListValues):
    # return index
    xfind    = []
    for j in range(len(findListValues)):
        findValue = findListValues[j]
        for i in range(len(valueList)):
            if (valueList[i] == findValue) and (findValue not in xfind):
                xfind.append(i)
            else:
                pass

    return xfind

def getIndexedValues(listIndexs, listAxis):
    """
        Return indexed elements from "listAxis" list.
    """
    outList = []
    for i in listIndexs:
        outList.append(listAxis[int(i)])
    return outList

def count_reversals_HighLow(contrastList):
    averageList = []
    TrialNumber = []
    countR      = 0
    checkSign   = 0
    signValue   = -1
    
    j           = 0
    
    if len(contrastList) > 1:
        if contrastList[j+1]-contrastList[j] < 0: # i.e. N12 - N10, then get  j-1
            signValue = -1
        elif contrastList[j+1]-contrastList[j] > 0: 
            signValue = 1
        elif contrastList[j+1]-contrastList[j] == 0: 
            pass
            
        checkSign = signValue
        j         = 0
        while j < len(contrastList)-2:
            j+=1
            if contrastList[j+1]-contrastList[j] < 0: # i.e. N12 - N10, then get  j-1
                signValue = -1
            elif contrastList[j+1]-contrastList[j] > 0: 
                signValue = 1
            elif contrastList[j+1]-contrastList[j] == 0: 
                pass
            
            if checkSign != signValue:
                checkSign = signValue 
                countR   += 1
                averageList.append(contrastList[j])
                TrialNumber.append(j)
            elif checkSign == signValue:
                pass
    elif len(contrastList) <= 1:
        TrialNumber=[1]
        TrialNumber=[0]
        averageList=[0] 
        
    return len(TrialNumber), TrialNumber, averageList 

def makeList(xArray, func, **kwargs):
    newList = []
    for x in xArray:
        newList.append(func(x, **kwargs))
    return newList

def greenFUNC_phase(x, x0, x1, x2, x3): # , x7, x8, x9, x10, x11):
    A,phase, kn, lam = x0, x1, x2, x3
    xf         = -0.5
    X          = x-xf
    value      = A*(np.exp(-lam*abs(X))*(np.cos((X)*kn+(phase))-np.sign(X)*np.sin((X)*kn+(phase))))
    return value

def find_peak_vally_trio(peaks, vallys, data):
    trioList = [[] for x in range(len(peaks))]
    vPrev    = min(data)
    vallysCOPY = sorted(vallys.copy())
    
    # from left to right
    for i, mode in enumerate(sorted(peaks)):
        P1 = mode
        V1 = vallysCOPY[i]
        
        if (i+1) < len(vallysCOPY):
            V2 = vallysCOPY[i+1]
        else:
            V2 = max(data)
        
        if (V1 < P1):
            trioList[i].extend([V1,P1,V2])
        elif (V1 > P1):
            trioList[i].extend([vPrev,P1,V1])
        vPrev = V1    
    return trioList

def largestNvalues(listValues, Nreturn):
    maxList = []
    listAbs = []
    for i in listValues:
        listAbs.append(abs(i))
    
    copyList=sorted(listAbs)
    for j in range(len(copyList)-1,len(copyList)-Nreturn-1,-1):
        maxList.append(copyList[j])
    
    outIndicies = findCorrespondingIndices(listAbs, maxList)    
    return outIndicies

def goodnessFit_Residual_Analysis(x_original,y_original,Nstd,fitFunction, **fitParams):
    from scipy.optimize import curve_fit
    X = x_original
    Y = np.array(y_original)

    # Compute residuals
    residuals = []
    for j in range(len(X)):
        value = Y[j] - fitFunction(X[j], **fitParams)
        residuals.append(value)
  
    # Identify outliers using standard deviation method
    threshold = Nstd * np.std(residuals)
    outliers  = Y[np.abs(residuals) > threshold]
    npres     = np.array(residuals)
    resOut    = npres[np.abs(residuals) > threshold]
    
    outlier_indices = np.where(np.abs(residuals) > threshold)[0]
    xIndex=[]
    for j in outlier_indices.tolist():
        xIndex.append(x_original[j])
    
    #findCorrespondingIndex(outliers)
    return xIndex,resOut

def nLargest_elements(listE, nLarge):
    """
        Finf the n many largest valued elements in list in decending order.
    """
    listCopy = []
    for j in range(len(listE)-1):
        ki = listE[j]
        kj = listE[j+1]
        if ki>kj:
            if ki not in listCopy:
                listCopy.append(ki)
        elif ki<kj:
            if kj not in listCopy:
                listCopy.append(kj)
        else:
            pass
    
    outList = [max(listCopy)]
    
    for i in range(nLarge-1):
        if len(listCopy)>0:
            outList.append(max(listCopy))
            listCopy.remove(max(listCopy))
        elif len(listCopy)==0:
            pass
    return outList

def get_peaks_vally_distributions(list_histogram, nPeaks, bandwidthScale=1.0, nBins='auto', kindInterp='linear', bw_method='scott', lam=0.02):
    from matplotlib import pyplot as plt
    from scipy import interpolate
    from scipy.stats import gaussian_kde
    from scipy.signal import find_peaks
    #import seaborn as sns 
    from scipy.interpolate import make_smoothing_spline
    import statistics
    
    rows, cols      = 1, 1
    fig, axs        = plt.subplots(rows, cols,figsize=(6, 6), sharex=False,sharey=False)  

    ## ==================================== ## 
    data = np.array([v for v in list_histogram if v < 10])
    #data                                = np.array(list_histogram)
    ### Use KDE estimate to find modes
    x_grid  = np.linspace(1, 8, 1000)
    kde     = gaussian_kde(data, bw_method=bw_method) #  silverman or scott
    bw      = kde.factor
    kde_bw  = gaussian_kde(data, bw_method=bw*bandwidthScale)
    y_kde   = kde_bw(x_grid)

    
    # xRange, y_fitted_interp = x_grid, y_kde #
    
    
    hist_values, nBins, _               = axs.hist(data, bins=nBins, alpha=0.0, color='darkgreen', edgecolor='dodgerblue', density=True)
    bin_centers                         = (nBins[:-1] + nBins[1:]) / 2 
    for _ in range(1):
        hist_values                         = np.append(hist_values,np.array([0.0]))
        bin_centers                         = np.append(bin_centers,np.array([bin_centers[-1]+bin_centers[-1]]))
    hist_values                         = np.insert(hist_values, 0, 0)
    bin_centers                         = np.insert(bin_centers, 0, 0)
    

    interp_func                         = interpolate.interp1d(bin_centers, hist_values, kind=kindInterp) # cubic, quadratic, linear
    xRange                              = np.linspace(0, max(bin_centers), 1000) # 0.2
    y_counts                            = interp_func(xRange)
    y_fitted_interp                     = y_counts.copy()
    
    
    spl             = make_smoothing_spline(np.array(bin_centers), np.array(hist_values), lam=lam)
    y_fitted_interp = spl(xRange)
    #print('y_fitted_interp', y_fitted_interp)
    #plt.plot(xRange, y_fitted_interp)
    #plt.show()
    ## ============== seaborn END ==============
    #axs1.plot(xRange,y_fitted_interp, c='red')
    _, rev_Index, rev_BooT_yaxis        = count_reversals_HighLow(y_fitted_interp)
    x_list                              = getIndexedValues(rev_Index, xRange)
    
    try:
        fitted_c                            = nLargest_elements(rev_BooT_yaxis, nPeaks)
        indexList                           = findCorrespondingIndices(rev_BooT_yaxis, fitted_c)
        peaks                               = getIndexedValues(indexList, x_list)
        
        indexList.pop(0)
        peaks.pop(0)
        

        indexCopy  = sorted(indexList.copy())
        indexVally = []
        if indexCopy[0]-1 >= 0:
            indexVally.append(indexCopy[0]-1)
        elif indexCopy[0]-1 <0:
            pass
        for j in range(len(indexCopy)):
            indexVally.append(indexCopy[j]+1)
        
        if max(indexVally) >= len(x_list):
            indexVally.pop(len(indexVally)-1)
        
        #fitted_c                            = nSmallest_elements(rev_BooT_yaxis, nPeaks)
        #indexList                           = findCorrespondingIndices(rev_BooT_yaxis, fitted_c)
        vallyes                             = getIndexedValues(indexVally, x_list)     
    except ValueError:
        peaks = [statistics.mode(list_histogram)]
        vallyes = [0]
 
        
    plt.close()
    
    if max(vallyes) < max(peaks):
        vallyes.append(max(data))
    elif min(vallyes) > min(peaks):
        vallyes.append(min(data))
    else:
        pass
    
    #for f in interp_f:
    #    axs1.axvline(f, c='dodgerblue', linestyle='dashed', label='f:'+str(round(f,2)))
    ## ==================================== ## 
    #peaks, properties   = find_peaks(data, prominence=(None, 0.6))
    #vallyes, properties = find_peaks(np.array(data)*-1, prominence=(None, 0.6))

    return peaks, vallyes

def frequency_match_score(X_original,Y_original,X_fit, Y_fit):
    from scipy.fft import rfft, rfftfreq
    from scipy import interpolate
    # Interpolate fitted data to original x points
    interp_func     = interpolate.interp1d(X_fit, Y_fit)#, kind='linear')
    y_fitted_interp = interp_func(X_original)
    
    # Compute FFT magnitudes
    freq_data = np.abs(rfft(Y_original))
    freq_fit  = np.abs(rfft(y_fitted_interp))
    
    # Compute frequency axis (only needed if plotting magnitude of found frequencies and their SF)
    # freqs = rfftfreq(len(X_original), d=(X_original[1] - X_original[0]))

    # Normalize
    freq_data /= np.max(freq_data)
    freq_fit  /= np.max(freq_fit)

    # Compute similarity score (cosine similarity or correlation)
    similarity = np.dot(freq_data, freq_fit) / (np.linalg.norm(freq_data) * np.linalg.norm(freq_fit))
    
    return similarity  # Closer to 1 means better frequency match

def getIndex_sortDataRegion(data,a,b):
    indexList   = []
    for i in range(len(data)): 
        if (data[i]>= a) and (data[i]<= b):   
            indexList.append(i)
        elif (data[i]< a) or (data[i]> b):   
            pass
    return indexList

def fit_to_ISF(xNew, yNew, fitFunction,guessParams, boundDW,boundUP, n_fittingTries, fitRULE):
    try:
        xy   = fitFunctionLimit([-0.5,0.5], xNew, yNew,fitFunction, guessParams,boundDW,boundUP, fitRULE, n_fittingTries)
        xFit        = xy[0]
        yFit        = xy[1]
        paramsFIT   = xy[2]
        ADM_kn_fit  = np.round(xy[2]['x2'],2)/(2*np.pi)
        ADM_lam_fit = np.round(xy[2]['x3'],2)
    except RuntimeError or RuntimeWarning:
        pass
    else:
        pass
    
    return xFit,yFit,paramsFIT,ADM_kn_fit

def fitFunctionLimit(xLimMax,xdata, ydata,function_fit,guessParams, boundDW, boundUP, fitm, fitNumber=10_000):
    from scipy.optimize import curve_fit
    
    maxX, minX      =  xLimMax[1], xLimMax[0]
    xStep           = (maxX - minX)/1000
    xdataContinuous = list(np.arange(minX, maxX,xStep) )
    params_totalFit = {}
    #try:  # lm, trf, dogbox 
    # fitm = 'lm'
    Rss, Tss, Rsquared= 0,0,0     
    if len(guessParams) > 1 and len(boundUP) > 1: 
        popt, pcov    = curve_fit(function_fit, xdata, ydata,p0=guessParams,method=fitm, maxfev=fitNumber, bounds=(boundDW, boundUP)) ## 20_000
    elif len(guessParams) <=1 and len(boundUP) > 1:
        popt, pcov    = curve_fit(function_fit, xdata, ydata,method=fitm, maxfev=fitNumber, bounds=(boundDW, boundUP)) ## 20_000
    elif len(guessParams) <=1 and len(boundUP) <= 1:
        popt, pcov    = curve_fit(function_fit, xdata, ydata,method=fitm, maxfev=fitNumber) 
    #nan_policy='omit'    
    for indexKey in range(len(popt)):
        params_totalFit.update({'x'+str(indexKey):popt[indexKey]})
        
    fitfunction_array = makeList(xdataContinuous, function_fit, **params_totalFit)


    return xdataContinuous, fitfunction_array, params_totalFit, pcov

def neg_loglik(params, model, x, y, sigma_vec):
    *theta, = params   # parameters of the model
    y_pred = model(x, *theta)
    resid = y - y_pred

    # --- Clip residuals and sigma to safe ranges ---
    resid_safe = np.clip(resid, -1e6, 1e6)
    sigma_safe = np.clip(sigma_vec, 1e-3, None)

    # Gaussian log-likelihood
    ll = -0.5 * np.sum(
        np.log(2*np.pi*sigma_safe**2) + (resid_safe**2) / (sigma_safe**2)
    )
    return -ll

def MLE_range_frequency(model, xData, yData, thetaDict, sigma, bounds, range_f,method='L-BFGS-B', fitNumber=10_000, idx='f',jac=None):
    from scipy.optimize import minimize
    import numpy as np
    from itertools import product

    
    # Convert input dict to clean, immutable array
    theta_master = np.array(list(thetaDict.values()), dtype=float)

    # Precompute mapping of parameter name → index
    idx_map = {key: i for i, key in enumerate(thetaDict)}
    #idx_f   = idx_map[idx]
    
    if isinstance(idx, list):
        idx_f = [idx_map[k] for k in idx]
    else:
        idx_f = [idx_map[idx]]

    nll_list = []
    #res_list = []

    f_combinations = list(product(range_f, repeat=len(idx_f)))

    #for f in range_f:
    for f_vals in f_combinations:
        # --- IMPORTANT: do NOT modify theta_master ---
        theta_local = theta_master.copy()
        #theta_local[idx_f] = f  # safe
        for i, f in zip(idx_f, f_vals):
            theta_local[i] = f
        
        if jac is not None:
            res = minimize(
                neg_loglik,
                theta_local,
                args=(model, xData, yData, sigma),
                method=method,
                bounds=bounds,
                jac=jac, 
                options = ({
                #'maxiter': 200,
                'maxfev': fitNumber,
                "maxiter": 15,
                "xtol": 1e-3,
                "ftol": 1e-3,
                })
            )
        else:
            res = minimize(
                neg_loglik,
                theta_local,
                args=(model, xData, yData, sigma),
                method=method,
                bounds=bounds,
                options = ({
                #'maxiter': 200,
                'maxfev': fitNumber,
                "maxiter": 15,
                "xtol": 1e-3,
                "ftol": 1e-3,
                })
            )
        # If minimizer fails or returns NaN, treat as bad
        if (not res.success) or np.isnan(res.fun):
            nll_list.append(np.inf)
        else:
            nll_list.append(res.fun)

    # Pick the best frequency (smallest NLL)
    best_idx = np.argmin(nll_list)
    #best_f   = range_f[best_idx]

    # --- Second pass: refine from best frequency ---
    theta_start = theta_master.copy()
    #theta_start[idx_f] = best_f
    
    best_f_vals    = f_combinations[best_idx]
    for i, f in zip(idx_f, best_f_vals):
        theta_start[i] = f

    #print('theta_start', theta_start)
    if jac is not None:
        res_final = minimize(
            neg_loglik,
            theta_start,
            args=(model, xData, yData, sigma),
            method=method,
            bounds=bounds,
            jac=jac, 
            options = ({
                #'maxiter': fitNumber
                "maxiter": 40,
                "xtol": 1e-6,
                "ftol": 1e-6,
                })
        )
    else:
        res_final = minimize(
            neg_loglik,
            theta_start,
            args=(model, xData, yData, sigma),
            method=method,
            bounds=bounds,
            options = ({
                #'maxiter': fitNumber
                "maxiter": 40,
                "xtol": 1e-6,
                "ftol": 1e-6,
                })
        )
    return res_final

def fitFunction(xdata, ydata,function_fit,guessParams, boundDW, boundUP, fitm):
    from scipy.optimize import curve_fit
    
    maxX, minX      = max(xdata), min(xdata)
    xStep           = (maxX - minX)/1000
    xdataContinuous = list(np.arange(minX, maxX,xStep) )
    params_totalFit = {}
    #try:  # lm, trf, dogbox 
    # fitm = 'lm'
    Rss, Tss, Rsquared= 0,0,0     
    if len(guessParams) > 1 and len(boundUP) > 1: 
        popt, pcov    = curve_fit(function_fit, xdata, ydata,p0=guessParams,method=fitm,max_nfev=5000, maxfev=5000, bounds=(boundDW, boundUP)) ## 20_000, 100_000
    elif len(guessParams) <=1 and len(boundUP) > 1:
        popt, pcov    = curve_fit(function_fit, xdata, ydata,method=fitm,max_nfev=5000, maxfev=5000, bounds=(boundDW, boundUP)) ## 20_000
    elif len(guessParams) <=1 and len(boundUP) <= 1:
        popt, pcov    = curve_fit(function_fit, xdata, ydata,method=fitm,max_nfev=5000, maxfev=5000) 
    #nan_policy='omit'    
    for indexKey in range(len(popt)):
        params_totalFit.update({'x'+str(indexKey):popt[indexKey]})
        
    fitfunction_array = makeList(xdataContinuous, function_fit, **params_totalFit)


    return xdataContinuous, fitfunction_array, params_totalFit
