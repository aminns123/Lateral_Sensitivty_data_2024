"""Selected author helper functions. See config/function_provenance.json."""

import numpy as np

import os, itertools, math, statistics

from itertools import combinations, permutations, product

from scipy.optimize import curve_fit

from scipy.stats import norm, skewnorm, binom, gaussian_kde, mode

def readText_toList(path):
    file        = open(path, "r")
    contents    = file.read()
    file.close()
    contents = contents.replace("\t",",")
    contLists= contents.split('\n')


    ListRows  = []
    for j in range(len(contLists)):
        contLists_j = contLists[j].replace("]","")
        contLists_j = contLists_j.split(",")
        ListRows.append([float(elements) for elements in contLists_j])
    
    outColumns = [[] for i in range(len(ListRows[0]))]
    
    for j in range(len(ListRows)):
        for i in range(len(ListRows[0])):
            outColumns[i].append(ListRows[j][i]) 
    
    return outColumns

def PositionCondition(valueIN, condition1, condition2, condition3, condition4):
    if condition4 == False:
        if (condition1 == True) and (condition2 == False) and (condition3 == False):
            value = abs(valueIN)
        elif (condition1 == False) and (condition2 == True) and (condition3 == False):
            value = -abs(valueIN)
        elif (condition1 == False) and (condition2 == False) and (condition3 == True):
            value = (valueIN)
    elif condition4 == True:
        value = abs(valueIN)
    return value 

def accendingOrder( xdata: [],ydata: [], sortAgainst_x: []): 
    # match xdata with list with, in accending order , small to large. 
    sortAgainst = sorted(sortAgainst_x)
    sortedY = []
    sortedX = []

    for i in range(len(sortAgainst_x)):
        for j in range(len(xdata)):
            
            if xdata[j] == sortAgainst[i]:
                sortedY.append(ydata[j])
                sortedX.append(xdata[j])
            elif xdata[j] != sortAgainst[i]:
                continue 
    return sortedX,sortedY

def condition1(lisT, n):
	for j in range(len(lisT)-1, 0, -1):
		if lisT[j] >= n:
			lisT[j-1] += 1
			if lisT[j-1]+1 < n:  # <=
				lisT[j] = lisT[j-1]+1
			elif lisT[j-1]+1 >= n:  # >
				lisT[j] = lisT[j]
		elif lisT[j] != n:
			pass
	return lisT

def condition2(lisT, n):
	j = len(lisT)-1
	found = False
	if lisT[len(lisT)-1] >= n-1:  # ==, >= n
		while found == False:
			if lisT[j] <= lisT[j-1]+1:
				pass
			elif lisT[j] > lisT[j-1]+1:
				lisT[j-1] += 1
				lisT[j] = lisT[j-1]+1
				found = True
			j += -1
	elif lisT[len(lisT)-1] != n-1:  # n
		pass
	return lisT

def condition3(lisT, n):
	"""
	 If any element is greater than n, then reset value.
	"""
	for x in range(1, len(lisT), 1):
		if lisT[x] >= n-1:  # >= n
			lisT[x] = lisT[x-1]+1
		elif lisT[x] < n-1:  # < n
			pass
	return lisT

def endCondition(lisT, n):
	count = 1

	if lisT[len(lisT)-1] >= n-1:
		count = 0
		for i in range(len(lisT)-1, 0, -1):
			if (lisT[i] == lisT[i-1]+1):
				count += 0
			elif (lisT[i] != lisT[i-1]+1):
				count += 1
	elif lisT[len(lisT)-1] != n-1:
		pass

	return count

def removeErrorIndex(listIndex, n):
    valuePops = []
    for j in range(len(listIndex)):
        # print(listIndex[j])
        if n in listIndex[j]:
            # print('pop:',listIndex[j][:])
            valuePops.append(listIndex[j][:])
        elif n not in listIndex[j]:
            pass
    # print(valuePops)
    for i in range(len(valuePops)):
        listIndex.remove(valuePops[i])
    # print(len(listIndex))
    return listIndex

def conditionP1End(pIndexList, n):
    count = 1
    if pIndexList[0] == n-1:
        count = 0
    elif pIndexList[0] < n-1:
        count = 1
    return count

def returnPositionList(dataN):
        positionList        = []
        for j in range(0, len(dataN), 1): 
            for m in range(0, len(dataN[j]), 1):
                positionList.append(dataN[j][m][1])

        startIndex          = 0
        checked_pos         = check_values(positionList)
        n_types_pos, totaln = count_values_type(positionList, checked_pos, startIndex)
        n_types_pos         = count_type(n_types_pos)
        arrayPOS            = appendTrials(n_types_pos)
        return arrayPOS, n_types_pos

def listCombitorial(n, r):
    pipList = []
    # testA   = find_indexPop(n, r, pipList)
    testA   = find_indexPop_Complete(n, r, pipList)
    testA   = removeErrorIndex(testA, n)  
    return testA

def find_indexPop_Complete(n,r, popList):
    if int(n-r) == 1:
        p = n-r
        n = n+1
        pIndexList = list(np.arange(0, p, 1))
        popList.append(pIndexList[:])
        # print(pIndexList)
        pS  = len(pIndexList)
        end = False
        while end == False:
            pIndexList[pS-1] += 1
            popList.append(pIndexList[:])
            endValue = conditionP1End(pIndexList, n)
            
            if endValue == int(0):
                end = True
            elif endValue > 0:
                pass
    elif int(n-r) > 1:
        p = n-r
        n = n+1
        pIndexList = list(np.arange(0, p, 1))
        popList.append(pIndexList[:])
        # print(pIndexList)
        pS = len(pIndexList)
    
        end = False
        while end == False:
            pIndexList[pS-1] += 1
            pIndexList = condition1(pIndexList, n)
            pIndexList = condition2(pIndexList, n)
            pIndexList = condition3(pIndexList, n)
    
            # print(pIndexList)
            # popList.append(pIndexList[:])

            if n in pIndexList[:]:
                continue
            elif n not in pIndexList[:]:
                popList.append(pIndexList[:])

            endValue = endCondition(pIndexList, n)
            if endValue == int(0):
                end = True
            elif endValue > 0:
                pass
    elif int(n-r) == 0:
        p          = n-r
        pIndexList = list(np.arange(0, p, 1))
        popList.append(pIndexList[:])
        
    return popList

def findIndexFloat(arrayList, indexNow):
    indexPOSfound = 0
    for i in range(len(arrayList)):
        if (indexNow) == (arrayList[i]):
            indexPOSfound = i
            break
        elif (indexNow) != (arrayList[i]):
            pass
    return indexPOSfound

def count_type(count_type_array):
    checkedValues = []
    checkekCounts = []
    for i in range(len(count_type_array)):
        if (count_type_array[i][0]) not in checkedValues:
            checkedValues.append((count_type_array[i][0]))
            checkekCounts.append(count_type_array[i][1])
        elif (count_type_array[i][0]) in checkedValues:
            indexfound                = findIndexFloat(checkedValues, (count_type_array[i][0])) 
            checkekCounts[indexfound] = checkekCounts[indexfound] + count_type_array[i][1]
    
    returnType = []
    for i in range(len(checkedValues)):
        e1 = checkedValues[i]
        e2 = checkekCounts[i]
        returnType.append([e1, e2])
    return returnType

def count_values_type(list_of_values, value_array, startIndex):
    count_type_array = []
    total            = 0
    for z in range(len(value_array)):
        count = 0.0
        for i in range(startIndex, len(list_of_values), 1):
            if list_of_values[i] == value_array[z]:
                count += 1
            elif list_of_values[i] != value_array[z]:
                continue
        count_type_array.append([value_array[z], count])
    for index in range(len(count_type_array)):
        total += count_type_array[index][1]
    return count_type_array, total

def check_values(value_array):
    new_array = []

    for i in range(len(value_array)):
        if  value_array[i] not in new_array:
            new_array.append(value_array[i])
        elif value_array[i] in new_array:
            continue
    return new_array

def weberContrast(contrast_cpu: [], background_crt,max_crt, logFunction, logStatment:bool):
    maxC        = max_crt # max(contrast_cpu)
    contrast_W  = []
    if logStatment == False:
        for i in range(len(contrast_cpu)):
            contrast_W.append((contrast_cpu[i]-background_crt)/(maxC-background_crt))
    elif logStatment == True:
        for i in range(len(contrast_cpu)):
            contrast_W.append(logFunction((contrast_cpu[i]-background_crt)/(maxC-background_crt)))
    return contrast_W

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

def returnOnce(lisT):
    newList = []
    for i in range(len(lisT)):
        value = abs(lisT[i])
        if value not in newList:
            newList.append(value)
        elif value in newList:
            continue
    return newList

def appendTrials(trails):
    arrayT = []
    for i in range(0, len(trails), 1):
        arrayT.append(trails[i][0])
    return arrayT  

def find_ADM_data_removedFiles(folder, removeArray, TrialLength):
    lengthDir       = len(os.listdir(folder))
    files_FULL      = []
    full_contrast   = []
    file_ID         = 0
    startIndex      = 0 
    indexN          = 0

    dir_list = os.listdir(folder)

    for i in range(len(dir_list)):
        file_ID  = i
        filename = folder+"/"+dir_list[i]

        if os.path.isfile(filename) == True:
            #print(filename)
            dataMatrix          = readText_toList(filename)

            checked_ids         = check_values(dataMatrix[4])
            n_types_ids, totaln = count_values_type(dataMatrix[4], checked_ids, startIndex)
            arrayID             = appendTrials(n_types_ids)
            
            checked_pos         = check_values(dataMatrix[0])
            n_types_pos, totaln = count_values_type(dataMatrix[0], checked_pos, startIndex)
            arrayPOS            = appendTrials(n_types_pos)
            arrayPOS            = returnOnce(arrayPOS)
            
            for indexPOS in range(len(arrayPOS)):
                pos_step = arrayPOS[indexPOS]
                # ==========================
                for indexID in range(len(arrayID)):
                    L_or_R_value        = []
                    probe_position      = []
                    Huamn_value         = []
                    probe_contrast      = []
                    new_dataFull        = []
                    new_data            = []
                    absolute_Position   = []

                    id_step = arrayID[indexID]
                    indexN += 1
                    for j in range(0, len(dataMatrix[6]), 1):
                        probe_pos   = dataMatrix[0][j]
                        probe_LR    = dataMatrix[1][j]
                        human_LR    = dataMatrix[2][j]
                        probe_ct    = dataMatrix[6][j]
                        trialID     = dataMatrix[4][j]
                        
                        
                        if len(dataMatrix) >= 16:
                            flank_sf    = dataMatrix[15][j]
                            
                        elif len(dataMatrix) < 16:
                            flank_sf    = 0
                        else:
                            flank_sf    = 0
                        
                        if len(dataMatrix) >= 11:
                            correctTick = dataMatrix[10][j]
                        elif  len(dataMatrix) <= 10:
                            correctTick = 0
                            

                        if (int(id_step) == int(trialID)) and (abs(pos_step) == abs(probe_pos)):
                            
                            #file_ID  = indexN # indexN+indexPOS
                            
                            L_or_R_value.append(probe_LR)
                            probe_position.append(probe_pos)
                            Huamn_value.append(human_LR)
                            probe_contrast.append(probe_ct)
                            absolute_Position.append(abs(probe_pos))
                            new_dataFull.append([probe_LR, (probe_pos), human_LR, probe_ct, 
                                            file_ID, trialID, dir_list[i], correctTick, flank_sf])
                                            
                        elif (int(id_step) != int(trialID)) or (abs(pos_step) != abs(probe_pos)):
                            continue

                    # ==========================
                    trimmedData = removeData(new_dataFull, removeArray, TrialLength)
                    if (len(trimmedData) == 0):
                        continue
                    # ==========================
                    elif len(trimmedData) > 0:
                        files_FULL.append(trimmedData)
                        full_contrast.append(probe_contrast)
                # ==========================
        elif os.path.isfile(filename) == False:
            continue
    return files_FULL, full_contrast # absolute_Position

def permutationsBootstrapNp(num_permutations: int, value_range: int, list_size: int, sizeList: int):
    """
        Generate a list of unique random permutations (with possible repeated values in each row).

        Parameters:
        ----------
        num_permutations : int
            Total number of unique rows (permutations) to generate.
            
        value_range : int
            Values in each row will be randomly chosen in the range [0, value_range).
            i.e., up to but not including value_range.

        list_size : int
            Length of each row (number of elements per permutation).

        sizeList : int
            Final number of permutations to return (must be <= num_permutations).
            This is useful if you want to generate a large pool of unique rows but only return a subset.

        Notes:
        -----
        - A `set()` is used to ensure all rows are unique.
          Since sets only allow unique items, adding a row (as a tuple) that already exists will be ignored.
        - Each row may contain repeated values (e.g., [1, 1, 3]), but no two rows will be exactly the same.
        - If not enough unique permutations can be generated given the constraints, the function raises an error.
    """
    unique_rows = set() 
    attempts = 0
    max_attempts = 10 * num_permutations  # prevent infinite loops

    while len(unique_rows) < num_permutations and attempts < max_attempts:
        row = tuple(np.random.randint(0, value_range, size=list_size))
        unique_rows.add(row)
        attempts += 1

    if len(unique_rows) < sizeList:
        sizeList = len(unique_rows)-1
        print('actual permutation: ', sizeList)
        # raise ValueError("Could not generate enough unique rows. Try increasing value_range or list_size.")

    # Convert to NumPy array
    return np.array(list(unique_rows)[:sizeList])

def factorialN(x):
    value = 1
    for i in range(1,int(x)+1, 1):
        value = value*i
    return value

def nChooseR_list(listData: [], sample_Size):
    n = len(listData)
    r = sample_Size

    value = factorialN(n)/(factorialN(r)*(factorialN(n-r)))
    return value

def singleListpop(arrayToPop, arrayIndex):
    newList  = []
    tickList = list(np.zeros(len(arrayToPop)))
    for i in range(len(arrayIndex)):
        j            = arrayIndex[i]
        tickList[j] += 1

    for m in range(len(tickList)):
        if tickList[m] == 0:
            newList.append(arrayToPop[m])
        elif  tickList[m] > 0:
            pass
    return newList

def createCombitorialList(listData: [], sample_Size):
    Cnr         = nChooseR_list(listData, sample_Size)
    r           = sample_Size
    n           = len(listData)
    newList     = [[] for x in range(int(Cnr))]
    ArrayWhole  = [listData for x in range(int(Cnr))]
    popList     = listCombitorial(n, r)
    # e.g. n = 4, r = 2, Cnr = 6. Therefore j -> 6, i-> 2, index->4 

    for j in range(len(popList)):
        ArrayNew = ArrayWhole[j]
        popArray = popList[j]
        ArrayNew = singleListpop(ArrayNew, popArray)
        newList[j].append(ArrayNew)

    # newList = newList[0]        
    # outList = listReturn_sequence(newList, Cnr)
    return newList

def fileMethodPair(fileMethod_pair: []):
	fileId = []
	admId = []
	for i in range(len(fileMethod_pair)):
		fileId.append(fileMethod_pair[i][0])
		admId.append(fileMethod_pair[i][1])
	return fileId, admId

def checkRemoveList(valuePair: [], removeList1: [], removeList2: []):
    count = 0
    v1, v2 = int(valuePair[0]), int(valuePair[1])
    if len(removeList2) != 0:
        for i in range(len(removeList1)):
            r1 = int(removeList1[i])
            r2 = int(removeList2[i])
            if (v1 == r1) and (v2 == r2):
                count += 1
            elif (v1 != r1) and (v2 != r2):
                count += 0
    elif len(removeList2) == 0:
        for i in range(len(removeList1)):
            r1 = int(removeList1[i])
            if (v1 == r1):
                count += 1
            elif (v1 != r1):
                count += 0
    return count

def removeData(data: [], fileMethod_pair: [], TrialLength: int): # <- need to rewrite, data has to be sorted data
    human_comb      = []
    probe_comb      = []
    contrastComb    = []
    experimentID    = []
    experimentADMID = []
    fileId_remove, adm_remove = fileMethodPair(fileMethod_pair)

    newData = []
    
    for n in range(0, len(data), 1):
      contrast_p    = data[n][3]
      human_lr      = data[n][2]
      probe_pos     = data[n][1]
      probe_lr      = data[n][0]
      fileID        = data[n][4]
      methodID      = data[n][5]
      fileName      = data[n][6]
      correctTick   = data[n][7]
      flank_sf      = data[n][8]
      
      count = int(checkRemoveList([fileID, methodID], fileId_remove, adm_remove))

      if count == 0 and n < TrialLength:
        contrastComb.append(contrast_p)
        human_comb.append(human_lr)
        probe_comb.append(probe_lr)
        experimentID.append(fileID)
        experimentADMID.append(methodID)
        newData.append([probe_lr, probe_pos, human_lr, contrast_p, fileID, methodID, fileName, correctTick, flank_sf])
      elif count > 0 or n >= TrialLength:
        continue
    return newData
