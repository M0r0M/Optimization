""" Univariable Polynomial Chaos implementation """

import numpy as np # necessary to find the roots of polynomials

# ------------------------------------------------------
# hermite chaos

def HermiteCoeffProbaOld(n):
    """ DEPRECIATED """
    """ recursively compute the
	probabilist Hermite polynomial
	of order n
	i.e. the Hermite Chaos polynomial """
    if n is 0: 
        return [1]
	if n is 1: return [1,0]
	else:
		result = HermiteCoeffProba(n-1)
		result.append(0)
		i = 0
        temp = HermiteCoeffProba(n-2)
        while i <= n-2:
			result[i+2] -= (n-1) * temp[i]
			i += 1
        return result

def HermiteCoeffProbaList(n):
    """ computes the list of 
    probabilist Hermit polynomial
    up to index n via their recusive definition
    ------------------------------------------------------
   input: n, an integer
    ------------------------------------------------------
    output: a list, list of the coefficients of the n first Hermite polynomials
    ------------------------------------------------------
    """
    listHermit = list()
    for i in range(n+1):
        if i == 0:
            listHermit.append([1])
        elif i == 1:
            listHermit.append([1,0])
        elif i > 1:
            nextHermite = 0
            nextHermite = copyList(listHermit[i-1])
            nextHermite.append(0)
            for j in range(i-1):
                nextHermite[j+2] -= (i-1) * listHermit[i-2][j]
            listHermit.append(nextHermite)
    return listHermit

def HermiteCoeffProba(n):
    """renaming of the HermiteCoeffProbaList function
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, list of coefficients (decreasing order) of the polynomial chaos of order n 
    ------------------------------------------------------
    """
    return HermiteCoeffProbaList(n)[n]    

def PolyChaos(x,n):
    # !!! USES NUMPY !!!
	""" evaluate the probabilist
	Hermite polynomial of order n
	in x 
    ------------------------------------------------------
    input: n, an integer
    input: x, a float
    ------------------------------------------------------
    output: a float, the evaluation of the polynmial chaos of order n in x
    ------------------------------------------------------
    """
	return np.polyval(HermiteCoeffProba(n),x)

def Evaluate(Polynomial,x):
    # !!! USES NUMPY !!!
	""" evaluate a list of
	coefficients of a
	polynomial in x 
    ------------------------------------------------------
    input: Polynomial, a list of coefficients (in decreasing order)
    input: x, a float
    ------------------------------------------------------
    output: a float, the evaluation of the polynomial in x
    ------------------------------------------------------
    """
	return np.polyval(Polynomial,x)



def HermiteCoeff(n):
	""" recursively compute the list of coefficients 
	of the physicist Hermite polynomial
	of order n 
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, the list of Hermite Physist ploynomial (list of list of float)
    ------------------------------------------------------
    """
	if n is 0:
		return [1]
	if n is 1:
		return [2,0]
	else:
		result = list()
		for coeff in HermiteCoeff(n-1):
			result.append(2 * coeff)
		result.append(0)
		i = 0
		while i <= n-2:
			result[i+2] -= 2*(n-1) * HermiteCoeff(n-2)[i]
			i += 1
		return result

def Hermite(x,n):
    # !!! USES NUMPY !!!
	""" return the evaluation in x 
	of the Hermite polynomial of
	order n 
    ------------------------------------------------------
    input: x, a float
    input: n, an integer
    ------------------------------------------------------
    output: a float, the evaluation of the hermite polynomial of order n in x
    ------------------------------------------------------
    """
	return np.polyval(HermiteCoeff(n),x)

def RatioHn(x,n):
	"""returns the ratio of Hn/Hn' 
    ------------------------------------------------------
    input: x, a float
    input: n, an integer
    ------------------------------------------------------
    output: a float, the quotient of the hermite polynomial of order n and its derivate at x
    ------------------------------------------------------
    """
	return Hermite(x,n)/(2 * n * Hermite(x,n-1))

def RootsHn(n):
    # !!! USES NUMPY !!!
	""" Gives the list of roots 
	of the Hermite polynomial 
	of order n 
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, the roots of the hermite polynomial of order n (list of float)
    ------------------------------------------------------
    """
	return np.roots(HermiteCoeff(n))

def HermiteXiWi(n):
	""" gives a list of lists
	xi: collocation points
	wi: Hermite weights 
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, the collocation point and the weight of the Hermite quadrature of order n
    ------------------------------------------------------
    """
	result = list()
	rootList = RootsHn(n)
	for root in rootList:
		xi = root
		wi = ((2**(n-1)) * fact(n-1) * (np.pi**0.5)) / (n*Hermite(xi,n-1)**2)
		result.append([xi,wi])
	return result

def tupleHermiteXi(n,multiplier=1):
    """
    ------------------------------------------------------
    input: n, an integer
    input: multiplier, a float
    ------------------------------------------------------
    output: a tuple, the collocation point of hermite quadrature of order n multiplied by multiplier
    ------------------------------------------------------
    """
    XiWi = HermiteXiWi(n)
    Xi = list()
    for point in XiWi:
        Xi.append(multiplier*point[0])
    Xi = tuple(Xi)
    return Xi

def DeterministicSamples(Mean,stdDev,PCorder,multiplier=np.sqrt(2),XiWi=[]):
    """
    ------------------------------------------------------
    input: Mean, a float
    input: stdDev, a float
    input: PCoder, an integer
    input: multiplier, a float
    XiWi: a list
    ------------------------------------------------------
    output: a list, the list of deterministic sample inputs (list of floats)
    """
    if XiWi == []:
        XiWi = HermiteXiWi(PCorder)
    samples = list()
    for i in range(PCorder):
        samples.append(Mean + stdDev * XiWi[i][0] * multiplier)
    return samples
    
def HermiteGaussQuadrature(n,polynomial):
    """ DEPRECIATED """
    """ Compute the gauss Hermite quadrature of order n """
    """ due to approximation errors, two orthogonal PC 
	have a small gauss quadrature (~1e-14) but not 0 
    ------------------------------------------------------
    input: n, an integer
    input: polynomial, a list of float
    ------------------------------------------------------
    output: a float, the Hermite quadrature of order n of the polynomial
    """
    result = 0
    sampleList = HermiteXiWi(n)
    for sample in sampleList:
		# sample[0] = xi and sample[1] = wi
		result += sample[1] * Evaluate(polynomial,sample[0]) 
    return result
    
def HermiteGaussQuadratureCorrected(n,polynomial):
	""" Compute de the Hermite Gauss quadrature of order n """
	""" correct to quadrature the PC of order 0 to 1;
	and the other PC of higher order to 0 
    ------------------------------------------------------
    input: n, an integer
    input: polynomial, a list of float
    ------------------------------------------------------
    output: a float, the Hermite quadrature of order n of the polynomial
    """
	result = 0
	sampleList = HermiteXiWi(n)
	for sample in sampleList:
		# sample[0] = xi and sample[1] = wi
		result += sample[1] * Evaluate(polynomial,sample[0]*np.sqrt(2)) 
	result = result / np.sqrt(np.pi)
	
	return result

# ------------------------------------------------------
# legendre chaos

def LegendreCoeffList(n):
    """ compute the list of
    Legendre polynomials up to
    index n 
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, list of the n first Legendre polynomials (list of list of float)
    ------------------------------------------------------
    """
    listLegendre = list()
    for i in range(n+1):
        if i == 0: # L0 = 1
            listLegendre.append([1.0])
        elif i == 1: # L1 = X
            listLegendre.append([1.0,0.0])
        elif i > 1: # n Ln = (2*n-1)XL(n-1) - (n-1)L(n-2)
            nextLegendre = copyList(listLegendre[i-1])
            for j in range(i):
                nextLegendre[j] *= 2*i-1
            nextLegendre.append(0)
            for j in range(i-1):
                nextLegendre[j+2] -= (i-1)*listLegendre[i-2][j]
            for j in range(i+1):
                nextLegendre[j] = nextLegendre[j] / i
            listLegendre.append(nextLegendre)
    return listLegendre

def LegendreCoeff(n):
    """
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, list of the legendre polynomial of order n (list of float)
    ------------------------------------------------------
    """
    return LegendreCoeffList(n)[n]
            
def Legendre(x,n):
    # !!! USES NUMPY !!!
    """ return the evaluation in x
    of Legendre polynomial of order n 
    ------------------------------------------------------
    input: n, an integer
    input: x, a float
    ------------------------------------------------------
    output: a float, the evaluation of the Legendre polynomial of order n in x
    ------------------------------------------------------
    """
    return np.polyval(LegendreCoeff(n),x)      

def DeriveLegendreCoeff(n):
    """ return the list of coefficients
    of the derivate of the Legendre polynomial 
    of order n 
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, the list of coefficients of the derviate of the Legendre polynomial of order n
    ------------------------------------------------------
    """
    legendre = LegendreCoeff(n)
    legendre.reverse()
    for i in range(n+1):
        legendre[i] *= i
    del legendre[0]
    legendre.reverse()
    return legendre

def DeriveLegendre(x,n):
    # !!! USES NUMPY !!!
    """ return the evaluation in x
    of the derivate of the Legendre
    polynomial of order n 
    ------------------------------------------------------
    input: n, an integer
    input: x, a float
    ------------------------------------------------------
    output: a float, the evaluation of the derivate of Legendre polynomial of order n
    ------------------------------------------------------
    """
    return np.polyval(DeriveLegendreCoeff(n),x)

def RootsLn(n):
    # !!! USES NUMPY !!!
    """ return the list of roots of 
    the Legendre polynomial of degree n 
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, list of the roots of the Legendre polynomial of order n
    ------------------------------------------------------
    """
    return np.roots(LegendreCoeff(n))
    
def LegendreXiWi(n):
    """ gives a list of lists 
    xi: collocation points 
    wi: Legendre weights 
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, the list of the couples of collocation point and weight [xi,wi] (a list of list)
    ------------------------------------------------------
    """
    result = list()
    rootList = RootsLn(n)
    for root in rootList:
        xi = root
        wi = 2 / ((1-xi**2) * DeriveLegendre(xi,n)**2)
        result.append([xi,wi])
    return result

def LegendreGaussQuadrature(n,polynomial,offset=0,scale=1):
    """ compute the Legendre Gauss quadrature of order n 
    in an interval [-1,1] 
    ------------------------------------------------------
    input: n, an integer
    input: polynomial, a list of floats
    input: offset, a float
    input: scale, a float
    ------------------------------------------------------
    output: a float, the gauss quadrature of order n of the polynomial 
    ------------------------------------------------------
    """
    result = 0
    sampleList = LegendreXiWi(n)
    for sample in sampleList:
        # sample[0] = xi and sample[1] = wi
        result += sample[1] * Evaluate(polynomial,scale*sample[0]+offset)
    result = result * scale
    
    return result

def LegendreGaussQuadratureInterval(n,polynomial,a=-1,b=1):
    """ compute the Legendre Gauss quadrature of order n 
    in a interval [a,b] 
    ------------------------------------------------------
    input: n, an integer
    input: polynomial, a list of floats
    input: a, a float
    input: b, a float
    ------------------------------------------------------
    output: a float, the gauss quadrature of order n of the polynomial between a and b
    ------------------------------------------------------
    """
    offset = (a+b)/2.0
    scale = (b-a)/2.0
    return LegendreGaussQuadrature(n,polynomial,offset,scale)
    
# ------------------------------------------------------
# laguerre guass quadrature

def LaguerreCoeffList(n):
    """ compute the list of
    Legendre polynomials up to
    index n
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, list of the n first Laguerre polynomials (list of list of float)
    ------------------------------------------------------
    """
    listLaguerre = list()
    for i in range(n+1):
        if i == 0: # La0 = 1
            listLaguerre.append([1.0])
        elif i == 1: # La1 = 1 - x
            listLaguerre.append([-1.,1.])
        elif i > 1: # k Lak = -xLa(k-1) + (2k-1)*La(k-1) - (k-1)*La(k-1)
            nextLaguerre = copyList(listLaguerre[i-1])
            nextLaguerre.append(0)
            for j in range(i+1):
                nextLaguerre[j] = - nextLaguerre[j]
            for j in range(i):
                nextLaguerre[j+1] += (2*i-1)*listLaguerre[i-1][j]
            for j in range(i-1):
                nextLaguerre[j+2] -= (i-1)*listLaguerre[i-2][j]
            for j in range(i+1):
                nextLaguerre[j] = nextLaguerre[j] / i
            listLaguerre.append(nextLaguerre)
    return listLaguerre

def LaguerreCoeff(n):
    """
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, list of the coefficients of the Laguerre polynomial of order n (list of float)
    ------------------------------------------------------
    """
    return LaguerreCoeffList(n)[n]

def Laguerre(x,n):
    # !!! USES NUMPY !!!
    """
    ------------------------------------------------------
    input: x, a float
    input: n, an integer
    ------------------------------------------------------
    output: a float, the evaluation of the Laguerre polynomial of order n in x
    ------------------------------------------------------
    """
    return np.polyval(LaguerreCoeff(n),x)
    
def DeriveLaguerreCoeff(n):
    """ return the list of coefficients
    of the derivate of the Legendre polynomial 
    of order n 
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, the list of coefficients of the derivate of the Laguerre polynomial of order n (list of float)
    ------------------------------------------------------
    """
    laguerre = LaguerreCoeff(n)
    laguerre.reverse()
    for i in range(n+1):
        laguerre[i] *= i
    del laguerre[0]
    laguerre.reverse()
    return laguerre

def DeriveLaguerre(x,n):
    # !!! USES NUMPY !!!
    """
    ------------------------------------------------------
    input: x, a float
    input: n, an integer
    ------------------------------------------------------
    output: a float, the derivate of the Lagurre polynomial of order n
    ------------------------------------------------------
    """
    return np.polyval(DeriveLaguerreCoeff(n),x)
    
def RootsLan(n):
    # !!! USES NUMPY !!!
    """
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, the list of roots of the Laguerre polynomial of order n (list of float)
    ------------------------------------------------------
    """
    return np.roots(LaguerreCoeff(n))
    
def LaguerreXiWi(n):
    """
    ------------------------------------------------------
    input: n, an integer
    ------------------------------------------------------
    output: a list, the list of the couples collocation points and weights of the Laguerre quadrature (list of list of float)
    ------------------------------------------------------
    """
    result = list()
    rootList = RootsLan(n)
    for root in rootList:
        xi = root
        wi = 1 /(root*DeriveLaguerre(root,n+1)**2)
        result.append([xi,wi])
    return result

def LaguerreGaussQuadrature(n,polynomial):
    """
    ------------------------------------------------------
    input: n, an integer
    input: polynomial, a list of floats
    ------------------------------------------------------
    output: a float, the Laguerre quadrature of order n of the polynomial
    ------------------------------------------------------
    """
    result = 0
    sampleList = LaguerreXiWi(n)
    for sample in sampleList:
        result += sample[1] * Evaluate(polynomial,sample[0])
    return result


# ------------------------------------------------------
# Non-Intrusive PC (by S. Hosder), To be tested

def NIPC(n,listOfSamples,listOfResults,scheme='Hermite'): 
    """ it is assumed that len(listOfSamples) == len(listOfResults) == n """
    
    if scheme == 'Legendre':
        CoeffList = LegendreCoeffList(n)
    elif scheme == 'Laguerre':
        CoeffList = LaguerreCoeffList(n)
    else:
        CoeffList = HermiteCoeffList(n)
    matrixCoeff = list()
    for i in range(n+1):
        matrixLine = list()
        for PCexpansion in CoeffList:
            matrixLine.append(Evaluate(PCexpansion,listOfSamples[i]))
        matrixCoeff.append(matrixLine)
    PhiMatrix = np.matrix(matrixCoeff) # create the proper matrix object
    try:
        InvPhiMatrix = PhiMatrix.getI() # inverse the matrix
    except np.linalg.linalg.LinAlgError:
        print "Matrix is not inversible"
        # do something here
        return None
    MatrixResult = np.matrix(listOfResults).getT()
    MatrixOfAlpha = InvPhiMatrix.dot(MatrixResult).getT().getA()[0]
    listOfAlpha = list()
    for i in range(len(MatrixOfAlpha)): # transform the matrix into a list
        listOfAlpha.append(MatrixOfAlpha[i])
    return listOfAlpha

# ------------------------------------------------------
# auxiliary functions

def fact(n):
	""" recursively compute the 
	factorial of n """
	if n < 2:
		return 1
	else:
		return n * fact(n-1)

def max(a,b):
    """ DEPRECIATED """
    """ already implemented in python """
    """ return the max of
	a and b """
    if a > b: 
		return a
    else: 
		return b

def add(P,Q): # add two polynomials
	""" add two lists of coefficients
	representing polynomials """
	result = list()
	if len(P) > len(Q):
		nmin = len(Q)
		nmax = len(P)
		polyMax = list(P)
		polyMin = list(Q)
	else:
		nmin = len(P)
		nmax = len(Q)
		polyMax = list(Q)
		polyMin = list(P)
	polyMin.reverse()
	polyMax.reverse()
	for i in range(nmin):
		result.append(polyMax[i]+polyMin[i])
	for i in range(nmin,nmax):
		result.append(polyMax[i])
	result.reverse()
	return result

def multiply(P,Q): # multiply two polynomials
	""" multiply two lists of coefficients
	representing polynomials """
	result = list()
	p = len(P)
	q = len(Q)
	for k in range(1,p+q):
		result.append(0)
	for i in range(p):
		for j in range(q):
			result[i+j] += P[i]*Q[j]
	return result

def copyList(nList):
    """ DEPRECIATED """
    """ use the moduel 'copy.deepcopy' instead """
    """ recursive function, 
	return a copy of a list """
    if nList == []:
		return []
    if type(nList[0]) is list:
		result = list()
		for subList in nList:
			result.append(copyList(subList))
		return result
    else:
		return list(nList)
        
def product(*scalars):
    """ return the product of the scalars """
    result = 1
    for scalar in scalars:
        result *= scalar
    return result

def increment(aList,listOfValues=list()):
    """ increment recursivly a list with predefined values """
    if listOfValues == []:
        for i in range(len(aList)):
            listOfValues.append(i)
    max = listOfValues[len(listOfValues)-1]
    head = aList.pop(0)
    tail = aList
    isTailMaxed = True
    for element in tail:
        if element != max:
            isTailMaxed = False
    if isTailMaxed:
        newTail = list()
        for element in tail:
            newTail.append(listOfValues[0])
        tail = newTail
        currentPosition = listOfValues.index(head) # can return ValueError if head not in listOfValues
        if head != max:
            head = listOfValues[currentPosition+1] # can return IndexError if head is the max element
    else:
        tail = increment(tail,listOfValues)
    aList = list(tail)
    aList.insert(0,head)
    return aList

def Mean(listCoeff,listPC):
    return listCoeff[0]

def StandardDeviation(listCoeff,listPC):
    variance = 0
    for i in range(1,len(listPC)):
        variance += (listCoeff[i]**2) * fact(i) # fact(i)==GHQ(listPC[i]**2)
    return np.sqrt(variance)

def PrintXiWi(n=5,write=False):
    result = list()
    for i in range(n+1):
        result.append(HermiteXiWi(i))
    content = str()
    if not write:
        content+='Hermite Xi Wi\n'
    for i in range(1,len(result)):
        content+=str(i)+'\n'
        for j in range(len(result[i])):
            content+=str(result[i][j][0]) + '\t' + str(result[i][j][1]) + '\n'
    if write:
        file = open('HermiteXiWi.txt','w')
        file.write(content)
        file.close()
    
    result = list()
    for i in range(n+1):
        result.append(LaguerreXiWi(i))
    if write:
        content = str()
    else:
        content+='Laguerre Xi Wi\n'
    for i in range(1,len(result)):
        content+=str(i)+'\n'
        for j in range(len(result[i])):
            content+=str(result[i][j][0]) + '\t' + str(result[i][j][1]) + '\n'
    if write:
        file = open('LaguerreXiWi.txt','w')
        file.write(content)
        file.close()
    
    result = list()
    for i in range(n+1):
        result.append(LegendreXiWi(i))
    if write:
        content = str()
    else:
        content+='Legendre Xi Wi\n'
    for i in range(1,len(result)):
        content+=str(i)+'\n'
        for j in range(len(result[i])):
            content+=str(result[i][j][0]) + '\t' + str(result[i][j][1]) + '\n'
    if write:
        file = open('LegendreXiWi.txt','w')
        file.write(content)
        file.close()
    else:
        print content

if __name__ == '__main__':
    PrintXiWi()