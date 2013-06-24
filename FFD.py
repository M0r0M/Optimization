# ------------------------------------------------------
# importations
#import matplotlib.pyplot as plt # for ploting data, optional

# ------------------------------------------------------
# class declaration

class Point(object):
	""" 3D point """
	
	def __init__(self,x,y,z):
		self.x = x
		self.y = y
		self.z = z
		self.s = x # x coord relative to control points
		self.t = y # y coord relative to control points
		self.u = z # z coord relative to control points
	
	def compute_s(self,p0,S,T,U):
		""" compute the s coordinate of the point
		p0: origin point
		S,T,U: coordinates of the base vectors """
		vect = Vector(self.x-p0.x,self.y-p0.y,self.z-p0.z)
		TxU = vectorialProduct(T,U)
		nominateur = innerProduct(TxU,vect)
		denominateur = innerProduct(TxU,S)
		self.s = nominateur / denominateur
		
	def compute_t(self,p0,S,T,U):
		""" compute the t coordinate of the point
		p0: origin point
		S,T,U: coordinates of the base vectors """
		vect = Vector(self.x-p0.x,self.y-p0.y,self.z-p0.z)
		SxU = vectorialProduct(S,U)
		nominateur = innerProduct(SxU,vect)
		denominateur = innerProduct(SxU,T)
		self.t = nominateur / denominateur
		
	def compute_u(self,p0,S,T,U):
		""" compute the u coordinate of the point
		p0: origin point
		S,T,U: coordinates of the base vectors """
		vect = Vector(self.x-p0.x,self.y-p0.y,self.z-p0.z)
		SxT = vectorialProduct(S,T)
		nominateur = innerProduct(SxT,vect)
		denominateur = innerProduct(SxT,U)
		self.u = nominateur / denominateur
		
class CtrlPoint(Point):
	""" FFD control point """
	
	def __init__(self,x,y,z,i,j,k,l,m,n):
		self.x = x
		self.y = y
		self.z = z
		self.s = i/l # x coord relative to control points
		self.t = j/m # y coord relative to control points
		self.u = k/n # z coord relative to control points
		self.i = i
		self.j = j
		self.k = k
		self.l = l
		self.m = m
		self.n = n
	
class Vector:
	""" a 3D vector """
	
	def __init__(self,X,Y,Z):
		self.X = X # component in X
		self.Y = Y # component in Y
		self.Z = Z # component in Z

# ------------------------------------------------------
# methods declaration

def comb(a,b):
	""" combination of 2 scalars, b among a """
	if b == 0:
		return 1
	if a == 0:
		return 0
	return (comb(a-1,b-1) + comb(a-1,b))

def vectorialProduct(v1,v2):
	""" compute the vectorial product of 2 vectors """
	result = Vector(0,0,0)
	result.X = v1.Y*v2.Z - v1.Z*v2.Y
	result.Y = v1.Z*v2.X - v1.X*v2.Z
	result.Z = v1.X*v2.Y - v1.Y*v2.X
	return result

def innerProduct(v1,v2):
	""" compute the inner product of 2 vectors """
	scalar = v1.X*v2.X + v1.Y*v2.Y + v1.Z*v2.Z
	return scalar

def compute_s(p1,p0,S,T,U):
    """ DEPRECIATED """
    """ use the object oriented version instead """
    """ compute the s coefficient of a point """
    vect = Vector(p1.x-p0.x,p1.y-p0.y,p1.z-p0.z)
    TxU = vectorialProduct(T,U)
    nominateur = innerProduct(TxU,vect)
    denominateur = innerProduct(TxU,S)
    s = nominateur / denominateur
    return s

def compute_t(p1,p0,S,T,U):
    """ DEPRECIATED """
    """ use the object oriented version instead """
    """ compute the u coefficient of a point """
    vect = Vector(p1.x-p0.x,p1.y-p0.y,p1.z-p0.z)
    SxU = vectorialProduct(S,U)
    nominateur = innerProduct(SxU,vect)
    denominateur = innerProduct(SxU,T)
    t = nominateur / denominateur
    return t
	
def compute_u(p1,p0,S,T,U):
    """ DEPRECIATED """
    """ use the object oriented version instead """
    """ compute the u coefficient of a point """
    vect = Vector(p1.x-p0.x,p1.y-p0.y,p1.z-p0.z)
    SxT = vectorialProduct(S,T)
    nominateur = innerProduct(SxT,vect)
    denominateur = innerProduct(SxT,U)
    u = nominateur / denominateur
    return u

def fromList2BoxCtrlPts(listCtrlPoints):
    """ transform a list of control points
    into a box of control points """
    boxCtrlPts = list()
    jPrev = 0 # j index of the previous control point
    ConstIList = list() # list of control points of constant i index
    for CtrlPt in listCtrlPoints:
        if CtrlPt.j == jPrev:
            ConstIList.append(CtrlPt)
        else:
            boxCtrlPts.append(ConstIList)
            ConstIList = [CtrlPt]
        jPrev = CtrlPt.j
    boxCtrlPts.append(ConstIList)
    return boxCtrlPts
    
def fromBox2ListCtrlPts(boxCtrlPts):
    """ transform a box of control points
    into a list of control points """
    listCtrlPts = list()
    for ConstIList in boxCtrlPts:
        for CtrlPt in ConstIList:
            listCtrlPts.append(CtrlPt)
    return listCtrlPts
    
def computeCoords(boxCtrlPts,listPoints):
	""" compute the relative coords of points in FFD representation """
	for points in listPoints:
		l = boxCtrlPts[0][0].l
		m = boxCtrlPts[0][0].m
		n = boxCtrlPts[0][0].n
		i = 0
		resultiX = 0
		resultiY = 0
		while i <= l:
			j = 0
			resultjX = 0
			resultjY = 0
			while j <= m:
				resultjX = resultjX + (comb(m,j) * ((1-points.t)**(m-j)) * ((points.t)**j) * boxCtrlPts[j][i].x)
				resultjY = resultjY + (comb(m,j) * ((1-points.t)**(m-j)) * ((points.t)**j) * boxCtrlPts[j][i].y)
				j = j + 1
			resultiX = resultiX + (comb(l,i) * ((1-points.s)**(l-i)) * ((points.s)**i)) * resultjX
			resultiY = resultiY + (comb(l,i) * ((1-points.s)**(l-i)) * ((points.s)**i)) * resultjY
			i = i + 1
		points.x = resultiX
		points.y = resultiY

# ------------------------------------------------------
# IO functions

def read(filePath):
	""" read and format the specified file to
	be modify by FFD,
	inputs:
	filePath: string of the relative path
	of the file to read
	outputs:
	listPts: list of Point defining the airfoil """
	
	file = open(filePath, 'r')
	content = file.readlines()
	file.close()
	
	listPts = list()
	i = 1
	for lines in content: # formating
		if i > 1:
			lines = lines.lstrip()
			lines = lines.strip('\n')
			lines = lines.split(' ')
			while '' in lines:
				lines.remove('')
			floatLines = list()
			for strCoord in lines:
				floatLines.append(float(strCoord))
			point = Point(floatLines[0],floatLines[1],0.0)
			listPts.append(point)
		i += 1
	
	return listPts

def writePtsXFOIL(filePath,listPts):
    """ write the x,y coordinates of a list of points
    compatible with XFOIL formating """
    
    coord = filePath.split('.')[0]
    for point in listPts:
        coord += '\n' + '    ' + str(point.x) + '      ' + str(point.y)
    file = open(filePath,'w')
    file.write(coord)
    file.close()

def writeRelativeCoordPts(filePath,listPts):
    """ update the s,t,u coordinates of a list of points
    from the filePath.pts file """
    
    fileRelCoord = filePath.split('.')[0] + '.rel'
    coord = str()
    for point in listPts:
        coord = coord + '\n' + '\t' + str(point.s) + '\t' + str(point.t) + '\t' + str(point.u)
    file = open(fileRelCoord,'w')
    file.write(coord)
    file.close()
    
def updateRelCoordFromFile(filePath,listPts):
    """ write the s,t,u coordinates of a list of points
    in a filePath.pts file """
    
    fileRelCoord = filePath.split('.')[0] + '.rel'
    try:
        file = open(fileRelCoord)
        file.readline()
        newListPts = list()
        for point in listPts:
            relCoord = file.readline().strip('\n').strip().split('\t')
            point.s = float(relCoord[0])
            point.t = float(relCoord[1])
            point.u = float(relCoord[2])
            newListPts.append(point)
        return newListPts
    except IOError:
        return listPts

def readCtrlPoints(filePath):
    """ read the control points defined for the geometry
    input: filePath of the geometry,
    will search for the file 'filePath.ctrl'
    if not found will take the default poistions
    output:
    listCtrlPoints: list of the control points in order
    S then T then U """
    
    fileCtrl = filePath.split('.')[0] + '.ctrl'
    try:
        file = open(fileCtrl, 'r')
    except IOError:
        CrtlPoint1 = CtrlPoint(0.00,-0.05,0,0,0,0,3,1,1)
        CrtlPoint2 = CtrlPoint(0.33,-0.05,0,1,0,0,3,1,1)
        CrtlPoint3 = CtrlPoint(0.66,-0.05,0,2,0,0,3,1,1)
        CrtlPoint4 = CtrlPoint(1.00,-0.05,0,3,0,0,3,1,1)
        CrtlPoint5 = CtrlPoint(0.00,0.05,0,0,1,0,3,1,1)
        CrtlPoint6 = CtrlPoint(0.33,0.05,0,1,1,0,3,1,1)
        CrtlPoint7 = CtrlPoint(0.66,0.05,0,2,1,0,3,1,1)
        CrtlPoint8 = CtrlPoint(1.00,0.05,0,3,1,0,3,1,1)
        listCtrlPts = [CrtlPoint1, CrtlPoint2, CrtlPoint3, CrtlPoint4, 
        CrtlPoint5, CrtlPoint6, CrtlPoint7, CrtlPoint8]
        return listCtrlPts
    
    listCtrlPts = list()
    headline = file.readline()
    l = int(headline.split('\t')[0])
    m = int(headline.split('\t')[1])
    n = int(headline.split('\t')[2])
    for k in range(n):
        for j in range(m+1):
            for i in range(l+1):
                line = file.readline()
                x = float(line.split('\t')[0])
                y = float(line.split('\t')[1])
                z = float(line.split('\t')[2])
                CtrlPt = CtrlPoint(x,y,z,i,j,k,l,m,n)
                listCtrlPts.append(CtrlPt)
    return listCtrlPts

def writeCtrlPoints(filePath,listCtrlPts):
    """ write the control points defined for the geometry
    inputs: 
    filePath of the geometry
    listCtrlPts to write
    """
    l = listCtrlPts[0].l
    m = listCtrlPts[0].m
    n = listCtrlPts[0].n
    
    fileCtrl = filePath.split('.')[0] + '.ctrl'
    file = open(fileCtrl, 'w')
    file.write(str(l) + '\t' + str(m) + '\t' + str(n))
    for CtrlPt in listCtrlPts:
        file.write('\n' + str(CtrlPt.x) + '\t' + str(CtrlPt.y) + '\t' + str(CtrlPt.z))
    file.close()
    
# ------------------------------------------------------
# main FFD function

def FFD(inputGeometry,outputGeometry,x0=0,y0=0,x1=0,y1=0,x2=0,y2=0,x3=0,y3=0,x4=0,y4=0,x5=0,y5=0,x6=0,y6=0,x7=0,y7=0):
    
    listPts = read(inputGeometry)
    listCtrlPts = readCtrlPoints(inputGeometry)
    boxCtrlPts = fromList2BoxCtrlPts(listCtrlPts)
    
    S = Vector(1,0,0)
    T = Vector(0,0.1,0)
    U = Vector(0,0,1)
    origin = Point(listCtrlPts[0].x,listCtrlPts[0].y,listCtrlPts[0].z)
    
    for point in listPts:
        point.compute_s(origin,S,T,U)
        point.compute_t(origin,S,T,U)
        point.compute_u(origin,S,T,U)
        #point.s = compute_s(point,origin,S,T,U)
        #point.t = compute_t(point,origin,S,T,U)
        #point.u = compute_u(point,origin,S,T,U)
    listPts = updateRelCoordFromFile(inputGeometry,listPts)
    
    listCtrlPts[0].x += x0
    listCtrlPts[0].y += y0
    listCtrlPts[1].x += x1
    listCtrlPts[1].y += y1
    listCtrlPts[2].x += x2
    listCtrlPts[2].y += y2
    listCtrlPts[3].x += x3
    listCtrlPts[3].y += y3
    listCtrlPts[4].x += x4
    listCtrlPts[4].y += y4
    listCtrlPts[5].x += x5
    listCtrlPts[5].y += y5
    listCtrlPts[6].x += x6
    listCtrlPts[6].y += y6
    listCtrlPts[7].x += x7
    listCtrlPts[7].y += y7
    
    boxCtrlPts = fromList2BoxCtrlPts(listCtrlPts)

    computeCoords(boxCtrlPts,listPts)
    
    writePtsXFOIL(outputGeometry,listPts)
    writeRelativeCoordPts(outputGeometry,listPts)
    writeCtrlPoints(outputGeometry,listCtrlPts)
    
def FFD_LEthickness(inputGeometry,outputGeometry,thickness=0):
    x0 -= thickness / 2
    x4 += thickness / 2
    FFD(inputGeometry,outputGeometry,x0,0,0,0,0,0,0,0,x4,0,0,0,0,0,0,0)

def FFD_ctrlPts(inputGeometry,outputGeometry,outputListCtrlPts):
    
    listPts = read(inputGeometry)
    listCtrlPts = readCtrlPoints(inputGeometry)
    boxCtrlPts = fromList2BoxCtrlPts(listCtrlPts)
    
    S = Vector(1,0,0)
    T = Vector(0,0.1,0)
    U = Vector(0,0,1)
    origin = Point(listCtrlPts[0].x,listCtrlPts[0].y,listCtrlPts[0].z)
    
    for point in listPts:
        point.compute_s(origin,S,T,U)
        point.compute_t(origin,S,T,U)
        point.compute_u(origin,S,T,U)
	    #point.s = compute_s(point,origin,S,T,U)
	    #point.t = compute_t(point,origin,S,T,U)
        #point.u = compute_u(point,origin,S,T,U)
    listPts = updateRelCoordFromFile(inputGeometry,listPts)
    
    boxCtrlPts = fromList2BoxCtrlPts(outputListCtrlPts)

    computeCoords(boxCtrlPts,listPts)
    
    writePtsXFOIL(outputGeometry,listPts)
    writeRelativeCoordPts(outputGeometry,listPts)
    writeCtrlPoints(outputGeometry,outputListCtrlPts)
    
# ------------------------------------------------------
# acutal start of the programm


if __name__ == '__main__':
    import matplotlib.pyplot as plt # for ploting data
    fileName = raw_input('Enter an airfoil file name: ')
    if fileName == '':
        fileName = 'naca0009.dat'
    
    listPts = read(fileName)
    listCtrlPts = readCtrlPoints(fileName)
    boxCtrlPts = fromList2BoxCtrlPts(listCtrlPts)

    S = Vector(1,0,0)
    T = Vector(0,0.1,0)
    U = Vector(0,0,1)
    origin = Point(listCtrlPts[0].x,listCtrlPts[0].y,listCtrlPts[0].z)

    for point in listPts:
	    point.s = compute_s(point,origin,S,T,U)
	    point.t = compute_t(point,origin,S,T,U)
	    point.u = compute_u(point,origin,S,T,U)
    listPts = updateRelCoordFromFile(fileName,listPts)

    plot = raw_input('plot ? ')
    if plot == 'y':
	    fig = plt.figure() # create figure
	    ax = fig.add_subplot(111) # create a subplot
	    ax.set_xlabel('X')
	    ax.set_ylabel('Y')
	    ax.set_xlim(0, 1)
	    ax.set_ylim(-0.5, 0.5)
	    x = []
	    for points in listPts:
	        x.append(points.x)
	    y = []
	    for points in listPts:
		    y.append(points.y)
	    xCtrl = []
	    for points in listCtrlPts:
		    xCtrl.append(points.x)
	    yCtrl = []
	    for points in listCtrlPts:
		    yCtrl.append(points.y)
	    ax.plot(x, y,xCtrl,yCtrl,'rs')
	    plt.show()

# ------------------------------------------------------	
# Free Form Deformation
    redo = 'y'
    while redo == 'y':
        change = raw_input('change control points ? ')
        if change == 'y':
	        deltaX1 = float(raw_input('x1 = '))
	        deltaY1 = float(raw_input('y1 = '))
	        deltaX2 = float(raw_input('x2 = '))
	        deltaY2 = float(raw_input('y2 = '))
	        deltaX3 = float(raw_input('x3 = '))
	        deltaY3 = float(raw_input('y3 = '))
	        deltaX4 = float(raw_input('x4 = '))
	        deltaY4 = float(raw_input('y4 = '))
        else:
	        deltaX1 = 0.0
	        deltaY1 = 0.0
	        deltaX2 = 0.0
	        deltaY2 = 0.0
	        deltaX3 = 0.0
	        deltaY3 = 0.0
	        deltaX4 = 0.0
	        deltaY4 = 0.0
	
        listCtrlPts[1].x += deltaX1
        listCtrlPts[1].y += deltaY1
        listCtrlPts[2].x += deltaX2
        listCtrlPts[2].y += deltaY2
        listCtrlPts[5].x += deltaX3
        listCtrlPts[5].y += deltaY3
        listCtrlPts[6].x += deltaX4
        listCtrlPts[6].y += deltaY4
    
        boxCtrlPts = fromList2BoxCtrlPts(listCtrlPts)
    
        computeCoords(boxCtrlPts,listPts)

# ------------------------------------------------------	
# plotting the FFD figure

        plot = raw_input('plot ? ')
        if plot == 'y':	
            fig = plt.figure() # create figure
            ax = fig.add_subplot(111) # create a subplot
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_xlim(0, 1)
            ax.set_ylim(-0.5, 0.5)
            x = []
            for points in listPts:
                x.append(points.x)
            y = []
            for points in listPts:
                y.append(points.y)
            xCtrl = []
            for points in listCtrlPts:
                xCtrl.append(points.x)
            yCtrl = []
            for points in listCtrlPts:
                yCtrl.append(points.y)
            ax.plot(x, y,xCtrl,yCtrl,'rs')
            plt.show()


# ------------------------------------------------------	
# Exporting the new geometry

        export = raw_input('export ? ')
        if export == 'y': 
            nameFile = raw_input('name of file: ')
            if nameFile == '':
                nameFile = 'FFDexport.dat'
            writePtsXFOIL(nameFile,listPts)
            writeRelativeCoordPts(nameFile,listPts)
            writeCtrlPoints(nameFile,listCtrlPts)
        
        redo = raw_input('redo ? ')
        