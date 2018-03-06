

#version 1
#
#
#Setup data structure
#Made timer that includes fps


from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import pyqtgraph as pg
import numpy as np
import random
#import time
from pyqtgraph.ptime import time
import functools

app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.show()
g = gl.GLGridItem()
w.addItem(g)


def nearFunction(mat,i,j,k):

	return mat[i+1,j,k-1] or mat[i,j+1,k-1] or mat[i,j,k-1] or \
		   mat[i-1,j,k] or mat[i,j-1,k] or mat[i,j,k-1] or \
		   mat[i+2,j,k] or mat[i,j+2,k] or \
		   mat[i-2,j,k] or mat[i,j-2,k] or mat[i,j,k-2]   


def makeSeedRand(mat):
	row, col, layer = mat.shape
	for i in range(2, row-2):
		for j in range(2, col-2):
			for k in range(2, layer-2):
				#p = 0.311
				p = 0.211
				randNum = random.uniform(0, 1)

				if(randNum <= p):
					mat[i,j,k] = 1
				#matB[i,j,k] = 1
				# if(1*(row/3) < i and i < 2*(row/3)): #middle third
				# 	if(1*(col/3) < j and j < 2*(col/3)): #middle third
				# 		if(k < 1*(layer/3)): 
				# 		#if(1*(layer/3) < k and k < 2*(layer/3)): #middle third
				# 			randNum = random.randint(0,25)

				# 			if(randNum <= 1):
				# 				mat[i,j,k] = 1
				# 			#matB[i,j,k] = 1

				# else:
				# 	randNum = random.randint(0,250)

				# 	if(randNum <= 1):
				# 		mat[i,j,k] = 1


	


def plantSeed(mat, numSeeds):
	#put in the middle third of box
	row, col, layer = mat.shape
	for i in range(numSeeds):
		rowRand = random.randint(2,row-2);
		colRand = random.randint(2,col-2);
		layerRand = random.randint(2,layer-2);
		mat[rowRand,colRand,layerRand] = 1

def iterateForwardVector():
	humCopy = hum.copy()
	actCopy = act.copy()
	cldCopy = cld.copy()
	row, col, lay = hum.shape

	hum[2:row-2, 2:col-2, 2:lay-2] = humCopy[2:row-2, 2:col-2, 2:lay-2] & (~ actCopy[2:row-2, 2:col-2, 2:lay-2])
	cld[2:row-2, 2:col-2, 2:lay-2] = np.logical_or(cldCopy[2:row-2, 2:col-2, 2:lay-2] , actCopy[2:row-2, 2:col-2, 2:lay-2])

	matR1 = np.roll(np.roll(act,-1,axis=0),1,axis=2) # mat[i+1,j,k-1] 
	matR2 = np.roll(np.roll(act,-1,axis=1),1,axis=2) # mat[i,j+1,k-1]
	matR3 = np.roll(act,1,axis=2) # mat[i,j,k-1]
	matR4 = np.roll(act,1,axis=0) # mat[i-1,j,k]
	matR5 = np.roll(act,1,axis=1) # mat[i,j-1,k]
	matR6 = np.roll(act,1,axis=2) # mat[i,j,k-1]
	matR7 = np.roll(act,-2,axis=0) # mat[i+2,j,k]
	matR8 = np.roll(act,-2,axis=1) # mat[i,j+2,k]
	matR9 = np.roll(act,2,axis=0) # mat[i-2,j,k]
	matR10 = np.roll(act,2,axis=1) # mat[i,j-2,k]
	matR11 = np.roll(act,2,axis=2) # mat[i,j,k-2]

	act[2:row-2, 2:col-2, 2:lay-2] = (~ actCopy[2:row-2, 2:col-2, 2:lay-2]) & humCopy[2:row-2, 2:col-2, 2:lay-2] & \
									np.logical_or(matR1[2:row-2, 2:col-2, 2:lay-2],
									   np.logical_or(matR2[2:row-2, 2:col-2, 2:lay-2],
									   np.logical_or(matR3[2:row-2, 2:col-2, 2:lay-2],
									   np.logical_or(matR4[2:row-2, 2:col-2, 2:lay-2],
									   np.logical_or(matR5[2:row-2, 2:col-2, 2:lay-2],
									   np.logical_or(matR6[2:row-2, 2:col-2, 2:lay-2],
									   np.logical_or(matR7[2:row-2, 2:col-2, 2:lay-2],
									   np.logical_or(matR8[2:row-2, 2:col-2, 2:lay-2],
									   np.logical_or(matR9[2:row-2, 2:col-2, 2:lay-2],
									   np.logical_or(matR10[2:row-2, 2:col-2, 2:lay-2],matR11[2:row-2, 2:col-2, 2:lay-2]))))))))))
								

lenI = 60
lenJ = 60
lenK = 60

hum = np.zeros((lenI, lenJ, lenK))
act = np.zeros((lenI, lenJ, lenK))
cld = np.zeros((lenI, lenJ, lenK))

hum = hum.astype(int)
act = act.astype(int)
cld = cld.astype(int)

makeSeedRand(hum)
plantSeed(act,2)


indexesFinal = np.array([[1,2,3]])
sp2 = gl.GLScatterPlotItem(pos=indexesFinal,size=1.5,pxMode=False)

w.addItem(sp2)

def resetVars():
	global hum, act, cld, indexesFinal
	hum = np.zeros((lenI, lenJ, lenK))
	act = np.zeros((lenI, lenJ, lenK))
	cld = np.zeros((lenI, lenJ, lenK))

	hum = hum.astype(int)
	act = act.astype(int)
	cld = cld.astype(int)

	makeSeedRand(hum)
	plantSeed(act,2)

	indexesFinal = np.array([[1,2,3]])


totalIterations = 80
numIteration    = 0

lastTime = time()
fps = None
def update():
	global numIteration, indexesFinal, lastTime, fps
	if(numIteration < totalIterations) :

		sp2.setData(pos=indexesFinal)
		indexes = np.where(cld==1)
		indexesFinal = np.array([[indexes[0][i],indexes[1][i],indexes[2][i]] for i in range(len(indexes[0]))])
		
		iterateForwardVector()
		numIteration+=1
	else:
		resetVars()
		numIteration = 0

	now = time()
	dt = now - lastTime
	lastTime = now
	if fps is None:
	    fps = 1.0/dt
	else:
	    s = np.clip(dt*3., 0, 1)
	    fps = fps * (1-s) + (1.0/dt) * s
	print('%0.2f fps' % fps)
	


t = QtCore.QTimer()
t.timeout.connect(update)
t.start(5)  

## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, PYQT_VERSION):
        QtGui.QApplication.instance().exec_()
