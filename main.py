import RF
import numpy as np
import math
import json
from tqdm import tqdm
import time


###
#	1. Running "main.py" for calculating a mode with one surface 
#	2. 



###
def main():
	#########################################
	##		Input file here
	####
	#			data input
	fSurface = open("kqQL91.txt", 'r')
	count = 0
	Surface = dict()
	pointList = list()
	while True:
		temp = [float(i) for i in fSurface.readline().split()]
		if len(temp) == 0:
			break
		if len(temp) == 1: 
			Surface[count] = pointList
			pointList.clear()
			count += 1
		else:
			pointList.append(temp)
	fSurface.close()
	print(Surface)
	#		props input
	fProps = open("props.txt","r")
	props = json.loads(fProps.read())
	fProps.close()

		###			props settings
	#############
	#		soil_props is permanently exist in the whole runtime !!!
	soil = RF.soil_props(float(props["soil_props"]["c"]),
						float(props["soil_props"]["phi"]),
						float(props["soil_props"]["gamma"]))
	##
	#		this instance use constant calculating props
	cal = RF.calculating_props(
		float(props["calculating_props"]["u"]),
		float(props["calculating_props"]["kW"]),
		float(props["calculating_props"]["A"]),
		float(props["calculating_props"]["D"]),
		float(props["calculating_props"]["omega"])
	)
	##
	

	######################################################################
	#--------------------------------------------
	FS_Surface = np.array([float(0)] * count)
	lamdaSurface = np.array([float(0)] * count)
	centerSurface = np.array([float(0)] * 2 * count).reshape(count, 2)
	RSurface = np.array([float(0)] * count)


	for i in tqdm(range(count),
				desc="Loading…",  
               	ascii=False, ncols=75):
		time.sleep(0.01)								#	Every Surface
		currData = Surface[i]						# this in list type
		innerData = RF.depth_converter(currData)	# this change to np arr

		###########################	DO SOMETHING TO GIVE FS AND BW	##############
		#		Model setting
		####

		numberCenter = innerData.shape[0] // 10
		isLeft, first = RF.index_first(innerData)
		dx = innerData[1, 0] - innerData[0, 0]
		centerArr = RF.center_defining(innerData,
									numberCenter,
									first,
									isLeft,
									dx)
		RL = RF.radius_lines_defining(innerData,
									numberCenter,
									np.amin(innerData[:, 0]),
									dx)
		R = np.array([RF.Radius(centerArr[numberCenter * i + 1, 1],
										RL[i]) for i in range(numberCenter)])
		FS = np.array([float(0)] * 2 * R.shape[0] * centerArr.shape[0]).reshape(centerArr.shape[0], R.shape[0], 2)
		
		###				
		#######
		##			lamda f(x) Tolerance setting:
		##
		fx = np.array([(k / innerData.shape[0]) * math.pi for k in range(1, innerData.shape[0])])
		setting = RF.setting(
			fx,
			np.array([(i + 1) * (1 / numberCenter) for i in range(numberCenter)]),
			float(props["setting"]["Tolerance"])
		)

		############


		for j in range(centerArr.shape[0]):

			###
			#		slide_props depends on the current surface state so must be declared in every time loop to innerData
			alpha = np.array([0] * (innerData.shape[0] - 1))
			beta = np.array([0] * (innerData.shape[0] - 1))
			a = np.array([0] * (innerData.shape[0] - 1))
			x = np.array([0] * (innerData.shape[0] - 1))
			W = np.array([0] * (innerData.shape[0] - 1))
			for jj in range (innerData.shape[0] - 1):
				
				alpha[jj] = math.atan((innerData[jj + 1, 1] - innerData[jj, 1]) / dx)
				beta[jj] = dx / (math.cos(alpha[jj]))
				a[jj] = 0
				x[jj] = abs(centerArr[j, 0] - innerData[jj, 0] - dx / 2)
				W[jj] = soil.gamma * (innerData[jj + 1, 1] - innerData[jj, 1]) * dx  #?#?/#?#


			slide = RF.slide_props(alpha, beta, a, x, W)
			####
			for k in range(R.shape[0]):
				FS[j, k] = RF.Calculating_FoS(innerData, setting.lamda, setting.fx, centerArr[j], R[k], soil, cal, slide, setting.Tolerance)
		
		FS_Surface[i] = np.amin(FS)
		index = np.where(FS == FS_Surface[i])
		RSurface[i] = R[index[0, 1]] 	
		centerSurface[i] = centerArr[index[0, 0]]
		
		#########################	END DO SOMETHING				##############
	BW = np.array([float(0)] * count)
	for i in range(count):
		if FS_Surface[i] < float(props["setting"]["FSCritical"]):
			currData = Surface[i]						# this in list type
			innerData = RF.depth_converter(currData)
			dx = innerData[1, 0] - innerData[0, 0]
			BW[i] = RF.Cal_BW(innerData, centerSurface[i], RSurface[i], dx)
	print(FS_Surface)
	print(BW)


	

if __name__ == "__main__":
	main()
