import math

def computeMI_outer(vectorX, vectorY):
    firstLoop = True
    return computeMI(vectorX, vectorY, 0, len(vectorX), 0, len(vectorY), firstLoop)/len(vectorX)+math.log(len(vectorX))

def computeMI(vectorX, vectorY, limX1, limX2, limY1, limY2, firstLoop):
    mi = 0
    pointsInQuadrant1 = 0
    pointsInQuadrant2 = 0
    pointsInQuadrant3 = 0
    pointsInQuadrant4 = 0

    # Pivots on half quadrants short
    splitX = ((limX2 + limX1) / (2 * 1.0))
    splitY = ((limY2 + limY1) / (2 * 1.0))

    libraryX1 = [0] * len(vectorX)
    libraryX2 = [0] * len(vectorX)
    libraryX3 = [0] * len(vectorX)
    libraryX4 = [0] * len(vectorX)
    libraryY1 = [0] * len(vectorX)
    libraryY2 = [0] * len(vectorX)
    libraryY3 = [0] * len(vectorX)
    libraryY4 = [0] * len(vectorX)

    i = 0
    while i < len(vectorX):
        if vectorX[i] <= splitX:
            if vectorY[i] <= splitY:
                libraryX1[pointsInQuadrant1] = vectorX[i]
                libraryY1[pointsInQuadrant1] = vectorY[i]
                pointsInQuadrant1 += 1
            else:
                libraryX2[pointsInQuadrant2] = vectorX[i]
                libraryY2[pointsInQuadrant2] = vectorY[i]
                pointsInQuadrant2 += 1
        else:
            if vectorY[i] > splitY:
                libraryX3[pointsInQuadrant3] = vectorX[i]
                libraryY3[pointsInQuadrant3] = vectorY[i]
                pointsInQuadrant3 += 1
            else:
                libraryX4[pointsInQuadrant4] = vectorX[i]
                libraryY4[pointsInQuadrant4] = vectorY[i]
                pointsInQuadrant4 += 1
        i += 1

    # Chi-square uniformity test
    expected = len(vectorX)/4.0
    tst = (math.pow(pointsInQuadrant1 - expected, 2) + math.pow(pointsInQuadrant2 - expected, 2) + math.pow(pointsInQuadrant3 - expected, 2) + math.pow(pointsInQuadrant4 - expected, 2))/expected

    # 7.815 is the hard-coded threshold, below which you pass the uniformity Chi-Square test with confidence >95%
    if tst > 7.815 or firstLoop:
        firstLoop = False

	    # Operations to put points in different quadrants
	    # TODO the splitX-X1 can be precomputed for future optimization
        if pointsInQuadrant1 > 2:
            x1 = libraryX1[:pointsInQuadrant1]
            y1 = libraryY1[:pointsInQuadrant1]

            mi += computeMI(x1, y1, limX1, splitX, limY1, splitY, firstLoop)

        else:
            mi += getInformation(pointsInQuadrant1, splitX-limX1, splitY-limY1)

        if pointsInQuadrant2 > 2:
            x2 = libraryX2[:pointsInQuadrant2]
            y2 = libraryY2[:pointsInQuadrant2]

            mi += computeMI(x2, y2, limX1, splitX, splitY, limY2, firstLoop)

        else:
            mi += getInformation(pointsInQuadrant2, splitX - limX1, limY2 - splitY)

        if pointsInQuadrant3 > 2:
            x3 = libraryX3[:pointsInQuadrant3]
            y3 = libraryY3[:pointsInQuadrant3]

            try:
                mi += computeMI(x3, y3, splitX, limX2, splitY, limY2, firstLoop)
            except:
                print('Error')
        else:
            mi += getInformation(pointsInQuadrant3, limX2 - splitX, limY2 - splitY)


        if pointsInQuadrant4 > 2:
            x4 = libraryX4[:pointsInQuadrant4]
            y4 = libraryY4[:pointsInQuadrant4]

            mi += computeMI(x4, y4, splitX, limX2, limY1, splitY, firstLoop)

        else:
            mi += getInformation(pointsInQuadrant4, limX2 - splitX, splitY - limY1)
    else:
        mi += getInformation(len(vectorX), limX2-limX1, limY2-limY1)

    return mi


def getInformation(pXY, pX, pY):
    if pXY == 0:
        return 0
    else:
        return float(pXY)*math.log(float(pXY)/float((pX*pY)))