

import matplotlib.pyplot as plt



def readTextFile(fileName):
    f = open(fileName, 'r')
    flag = 1


    for line in f.readlines():
        temp = line.split()
        if flag==1:
            numCols = len(temp)
            cols = [[] for i in range(numCols)]
            flag =0

        for i in range(numCols):

            cols[i].append(float(temp[i]))
    return cols




def main():
    fileName = "output.txt"
    cols = readTextFile(fileName)



    X = cols[0]
    Y = cols[1]
    Z = cols[2]

    uX = list(set(X))
    uY = list(set(Y))
    uZ = list(set(Z))

    #plt.contour(X, Y, Z)
    #plt.show()

  
    #n = 3 # vegetti regularization:
    n = 1 # zeroth:
    n = 2 # gradient;
    n = 3 # curvature;

    n = 3

    index = cols[n].index(max(cols[n]))
    for i in range(len(cols)):
        print cols[i][index]



if __name__=='__main__':
    main()
