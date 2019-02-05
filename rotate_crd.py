import numpy as np
import math

all_crd = []
center_crd = []
rotated_crd = []
prefix = []
appendix = []
matrixA = np.array([[],[],[]])
matrixB = np.array([[],[],[]])

def vector_gen(atomX,atomY):

    vector = [0]*3
    vector[0] = float(all_crd[atomY*3-3]) - float(all_crd[atomX*3-3])
    vector[1] = float(all_crd[atomY*3-2]) - float(all_crd[atomX*3-2])
    vector[2] = float(all_crd[atomY*3-1]) - float(all_crd[atomX*3-1])

    return vector


def import_crd(file_name,center_atm):

    _inp_crd = []
    inp_crd = []
    with open (file_name, 'r') as fref:
        for i, line in enumerate(fref.readlines()):
            if i == 0:
                prefix.extend(line)
            elif i == 1:
                prefix.extend(line)
                tmp = line.split()
                total_atom_num = int(tmp[0])
            elif (i < 2 + total_atom_num/2):
                _inp_crd.extend(line.split())
            else:
                appendix.extend(line)

    del _inp_crd[total_atom_num*3:]

    for i in range(len(_inp_crd)):
        tmp2 = float(_inp_crd[i]) - float(_inp_crd[int(center_atm)*3-3+i%3])
        inp_crd.append(tmp2)

    return inp_crd, _inp_crd[int(center_atm)*3-3:int(center_atm)*3] 


def determine_rotation(atomA,atomB,atomC):

    vec1 = vector_gen(atomA,atomB)
    vec2 = vector_gen(atomA,atomC)
    _vert = np.cross(vec1,vec2)
    vert = [0]*3
    for i in range(len(vert)):
        vert[i] = _vert[i]/math.sqrt(_vert[0]*_vert[0]+_vert[1]*_vert[1]+_vert[2]*_vert[2])
    beta = math.asin(vert[2])
    alpha = math.asin(vert[1]/math.cos(beta))
    _axis1 = np.cross([0,0,1],vert)
    axis1 = [0]*3
    for i in range(len(axis1)):
        axis1[i] = _axis1[i]/math.sqrt(_axis1[0]*_axis1[0]+_axis1[1]*_axis1[1]+_axis1[2]*_axis1[2])
    sina = math.sin(-alpha)
    cosa = math.cos(-alpha)
    matrixA = np.array([[cosa,-sina,0],[sina,cosa,0],[0,0,1]])
    inv_matrixA = np.linalg.inv(matrixA)
    sinb = math.sin(beta)
    cosb = math.cos(beta)
    x,y,z = axis1[0],axis1[1],axis1[2]
    matrixB = np.array([[(cosb+x*x*(1-cosb)),(x*y*(1-cosb)-z*sinb),(z*x*(1-cosb)+y*sinb)],[(x*y*(1-cosb)+z*sinb),(cosb+y*y*(1-cosb)),(y*z*(1-cosb)-x*sinb)],[(z*x*(1-cosb)-y*sinb),(y*z*(1-cosb)+x*sinb),(cosb+z*z*(1-cosb))]])  
    inv_matrixB = np.linalg.inv(matrixB)
    _temp = np.dot(matrixB,vert)
    temp = np.dot(matrixA,_temp)
    print (temp)
    return matrixA, matrixB


def rotate_crd():

    crd_mod = [0]*3
    center = np.array(center_crd)
    print center
    for i in range(len(all_crd)/3):
        crd = np.array([float(all_crd[i*3]),float(all_crd[i*3+1]),float(all_crd[i*3+2])])
        _temp = np.dot(matrixB,crd)
        temp = np.dot(matrixA,_temp)
        for i in range(len(center)):
            crd_mod[i] = float(center[i]) + temp[i]
        #crd_mod = center + temp
        rotated_crd.extend(crd_mod)
    crd = np.array([float(all_crd[53*3-3]),float(all_crd[53*3-2]),float(all_crd[53*3-1])])
    print (crd)
    _temp = np.dot(matrixB,crd)
    temp = np.dot(matrixA,_temp)
    print (temp)
    print (float(center[0]) + temp[0])
        
    print (rotated_crd[53*3-3])
    return rotated_crd
        

def output_crd():
    with open ("rmout_00210cycle_mod.crd","w") as fout:
        for i in range(len(prefix)):
            fout.write(prefix[i])
        for i in range(len(all_crd)/6):
            for j in range(6):
                rotated_crd[i*6+j] = "{0:.7f}".format(rotated_crd[i*6+j])
                rotated_crd[i*6+j] = rotated_crd[i*6+j].rjust(11)
                fout.write(" "+str(rotated_crd[i*6+j]))
            fout.write("\n")
        for i in range(len(appendix)):
            fout.write(appendix[i])
    

if __name__ == '__main__':
    all_crd, center_crd = import_crd("rmout_00210cycle.crd",53)
    print all_crd[53*3-3]
    matrixA, matrixB = determine_rotation(53,23,24)
    rotated_crd = rotate_crd()
    output_crd()

