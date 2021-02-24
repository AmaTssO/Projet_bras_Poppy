import numpy as np
import matplotlib.pyplot as plt
import radtodeg as tr


def dh_valeur(j1,j2,j3,j4,j5):
    j=np.array([[j1, j2, j3, j4, j5],
                [0, 0, 0, -0.01, 0],
                [0, 0.185, -0.148 , 0, 0.215],
                [0, -90, 90, -90, -90]])
    return j    

def dh(theta,d,a,alpha):
    dhm=np.array([[tr.cosd(theta),-tr.sind(theta)*tr.cosd(alpha),tr.sind(theta)*tr.sind(alpha),a*tr.cosd(theta)],
    [tr.sind(theta),tr.cosd(theta)*tr.cosd(alpha),-tr.cosd(theta)*tr.sind(alpha),a*tr.sind(theta)],
    [0,tr.sind(alpha),tr.cosd(alpha),d],
    [0,0,0,1]])
    return dhm

def matrice_transf(j): 
    T01=dh(j[0,0],j[1,0],j[2,0],j[3,0])
    T12=dh(j[0,1],j[1,1],j[2,1],j[3,1])
    T23=dh(j[0,2],j[1,2],j[2,2],j[3,2])
    T34=dh(j[0,3],j[1,3],j[2,3],j[3,3])
    T45=dh(j[0,4],j[1,4],j[2,4],j[3,4])
    
    T02=np.dot(T01,T12)
    T03=np.dot(T02,T23)
    T04=np.dot(T03,T34)
    T05=np.dot(T04,T45)

    return np.hstack((T01, T02, T03, T04, T05))

def mt_xyz(m): 
    Q1=np.array([0,m[0,3],m[0,7],m[0,11],m[0,15],m[0,19]])
    Q2=np.array([0,m[1,3],m[1,7],m[1,11],m[1,15],m[1,19]])
    Q3=np.array([0,m[2,3],m[2,7],m[2,11],m[2,15],m[2,19]])
    Q=np.vstack((Q1,Q2,Q3))
    return Q

def mgd(Q):
    X,Y,Z=Q[0,:],Q[1,:],Q[2,:]
    mgd=np.array([X[5],Y[5],Z[5]])
    return mgd


def jacobian(dhkine):
    Z0=np.array([[0],[0],[1]])
    O=np.array([[0],[0],[0]])
    O5=dhkine[0:3,19:]
    Jac1=np.cross(Z0,(O5-O), axis=0)

    Z1=dhkine[0:3,2:3]
    O1=dhkine[0:3,3:4]
    Jac2=np.cross(Z1,(O5-O1), axis=0)

    Z2=dhkine[0:3,6:7]
    O2=dhkine[0:3,7:8]
    Jac3=np.cross(Z2,(O5-O2), axis=0)

    Z3=dhkine[0:3,10:11]
    O3=dhkine[0:3,11:12]
    Jac4=np.cross(Z3,(O5-O3), axis=0)
    
    Z4=dhkine[0:3,14:15]
    O4=dhkine[0:3,15:16]
    Jac5=np.cross(Z4,(O5-O4), axis=0)
    return np.hstack((Jac1,Jac2,Jac3,Jac4,Jac5))

def PinvJac(jacobian):
    return np.linalg.pinv(jacobian)

def dist(ptarget, pstart):
    dXYZ=ptarget-pstart
    dX=dXYZ[0]
    dY=dXYZ[1]
    dZ=dXYZ[2]
    DisXY=np.sqrt(np.power(dX,2)+np.power(dY,2))#phytagoras x dan y target dan start sehingga didapat jarak x-y
    DisXYZ=np.sqrt(np.power(DisXY,2)+np.power(dZ,2))
    AngXY=np.rad2deg(np.arctan2(dY,dX)) #nilai sudut (drajat)
    AngZ=np.rad2deg(np.arctan2(dZ,DisXY))
    dx=tr.cosd(AngXY)*DisXY
    dy=tr.sind(AngXY)*DisXY
    dz=tr.sind(AngZ)*DisXYZ
    return np.vstack([dx,dy,dz,DisXY,DisXYZ,AngZ])

poppy_head_config = {
    'controllers': {
        'my_dxl_controller': {
            'port': '/dev/ttyACM0',  # Depends on your OS
            'sync_read': False,
            'attached_motors': ['head'],  # You can mix motorgroups or individual motors
        },
    },

    'motorgroups': {
        'head': ['head_y'],
    },

    'motors': {
        'head_y': {
            'id': 1,
            'type': 'AX-12',
            'orientation': 'INdirect',
            'offset': 20.0,
            'angle_limit': (-45, 6),
        },
    },
}

