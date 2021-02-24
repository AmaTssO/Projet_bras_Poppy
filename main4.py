import numpy as np
import model4 as fk
import matplotlib.pyplot as plt
import pypot
import numpy as np
import pypot.vrep.io as Vrep_io
import math

class bras_poppy :
    j1,j2,j3,j4,j5,x5,y5,z5,Tm=0,0,0,0,0,0,0,0,0
    Jac,Jac_Inv=0,0
    xt,yt,zt=0,0,0
    Poppy=0
    re_set=0
    
    
    def initState(self):
       # pypot.vrep.close_all_connections()
       #self.Poppy = Vrep_io.VrepIO ('127.0.0.1',19997,'/home/oussama/projet/test/test3/poppy_humanoid.ttt')
        #self.Poppy.start_simulation()
        self.re_set=1
        self.j1=90
        self.j2=90
        self.j3=0
        self.j4=0
        self.j5=0
        self.drawfk(self.j1,self.j2,self.j3,self.j4,self.j5)   
        

    def start_ik(self,x,y,z):
        self.xt=x
        self.yt=y
        self.zt=z
        self.re_set=0
        self.pinj_ik()

    def pinj_ik(self):
        if self.re_set==0:      
            x_start=self.x5
            y_start=self.y5
            z_start=self.z5
            
            ptarget=np.vstack([self.xt, self.yt, self.zt])
            pstart=np.vstack([x_start,y_start,z_start])
            
            delta=fk.dist(ptarget, pstart)
            step=0.01
            DisXYZ=delta[4]
            print("=============")
            #print(DisXYZ)
            if abs(DisXYZ)>0.05:
                pembagi_step=abs(DisXYZ)/step
                #print(pembagi_step)
                dXYZ=(delta[0:3]/pembagi_step)
                #print(dXYZ)
                print("=============")
            else:
                dXYZ=delta[0:3]
                
            if abs(DisXYZ)<=0.001:
              
                print("position atteint :\n",ptarget)
                print("X :",self.x5)
                print("Y :",self.y5)
                print("Z :",self.z5)
                
           #     self.Poppy.set_motor_position("l_shoulder_y",np.deg2rad(self.j1))
            #    self.Poppy.set_motor_position("l_shoulder_x",np.deg2rad(self.j2))
            #    self.Poppy.set_motor_position("l_arm_z",np.deg2rad(self.j3))
            #    self.Poppy.set_motor_position("l_elbow_y",np.deg2rad(self.j4))
                
            else:
                self.Jac=fk.jacobian(self.Tm)
                self.Jac_Inv=fk.PinvJac(self.Jac)
            #    print(self.Jac_Inv)
                dTheta=np.dot(self.Jac_Inv,dXYZ)
               
                dTheta1=np.rad2deg(dTheta[0])
                dTheta2=np.rad2deg(dTheta[1])
                dTheta3=np.rad2deg(dTheta[2])
                dTheta4=np.rad2deg(dTheta[3])
                dTheta5=np.rad2deg(dTheta[4])
                
                self.j1=self.j1
                self.j2=self.j2+dTheta2
                self.j3=self.j3+dTheta3
                self.j4=self.j4+dTheta4
                self.j5=self.j5+dTheta5
             #   self.Poppy.set_motor_position("l_shoulder_y",np.deg2rad(self.j1))
             #  self.Poppy.set_motor_position("l_shoulder_x",np.deg2rad(self.j2))
             #  self.Poppy.set_motor_position("l_arm_z",np.deg2rad(self.j3))
             #  self.Poppy.set_motor_position("l_elbow_y",np.deg2rad(self.j4))
                
                print(self.j1,self.j2,self.j3,self.j4,self.j5)
             
                self.drawfk(self.j1,self.j2,self.j3,self.j4,self.j5)
                self.pinj_ik()

    def drawfk(self,a,b,c,d,e): 
        dh=fk.dh_valeur(a,b,c,d,e)
       # print(dh)
        self.Tm=fk.matrice_transf(dh)
       # print(self.Tm)
        ee=fk.mt_xyz(self.Tm)
       # print(ee)
        X1,Y1,Z1=ee[0,0:2],ee[1,0:2],ee[2,0:2]
        X2,Y2,Z2=ee[0,1:3],ee[1,1:3],ee[2,1:3]
        X3,Y3,Z3=ee[0,2:4],ee[1,2:4],ee[2,2:4]
        X4,Y4,Z4=ee[0,3:5],ee[1,3:5],ee[2,3:5]
        X5,Y5,Z5=ee[0,4:6],ee[1,4:6],ee[2,4:6]
        
        self.x5=X5[1]
        self.y5=Y5[1]
        self.z5=Z5[1]

        fig = plt.figure(num=None, figsize=(14, 20), dpi=80, facecolor='w', edgecolor='k')
        ax = fig.add_subplot(111, projection='3d')
    
        ax.plot(X1,Y1,Z1, color='green', marker='o', linestyle='solid', linewidth=5, markersize=20)
        ax.plot(X2,Y2,Z2, color='green', marker='o', linestyle='solid', linewidth=5, markersize=20)
        ax.plot(X3,Y3,Z3, color='blue', marker='o', linestyle='solid', linewidth=5, markersize=20)
        ax.plot(X4,Y4,Z4, color='black', marker='o', linestyle='solid', linewidth=5, markersize=20)
        ax.plot(X5,Y5,Z5, color='red', marker='o', linestyle='solid', linewidth=5, markersize=20)
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_xlim([-0.50,0.50])
        ax.set_ylim([-0.50,0.50])
        ax.set_zlim([-0.50,0.50])
        ax.set_title('bras')
        ax.view_init(20, 40)
        plt.show()
   
       

X=bras_poppy()
X.initState()
X.start_ik(0.0367,0.336,0.21)
