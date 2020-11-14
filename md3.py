import numpy as np
import random as rnd
import matplotlib.pyplot as plt
import math



def kin(vx,vy,vz,mass,ntot):
	sumv=0.0
	for i in range(ntot):
		sumv=sumv+vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]
	ke=0.5*mass*sumv
	return ke


def pbc(r):
	if r>0:
		sign =1
	else:
		sign =-1
	r=r-sign
	
	return r


def force(x,y,z,ntot,lx,ly,lz,fx,fy,fz,E,ro,m,n):

	cc=m*n/(n-m)
	pe=0.0

	for i in range(0,ntot,1):
		for j in range(i+1,ntot,1):
			dx=x[j]-x[i]
			dy=y[j]-y[i]
			dz=z[j]-z[i]

			dx=dx/lx
			dy=dy/ly
			dz=dz/lz

			if(math.fabs(dx) > 0.5): dx = pbc(dx)
			if(math.fabs(dy) > 0.5): dy = pbc(dy)
			if(math.fabs(dz) > 0.5): dz = pbc(dz)

			dx=dx*lx
			dy=dy*ly
			dz=dz*lz

			r2=dx*dx+dy*dy+dz*dz
			r=np.sqrt(r2)

			sc=ro/r

			fx[i] = fx[i]+((cc*E*sc)*(pow(sc,m-1)-pow(sc,n-1))*(dx/r2))
			fy[i] = fy[i]+((cc*E*sc)*(pow(sc,m-1)-pow(sc,n-1))*(dy/r2))
			fz[i] = fz[i]+((cc*E*sc)*(pow(sc,m-1)-pow(sc,n-1))*(dz/r2))

			fx[j] = fx[j]-((cc*E*sc)*(pow(sc,m-1)-pow(sc,n-1))*(dx/r2))
			fy[j] = fy[j]-((cc*E*sc)*(pow(sc,m-1)-pow(sc,n-1))*(dy/r2))
			fz[j] = fz[j]-((cc*E*sc)*(pow(sc,m-1)-pow(sc,n-1))*(dz/r2))

			msc=pow(sc,m)
			nsc=pow(sc,n)
			pe = pe + ((E/(n-m))*((m*nsc)-(n*msc)))
	return pe



def vverlet(lx,ly,lz,x,y,z,fx,fy,fz,mass,dt,ntot,E,ro,m,n,vx,vy,vz):
	
	for i in range(0,ntot,1):
		x[i]=x[i] + vx[i]*dt + 0.5*fx[i]/mass*dt*dt
		y[i]=y[i] + vy[i]*dt + 0.5*fy[i]/mass*dt*dt
		z[i]=z[i] + vz[i]*dt + 0.5*fz[i]/mass*dt*dt


	fxnew=[]
	fynew=[]
	fznew=[]


	fxnew=[0 for i in range(0,ntot,1)]
	fynew=[0 for i in range(0,ntot,1)]
	fznew=[0 for i in range(0,ntot,1)]

	force(x,y,z,ntot,lx,ly,lz,fxnew,fynew,fznew,E,ro,m,n)

	for i in range(0,ntot,1):
		vx[i]=vx[i]+dt*0.5*((fx[i]/mass)+(fxnew[i]/mass))
		vy[i]=vy[i]+dt*0.5*((fy[i]/mass)+(fynew[i]/mass))
		vz[i]=vz[i]+dt*0.5*((fz[i]/mass)+(fznew[i]/mass))






def linmom(vx,vy,vz,ntot):
	sumvx=0.0
	sumvy=0.0
	sumvz=0.0

	for i in range(0,ntot,1):
		sumvx=sumvx+vx[i]
		sumvy=sumvy+vy[i]
		sumvz=sumvz+vz[i]

	for i in range(0,ntot,1):
		vx[i]=vx[i]-sumvx/ntot;
		vy[i]=vy[i]-sumvy/ntot;
		vz[i]=vz[i]-sumvz/ntot;






def TemScale(vx,vy,vz,ntot,Tdesire,Trand):
	for i in range(0,ntot,1):
		vx[i]=vx[i]*np.sqrt(Tdesire/Trand)
		vy[i]=vy[i]*np.sqrt(Tdesire/Trand)
		vz[i]=vz[i]*np.sqrt(Tdesire/Trand)




def Temp(vx,vy,vz,ntot,mass,kB):
	sumV=0.0
	for i in range(0,ntot,1):
		sumV=sumV+vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]

	temp=(mass*sumV)/(3*kB*ntot)

	return temp



def randomVel(vx,vy,vz,ntot):
	for i in range(0,ntot,1):
		vx.append(rnd.uniform(-0.5,0.5))
		vy.append(rnd.uniform(-0.5,0.5))
		vz.append(rnd.uniform(-0.5,0.5))





def createFcc(x,y,z,nx,ny,nz,ntot,nunit,a,lx,ly,lz):
	u=[[0 for i in range(3)] for j in range(nunit)]

	u[0][0]=0.0*a;u[0][1]=0.0*a;u[0][2]=0.0*a
	u[1][0]=0.5*a;u[1][1]=0.5*a;u[1][2]=0.0*a
	u[2][0]=0.5*a;u[2][1]=0.0*a;u[2][2]=0.5*a
	u[3][0]=0.0*a;u[3][1]=0.5*a;u[3][2]=0.5*a


	for i in range(0,nx,1):
		for j in range(0,ny,1):
			for k in range(0,nz,1):
				for l in range(0,nunit,1):
					x.append(u[l][0]+i*a)
					y.append(u[l][1]+j*a)
					z.append(u[l][2]+k*a)


	



if __name__ == "__main__" :
	a=3.6e-10
	nx=2
	ny=2
	nz=2
	nunit=4
	x=[]
	y=[]
	z=[]
	fx=[]
	fy=[]
	fz=[]
	vx=[]
	vy=[]
	vz=[]
	m=5
	n=8
	dt=1e-15
	ro=2.5487e-10
	epsilon=3401.1
	nstep=5000
	Tdesire=300.0
	mass=105.5206e-27
	kB=1.38064852E-23
	E=epsilon*kB
	ntot=nunit*nx*ny*nz

	
	lx=nx*a
	ly=ny*a
	lz=nz*a

	createFcc(x,y,z,nx,ny,nz,ntot,nunit,a,lx,ly,lz)
        
	randomVel(vx,vy,vz,ntot)

	linmom(vx,vy,vz,ntot)
	


	trand=Temp(vx,vy,vz,ntot,mass,kB)


	TemScale(vx,vy,vz,ntot,Tdesire,trand)

	fp=open("datafile","w+")

	time=0.0
	ksav=0
	while(time<nstep*dt):
		fx=[0.0 for i in range(ntot)]
		fy=[0.0 for i in range(ntot)]
		fz=[0.0 for i in range(ntot)]
		pot=force(x,y,z,ntot,lx,ly,lz,fx,fy,fz,E,ro,m,n)
		vverlet(lx,ly,lz,x,y,z,fx,fy,fz,mass,dt,ntot,E,ro,m,n,vx,vy,vz)
		ke=kin(vx,vy,vz,mass,ntot)
		Temperature=Temp(vx,vy,vz,ntot,mass,kB)
		if(ksav%50==0):
			fp.write("ITEM: TIMESTEP\n")
			fp.write('%e\n'% (time))
			fp.write("ITEM: NUMBER OF ATOMS\n")
			fp.write('%d\n'% (ntot))
			fp.write("ITEM: BOX BOUNDS pp pp pp\n")
			fp.write("0.0 %e\n"% (lx/a))
			fp.write("0.0 %e\n"% (ly/a))
			fp.write("0.0 %e\n"% (lz/a))
			fp.write("ITEM: ATOMS id type x y z\n")
			for i in range(ntot):
				fp.write("%d 1 %e %e %e\n"% (i+1,x[i]/a,y[i]/a,z[i]/a))
				print("%e %e %e %e %lf\n"%(time,pot,ke,pot+ke,Temperature))
		time=time+dt
		ksav=ksav+1
