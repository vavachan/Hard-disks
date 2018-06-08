import numpy as np
import sys
import random
P=np.loadtxt(sys.argv[1])
buff=5.
maxx=max(P[:,0])+buff
maxy=max(P[:,1])+buff
maxz=max(P[:,2])+buff
minx=min(P[:,0])-buff
miny=min(P[:,1])-buff
minz=min(P[:,2])-buff
#print "xmax=",maxx,";"
#print "ymax=",maxy,";"
#print "zmax=",maxz,";"
#print "xmin=",minx,";"
#print "ymin=",miny,";"
#print "zmin=",minz,";"
boxx=maxx-minx
boxy=maxy-miny
boxz=maxz-minz
dp=2
epsilon=1.5
def print_for_file():
    for i in range (0,dp+1):
        for j in range (0,dp+1):
            print minx+i*boxx/dp+random.uniform(-0.5,0.5)*epsilon,"\t",miny+j*boxy/dp+random.uniform(-0.5,0.5)*epsilon,"\t",maxz+random.uniform(-0.5,0.5)*epsilon,"\t",1.0,"\t",1.0
            print minx+i*boxx/dp+random.uniform(-0.5,0.5)*epsilon,"\t",miny+j*boxy/dp+random.uniform(-0.5,0.5)*epsilon,"\t",minz+random.uniform(-0.5,0.5)*epsilon,"\t",1.0,"\t",1.0
            #print "draw sphere {",minx+i*boxx/dp,miny+j*boxy/dp,maxz,"} radius 1.0 resolution 10"
            #print "draw sphere {",minx+i*boxx/dp,miny+j*boxy/dp,minz,"} radius 1.0 resolution 10"
    for i in range (0,dp+1):
        for j in range (1,dp):
            print minx+i*boxx/dp+random.uniform(-0.5,0.5)*epsilon,"\t",miny+random.uniform(-0.5,0.5)*epsilon,"\t",minz+j*boxz/dp+random.uniform(-0.5,0.5)*epsilon,"\t",1.0,"\t",1.0
            print minx+i*boxx/dp+random.uniform(-0.5,0.5)*epsilon,"\t",maxy+random.uniform(-0.5,0.5)*epsilon,"\t",minz+j*boxz/dp+random.uniform(-0.5,0.5)*epsilon,"\t",1.0,"\t",1.0
            #print "draw sphere {",minx+i*boxx/dp,miny,minz+j*boxz/dp,"} radius 1.0 resolution 10"
            #print "draw sphere {",minx+i*boxx/dp,maxy,minz+j*boxz/dp,"} radius 1.0 resolution 10"
    for i in range (1,dp):
        for j in range (1,dp):
            print minx+random.uniform(-0.5,0.5)*epsilon,"\t",miny+i*boxy/dp+random.uniform(-0.5,0.5)*epsilon,"\t",minz+j*boxz/dp+random.uniform(-0.5,0.5)*epsilon,"\t",1.0,"\t",1.0
            print maxx+random.uniform(-0.5,0.5)*epsilon,"\t",miny+i*boxy/dp+random.uniform(-0.5,0.5)*epsilon,"\t",minz+j*boxz/dp+random.uniform(-0.5,0.5)*epsilon,"\t",1.0,"\t",1.0
            #print "draw sphere {",minx,miny+i*boxy/dp,minz+j*boxz/dp,"} radius 1.0 resolution 10"
            #print "draw sphere {",maxx,miny+i*boxy/dp,minz+j*boxz/dp,"} radius 1.0 resolution 10"
def print_for_vmd():
    for i in range (0,dp+1):
        for j in range (0,dp+1):
          ##print minx+i*boxx/dp,"\t",miny+j*boxy/dp,"\t",maxz,"\t",1.0,"\t",1.0
          ##print minx+i*boxx/dp,"\t",miny+j*boxy/dp,"\t",minz,"\t",1.0,"\t",1.0
            print "draw sphere {",minx+random.uniform(-0.5,0.5)*epsilon+i*boxx/dp,miny+j*boxy/dp+random.uniform(-0.5,0.5)*epsilon,maxz+random.uniform(-0.5,0.5)*epsilon,"} radius 1.0 resolution 10"
            print "draw sphere {",minx+random.uniform(-0.5,0.5)*epsilon+i*boxx/dp,miny+j*boxy/dp+random.uniform(-0.5,0.5)*epsilon,minz+random.uniform(-0.5,0.5)*epsilon,"} radius 1.0 resolution 10"
    for i in range (0,dp+1):
        for j in range (1,dp):
          ##print minx+i*boxx/dp,"\t",miny,"\t",minz+j*boxz/dp,"\t",1.0,"\t",1.0
          ##print minx+i*boxx/dp,"\t",maxy,"\t",minz+j*boxz/dp,"\t",1.0,"\t",1.0
            print "draw sphere {",minx+random.uniform(-0.5,0.5)*epsilon+i*boxx/dp,miny+random.uniform(-0.5,0.5)*epsilon,minz+j*boxz/dp+random.uniform(-0.5,0.5)*epsilon,"} radius 1.0 resolution 10"
            print "draw sphere {",minx+random.uniform(-0.5,0.5)*epsilon+i*boxx/dp,maxy+random.uniform(-0.5,0.5)*epsilon,minz+j*boxz/dp+random.uniform(-0.5,0.5)*epsilon,"} radius 1.0 resolution 10"
    for i in range (1,dp):
        for j in range (1,dp):
          ##print minx,"\t",miny+i*boxy/dp,"\t",minz+j*boxz/dp,"\t",1.0,"\t",1.0
          ##print maxx,"\t",miny+i*boxy/dp,"\t",minz+j*boxz/dp,"\t",1.0,"\t",1.0
            print "draw sphere {",minx+random.uniform(-0.5,0.5)*epsilon,miny+i*boxy/dp+random.uniform(-0.5,0.5)*epsilon,minz+j*boxz/dp+random.uniform(-0.5,0.5)*epsilon,"} radius 1.0 resolution 10"
            print "draw sphere {",maxx+random.uniform(-0.5,0.5)*epsilon,miny+i*boxy/dp+random.uniform(-0.5,0.5)*epsilon,minz+j*boxz/dp+random.uniform(-0.5,0.5)*epsilon,"} radius 1.0 resolution 10" 
print_for_file()
#   print "draw line {",maxx,maxy,maxz,"} {",maxx,miny,maxz,"}"
#   print "draw line {",maxx,maxy,maxz,"} {",minx,maxy,maxz,"}"
#   print "draw line {",minx,miny,maxz,"} {",maxx,miny,maxz,"}"
#   print "draw line {",minx,miny,maxz,"} {",minx,maxy,maxz,"}"

#   print "draw line {",maxx,maxy,minz,"} {",maxx,miny,minz,"}"
#   print "draw line {",maxx,maxy,minz,"} {",minx,maxy,minz,"}"
#   print "draw line {",minx,miny,minz,"} {",maxx,miny,minz,"}"
#   print "draw line {",minx,miny,minz,"} {",minx,maxy,minz,"}"

#   print "draw line {",maxx,maxy,maxz,"} {",maxx,maxy,minz,"}"
#   print "draw line {",maxx,miny,maxz,"} {",maxx,miny,minz,"}"
#   print "draw line {",minx,maxy,maxz,"} {",minx,maxy,minz,"}"
#   print "draw line {",minx,miny,maxz,"} {",minx,miny,minz,"}"


