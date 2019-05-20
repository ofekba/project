import sys
import math 

x1 = float (input("Enter x1: ") )
x2 = float (input("Enter x2: ") )
x3 = float (input("Enter x3: ") )
y1 = float (input("Enter x1: ") )
y2 = float (input("Enter x2: ") )
y3 = float (input("Enter x3: ") )

r2=pow(y1-x1,2)+pow(y2-x2,2)+pow(y3-x3,2)
r=math.sqrt(r2)
print("the distance is: ",r)

