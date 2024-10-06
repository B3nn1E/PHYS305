# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:35:04 2024

@author: kietb
"""
"""
# Area calculation code
# Computes area of a rectangle

print("This program calculates the area of a rectangle")

Length = float(input("What is the length of the rectangle?"))
Width = float(input("What is the width of the rectangle?"))

area = Length * Width

print(area)
"""

# Area calculation code
# Computes the area of a rectangle

print("This program calculates areas.")
print()

print("What shape would you like to know the area for?")
print("1 Rectangle")
print("2 Circle")

shape = int(input(">"))

if shape == 1:
    
    Length = float(input("What is the length of the rectangle?"))
    Width = float(input("What is the width of the rectangle?"))
    area = Length * Width
    print("The area is", area)
    
elif shape ==2:
    
    Radius = float(input("What is the radius of the circle?"))
    area = Radius**2*3.1415
    print("The area is", area)
    
else:
    print("ERROR: Wrong shape selected.")
    