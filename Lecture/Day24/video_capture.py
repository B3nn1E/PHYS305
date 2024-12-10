# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 15:51:42 2024

@author: kietb
"""

import cv2
import numpy as np
import matplotlib.pyplot as plt

VidObj = cv2.VideoCapture('Rolling_can.mp4')
#VidObj.get(5)
ret,frame = VidObj.read()
print(frame[39,27,:])
frame[:,:,:] = frame[:,:,::-1]
cv2.imshow('image',frame)
plt.imshow(frame)

frame = frame[:,:,[2,1,0]]
plt.imshow(frame)
points = plt.ginput(2)
points
d = np.sqrt((points[1,0] - points[0,0]**2) + (points[1,1] - points[1,0])**2)
Pix2cm = 10.16/d 
red =frame
red[:,:,1] = 0 # sets all green to zero
red[:,:,2] = 0 # sets all blue to zero
plt.imshow(red)
gray = np.mean(frame, axis = 2, dtype = 'float')

plt.imshpw(gray, cmap = 'gray')
BminusGray = frame[:,:,2] - gray
plt.imshpw(BminusGray, cmap = 'gray')



