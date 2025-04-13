"""
pip install opencv-python
"""
import cv2
import numpy as np

jim = cv2.imread("0.jpg")
im = cv2.imread("0.png")
# Use this to if you read png image
im = cv2.cvtColor(im,cv2.COLOR_BGRA2RGB)
# Use this to if you read jpg image
jim = cv2.cvtColor(im,cv2.COLOR_BGR2YCrCb)

from dct import *

print("conv_check")
print("ORGB",list(im[0,0]))
print("CRGB",YCbCrtoRGB(list(jim[0,0])))
print("OYCbCr",list(jim[0,0]))
print("CYCbCr",RGBtoYCbCr(list(im[0,0])))


ty,tx = 0,150

ojim = jim
oim = im
jim = jim.tolist()
YCbCr = subsampl_4_2_0(jim)



jim = [[
    [
        YCbCr[0][i][k],
        YCbCr[2][round((i/2)-0.1)][round((k/2)-0.1)],
        YCbCr[1][round((i/2)-0.1)][round((k/2)-0.1)]
        ]
    for k in range(len(YCbCr[0][0]))]for i in range(len(YCbCr[0]))]

y1,y2=0,511
x1,x2=0,511
for i in range(y1,y2):
    for k in range(x1,x2):
        im[i,k] = np.array(YCbCrtoRGB(jim[i][k]))


print("YCbCr after subs ",jim[0][7])
print("RGB after subs",YCbCrtoRGB(jim[ty][tx]))
print("origin YCbCr",ojim[ty,tx])
print("origin RGB",oim[ty,tx])


img = np.array(im[y1:y2,x1:x2], np.uint8)
cv2.imwrite('output_image.jpg', img,)

    