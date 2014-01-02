# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <rawcell>

# Demonstration of OpenCV's cv2 API for manipulating images, creating contours, simplifying those contours and creating various bounding structures.  I started with Abid Rahman's Contours 1 and 2 blog posts.  Those posts do not directly give you all that you need to run them.  By virtue of being an IPython Notebook, these are the actual command that worked for me.  You will need wget, numpy, matplotlib, OpenCV and shapely.
# 
# -kurt schwehr 2013-Jan-31

# <codecell>

# BEGIN http://www.opencvpython.blogspot.com/2012/06/hi-this-article-is-tutorial-which-try.html

# <codecell>

!wget http://3.bp.blogspot.com/-a5blM3JLkIU/T9OXg1YhN0I/AAAAAAAAASQ/MbdfSG2oaYg/s200/test.jpg

# <codecell>

import numpy as np
import cv2
import cv
import shapely.geometry

# <codecell>

im = cv2.imread('test.jpg')

# <codecell>

whos

# <codecell>

imgray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
imshow(imgray)

# <codecell>

ret, thresh = cv2.threshold(imgray, 127, 255, 0)
imshow(thresh)

# <codecell>

contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
imshow(thresh)

# <codecell>

len(contours), hierarchy

# <codecell>

len(contours[0])

# <codecell>

cv2.drawContours(im, contours, -1, (0,255,0), 3)
imshow(im)

# <codecell>

cv2.drawContours(im, contours, -1, (0,255,0), -1)
imshow(im)

# <codecell>

cv2.drawContours(im, contours, 0, (0,255,0), 1)
imshow(im)

# <codecell>

!wget http://3.bp.blogspot.com/-1UtLXb7c73U/T9QZT3tpVjI/AAAAAAAAATE/Nyo7SFg8T1o/s1600/balls.png

# <codecell>

im = cv2.imread('balls.png')
imgray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)

# <codecell>

ret, thresh = cv2.threshold(imgray, 100, 255, 0)
imshow(thresh)

# <codecell>

contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

# <codecell>

for h, cnt in enumerate(contours):
    mask = np.zeros(imgray.shape, np.uint8)
    cv2.drawContours(mask, [cnt], 0, 255, -1)
    mean = cv2.mean(im, mask=mask)
    print h, mean

# <codecell>

# END http://www.opencvpython.blogspot.com/2012/06/hi-this-article-is-tutorial-which-try.html

# <headingcell level=1>

# BEGIN Contours - 2

# <codecell>

moments = cv2.moments(contours[0])
area = moments['m00']
moments

# <codecell>

cv2.contourArea(contours[0])

# <codecell>

perimeter = cv2.arcLength(contours[0], True) # True says that curve is closed
perimeter

# <codecell>

cnt = contours[0]
len(cnt)

# <codecell>

x = cnt[:,0, 0]; y = cnt[:,0,1]
plot(x,y)

# <codecell>

simpler = cv2.approxPolyDP(cnt, 2, True)
plot(simpler[:,0,0], simpler[:,0,1])

# <codecell>

hull = cv2.convexHull(cnt)
plot(hull[:,0,0], hull[:,0,1])

# <codecell>

simpler_hull = cv2.approxPolyDP(hull, 2.5, True)
plot(simpler_hull[:,0,0], simpler_hull[:,0,1])
plot(x,y)

# <codecell>

r_x, r_y, r_w, r_h = cv2.boundingRect(cnt)
r_x, r_y, r_w, r_h  # (150, 121, 103, 146)
plot((r_x, r_x, r_x+r_w, r_x+r_w, r_x), (r_y, r_y+r_h, r_y+r_h, r_y, r_y))
plot(x,y)

# <codecell>

rect = cv2.minAreaRect(cnt)
rect # ((202.134521484375, 192.14178466796875), (102.39618682861328, 140.3079376220703), -5.128190994262695)
box = cv2.cv.BoxPoints(rect)
box # ((157.41201782226562, 266.59124755859375), (144.87069702148438, 126.84494018554688), (246.85702514648438, 117.69232177734375), (259.3983459472656, 257.4386291503906))
# plot( [p[0] for p in box] + [box[0][0]], [p[1] for p in box] + [box[0][1]] )
box_list = list(box)
box_list.append(box[0])
ba = np.array(box_list) # Box array
plot(ba[:,0], ba[:,1])
plot(x,y)

# <codecell>

(c_x, c_y), radius = cv2.minEnclosingCircle(cnt)
c_x, c_y, radius # (197.0, 194.5, 82.92139434814453)
center = shapely.geometry.Point(c_x, c_y)
circle = np.array(center.buffer(radius).boundary.coords)
len(circle) # 66 points
plot(circle[:,0], circle[:,1])
plot(x,y)

# <codecell>

ellipse = cv2.fitEllipse(cnt)
ellipse # ((199.31251525878906, 185.9192352294922), (93.7149658203125, 138.58531188964844), 202.948486328125)
# TODO: what is a convenient way to get the coords for an ellipse?

# <codecell>


