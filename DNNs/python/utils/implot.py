import matplotlib.pyplot as plt
import cv2
def imshow(img):
    plt.imshow(img)
    plt.draw()
    plt.show()

 
def imshowoverlay(img1,img2):
    plt.imshow(img1, cmap='gray')
    plt.imshow(img2, cmap='jet', alpha=0.5)  
    plt.show()


 

import numpy as np
import matplotlib.pyplot as plt


class IndexTracker(object):
    def __init__(self, ax, X):
        self.ax = ax
        ax.set_title('use scroll wheel to navigate images')

        self.X = X
        rows, cols, self.slices = X.shape
        self.ind = self.slices//2

        self.im = ax.imshow(self.X[:, :, self.ind])
        self.update()

    def onscroll(self, event):
        
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        self.im.set_data(self.X[:, :, self.ind])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()

class IndexTrackerRGB(object):
    def __init__(self, ax, X):
        self.ax = ax
        ax.set_title('use scroll wheel to navigate images')

        self.X = X
        rows, cols, colc, self.slices = X.shape
        self.ind = self.slices//2

        self.im = ax.imshow(self.X[:, :, :,self.ind])
        self.update()

    def onscroll(self, event):
        
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        self.im.set_data(self.X[:, :, :,self.ind])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()
def imshow3DOverlay(X,M=None):
    # mainN = drawContour(X,M,1,(255,0,0))
    # 
    Irgb = np.zeros((M.shape[0],M.shape[1],3,M.shape[2])) 
    for i in range(0,M.shape[2]):
        Irgb[:,:,:,i] = drawContourCV(X[:,:,i],M[:,:,i])
        
    # if not M is None:
    #     X=np.concatenate((X,M),axis=1)
    # X=np.concatenate((X,M),axis=1)
    fig, ax = plt.subplots(1, 1)

# for i in range(0,M.shape[0])
    
    tracker = IndexTrackerRGB(ax, Irgb)


    fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
    plt.show()

def imshow3D(X):
    # mainN = drawContour(X,M,1,(255,0,0))
    # 

        
    # if not M is None:
    #     X=np.concatenate((X,M),axis=1)
    # X=np.concatenate((X,M),axis=1)
    fig, ax = plt.subplots(1, 1)

# for i in range(0,M.shape[0])
    
    tracker = IndexTracker(ax, X)


    fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
    plt.show()

 
 
def drawContourCV(img,seg,color=(1,0,0)):
    img = img.astype('float32')
    if len(img.shape)<3:
        img = cv2.cvtColor(img.astype('float32') ,cv2.COLOR_GRAY2BGR)

    # Dictionary giving RGB colour for label (segment label) - label 1 in red, label 2 in yellow
    RGBforLabel = { 1:(1,0,0), 2:(0,1,1) }

    # Find external contours,
    edged = cv2.Canny(seg.astype('uint8'), 0, 0)
    contours,_ = cv2.findContours(edged,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE)
    cv2.drawContours(img,contours,-1,color,1)
    return img

def drawLabelCV(img,seg):
    R=img[:,:,0]
    G=img[:,:,1]
    B=img[:,:,2]
    R[seg]=0
    G[seg]=0
    B[seg]=1
    img[:,:,0]=R
    img[:,:,1]=G
    img[:,:,2]=B
 
    return img
    