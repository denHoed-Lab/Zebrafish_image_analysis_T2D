import sys
#import matplotlib.pyplot as plt
import cv2
import Augment as Augmentation
import numpy as np
import glob, os
import torch
import torch.backends.cudnn as cudnn
import torch.nn as nn
from torch import optim
from modeling.deeplab import DeepLab
import csv
from skimage.morphology import skeletonize, binary_dilation, medial_axis
from scipy.spatial import distance as dist
def loadImage(imageFile):
    

    nrColorChannels = 3
    img = cv2.imread(imageFile, cv2.IMREAD_COLOR)
    img_width = img.shape[1]
    img_height = img.shape[0]
    img = img.astype(float)
    img = Augmentation.normalize_meanstd(img)
    img = Augmentation.normalize_float(img)
    #            img = cv2.resize(img, (self.img_width, self.img_height))
    img = np.transpose(img, (2,0,1))
    img = np.reshape(img, (1, nrColorChannels, img_height, img_width ))
    return img
def imshow(image):
    plt.imshow(image)
    plt.show()

def cvimshow(image):
    image = cv2.normalize(image,   0, 255, cv2.NORM_MINMAX)
    cv2.imshow('img', image)
    cv2.waitKey()

def midpoint(ptA, ptB):
	return ((ptA[0] + ptB[0]) * 0.5, (ptA[1] + ptB[1]) * 0.5)

def skeleton_endpoints(skel):
    # make out input nice, possibly necessary
    skel = skel.copy()
    skel[skel!=0] = 1
    skel = np.uint8(skel)

    # apply the convolution
    kernel = np.uint8([[1,  1, 1],
                       [1, 10, 1],
                       [1,  1, 1]])
    src_depth = -1
    filtered = cv2.filter2D(skel,src_depth,kernel)

    # now look through to find the value of 11
    # this returns a mask of the endpoints, but if you just want the coordinates, you could simply return np.where(filtered==11)
    out = np.zeros_like(skel)
    out[np.where(filtered==11)] = 1
    return out

def analyzeFolder(image_file_path,output_path):
    
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    net_weights = os.path.join(script_dir, 'weights', 'weights_body.pth')
    #print(net_weights)
    net = DeepLab(backbone='resnet', output_stride=16)
    device = 'cpu'
    #net.cuda()
    net.load_state_dict(torch.load(net_weights, map_location=torch.device('cpu')))
    #net.load_state_dict(torch.load(net_weights))
 
    
    net.eval()
 
    fileName = os.path.basename(image_file_path)
    print('Analysing image: '+fileName)
    fileName = ('').join(fileName.split('.')[0:-1])

    image = loadImage(image_file_path)
    image_torch = torch.from_numpy(image).float().to(device)
    masks_pred = net(image_torch)
    bw = np.argmax(masks_pred.cpu().data.numpy(), axis=1)
    bw = bw[0,:,:].astype(np.uint8)
  
        
    # Get all contours
    cnts, hierarchy = cv2.findContours(bw, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE) 
    cnt = sorted(cnts, key=cv2.contourArea)

    # Redraw results with only the largest contour
    segment_result = np.zeros_like(bw,np.uint8)
    cv2.fillPoly(segment_result,[cnt[-1]],1)

    # Get the contour of only the largest
    cnts, hierarchy = cv2.findContours(segment_result, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE) 
    cnt = sorted(cnts, key=cv2.contourArea)

    # Get the arera of the contour
    area = cv2.contourArea(cnt[-1])
    perimeter = cv2.arcLength(cnt[-1],False)
    # Read input image agian
    image_contour = cv2.imread(image_file_path, cv2.IMREAD_COLOR)
    
    # Get the skeleton
    skeleton = skeletonize(segment_result)
    skeleton = binary_dilation(skeleton)
    
    image_contourB=image_contour[:,:,0]
    image_contourB[skeleton]=255
    image_contour[:,:,0] = image_contourB

    # Box metric
    box = cv2.minAreaRect(cnt[-1])
    box = cv2.boxPoints(box)
    box = np.array(box, dtype="int")
    (tl, tr, br, bl) = box
    (tltrX, tltrY) = midpoint(tl, tr)
    (blbrX, blbrY) = midpoint(bl, br)
    (tlblX, tlblY) = midpoint(tl, bl)
    (trbrX, trbrY) = midpoint(tr, br)
    boxLength = max(dist.euclidean((tltrX, tltrY), (blbrX, blbrY)),dist.euclidean((tlblX, tlblY), (trbrX, trbrY)))
    boxWidth = min(dist.euclidean((tltrX, tltrY), (blbrX, blbrY)),dist.euclidean((tlblX, tlblY), (trbrX, trbrY)))

    # Draw contour on image
    cv2.drawContours(image_contour, [box.astype("int")], -1, (0, 255, 0), 2)
    cv2.drawContours(image_contour, cnt, -1, (0,0,255), 2)
    
    # Use proper path joining instead of manual string manipulation
    # cv2.imwrite(os.path.join(output_path, fileName+'.tif'), image_contour)  # Removed demo TIF file

    cv2.imwrite(os.path.join(output_path, fileName+'_seg.tif'), segment_result*255) 
    # print(output_path+fileName+'_seg.tif')


if __name__ == "__main__":

    image_file_path = sys.argv[1]
    output_path = sys.argv[2]
    
    # Clean the output path - remove any quotes or problematic characters
    output_path = output_path.strip().strip('"').strip("'")
    
    print("DEBUG: Received arguments:")
    print(f"  image_file_path: '{image_file_path}'")
    print(f"  output_path: '{output_path}'")
    print(f"  output_path exists: {os.path.exists(output_path)}")
    
    analyzeFolder(image_file_path, output_path)
 
  