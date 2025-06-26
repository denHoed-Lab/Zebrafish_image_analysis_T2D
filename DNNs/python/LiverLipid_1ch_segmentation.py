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
#from skimage.morphology import skeletonize, binary_dilation, medial_axis
#from scipy.spatial import distance as dist
import tifffile as tiff
# import utils.implot as imp
def loadImage(imageFile):
    

    nrColorChannels = 1
    img = tiff.imread(imageFile)
    img = img[:,:,None]
    img_width = 512
    img_height = 512
    img = img.astype(float)
    img = Augmentation.normalize_meanstd(img)
    for ch in range(0,img.shape[2]):
        img[:,:,ch] = Augmentation.normalize_float(img[:,:,ch])
    img= cv2.resize(img, ( img_width, img_height))
    img = img[:,:,None]
    img = np.reshape(img, (1, nrColorChannels, img_height, img_width ))
    return img
# def imshow(image):
#     plt.imshow(image)
#     plt.show()

def cvimshow(image):
    image = cv2.normalize(image,   0, 255, cv2.NORM_MINMAX)
    cv2.imshow('img', image)
    cv2.waitKey()

 
def analyzeFolder(image_file_path,output_path):
    
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    net_weights = os.path.join(script_dir, 'weights', 'weights_lipid.pth')
    net = DeepLab(backbone='resnet', output_stride=16,num_classes=3,in_channels=1)

    device = 'cpu'
    net.load_state_dict(torch.load(net_weights, map_location=torch.device('cpu')))

    net.eval()
 
    fileName = os.path.basename(image_file_path)
    print('Analysing image: '+fileName)
    fileName = ('').join(fileName.split('.')[0:-1])

    image = loadImage(image_file_path)
    image_torch = torch.from_numpy(image).float().to(device)
    masks_pred = net(image_torch)
    pred = masks_pred.argmax(1)#
    bw = np.uint8(pred.detach().cpu().numpy())[0,:,:]

    # Get all contours
    cnts, hierarchy = cv2.findContours((bw==1).astype(np.uint8), cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE) 
    cnt = sorted(cnts, key=cv2.contourArea)

    # Redraw results with only the largest contour
    roi_filled = np.zeros_like(bw,np.uint8)
    cv2.fillPoly(roi_filled,[cnt[-1]],1)
    
    roi_filled[bw==2]=2
    bw = roi_filled

    cv2.imwrite(os.path.join(output_path, fileName+'_seg.png'), bw)
    #print(output_path+fileName+'_seg.tif')

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
    
 
  