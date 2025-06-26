import sys
import os
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
import modeling.Unet2D as unet2D
import csv
#from skimage.morphology import skeletonize, binary_dilation, medial_axis
 #from scipy.spatial import distance as dist
 #import matplotlib.pyplot as plt

orignal_width = 2048
orignal_height = 2048

def loadImage(imageFile):
    

    nrColorChannels = 1
 
    img = cv2.imread(imageFile, cv2.IMREAD_ANYDEPTH)
    img_width = 512
    img_height = 512
    img = img.astype(float)
    img = Augmentation.normalize_meanstd_gray(img)
    img = Augmentation.normalize_float(img)
    
    img = cv2.resize(img, (img_width, img_height))
    img = np.reshape(img, (1, nrColorChannels, img_height, img_width ))
    return img
    
    
def imshow(image):
    plt.imshow(image)
    plt.show()
    

def cvimshow(image):
    image = cv2.normalize(image,   0, 255, cv2.NORM_MINMAX)
    cv2.imshow('img', image)
    cv2.waitKey()


def analyzeFolder(image_file_path, output_path):
    
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # NOTE: These are different neural network snapshots
    net_weights = os.path.join(script_dir, "weights", "weights_liver.pth")
    #net_weights = "./weights/Liver_weights_v1.00_20210211.pth" # change back when done
    
    in_channels = 1
    num_classes = 2
    net = unet2D.Unet(in_channels, num_classes)
    #print(net_weights)
    
 
    # if torch.cuda.is_available():
    #     net.load_state_dict(torch.load(net_weights))
    # else:
    #     net.load_state_dict(torch.load(net_weights, map_location={'cuda:0': 'cpu'}))
    # net.load_state_dict(torch.load(net_weights, map_location=torch.device('cpu')))
    
    #device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    device = 'cpu'
    # net.load_state_dict(torch.load(net_weights,map_location='cuda:0'))
    # net.cuda()
    net.load_state_dict(torch.load(net_weights, map_location=torch.device('cpu')))
 
    
    net.eval()
 
    fileName = os.path.basename(image_file_path)
    print('Analysing image: '+fileName)
    fileName = ('').join(fileName.split('.')[0:-1])

    image = loadImage(image_file_path)
    image_torch = torch.from_numpy(image).float().to(device)
    masks_pred = net(image_torch)
    bw = np.argmax(masks_pred.cpu().data.numpy(), axis=1)
    bw = bw[0,:,:].astype(np.uint8)
    bw = cv2.resize(bw, (orignal_width, orignal_height))
    
    
        
    # Get all contours
    cnts, hierarchy = cv2.findContours(bw, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE) 
    cnt = sorted(cnts, key=cv2.contourArea)
    
    # Redraw results with only the largest contour
    segment_result = np.zeros_like(bw,np.uint8)
    for i in range(0,len(cnt)):
        cv2.fillPoly(segment_result,[cnt[i]],1)
    
    image_in = cv2.imread(image_file_path, cv2.IMREAD_ANYDEPTH)
    image_out = Augmentation.normalize_float(image_in)*255
    # Get the skeleton
     
    image_out = image_out.astype(np.uint8)
    image_out = cv2.cvtColor(image_out, cv2.COLOR_GRAY2BGR)
    

    # Draw contour on image
    cv2.drawContours(image_out, cnt, -1, (0,0,255), 2)
    cv2.imwrite(os.path.join(output_path, fileName+'.tif'), image_out)
    print(os.path.join(output_path, fileName+'.tif'))
    cv2.imwrite(os.path.join(output_path, fileName+ '_seg.tif'), segment_result*255)
    

    # NOTE: Set this to True to write the debug images, False otherwise
    write_debug = False
    debug_path = r"D:/Pipeline_2022/debug_images/"
    if write_debug and os.path.exists(debug_path):
        cv2.imwrite(debug_path+fileName+ '_img.tif', image_in) 
        cv2.imwrite(debug_path+fileName+ '_seg.tif', segment_result*255) 

    

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
 
  