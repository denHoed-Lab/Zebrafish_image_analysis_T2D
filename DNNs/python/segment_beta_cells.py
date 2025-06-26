import sys
#import matplotlib.pyplot as plt
import cv2
import utils.Augment as Augmentation
import numpy as np
import glob, os
import torch
import torch.backends.cudnn as cudnn
import torch.nn as nn
from torch import optim
import tifffile
from modeling.Unet3D import UNet3D
import csv
from skimage.morphology import skeletonize, binary_dilation, medial_axis
from scipy.spatial import distance as dist
from scipy import ndimage as ndi
def loadImage(imageFile):
    
    nrColorChannels = 1
    z_size = 128
    img  = cv2.imreadmulti(imageFile)
    img = np.dstack(img[1])
    #img =img[:,:,20:]
    img_original = img
    img_size = img.shape
    input_resolution=[0.419, 0.419, 1.143]; # micrometer , z stack beta cells data
    output_resolution=1.143
    input_resolution[0]/output_resolution
    final_image_size = [round(x/output_resolution*img_size[ii]) for ii,x in enumerate(input_resolution)]
    final_image_size = [224,224,final_image_size[2]]
    img = cv2.resize(img, (final_image_size[0], final_image_size[1]))
 
 
    img = img.astype(float)
    img = Augmentation.normalize_meanstd_3D(img)
    img = Augmentation.normalize_float(img)
    if img.shape[2]>z_size:
        diff = img.shape[2]-z_size
        img = img[:,:,diff:]
    if img.shape[2]<z_size:

        img = np.pad(img,( (0, 0), (0, 0), ( z_size-img.shape[2],0)))
 
    #            img = cv2.resize(img, (self.img_width, self.img_height))
     
    img = np.reshape(img, (1, nrColorChannels, final_image_size[0], final_image_size[1],z_size ))
    return img, img_original


def analyzeFolder(image_file_path,output_path):
    
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    net_weights = os.path.join(script_dir, 'weights', 'weights_betacells.pth')
 
    #print(net_weights)
    
    net = UNet3D(in_channels=1, n_classes=2, base_n_filter=3)
    # if torch.cuda.is_available():
    #     net.load_state_dict(torch.load(net_weights))
    # else:
    #     net.load_state_dict(torch.load(net_weights, map_location={'cuda:0': 'cpu'}))
    # net.load_state_dict(torch.load(net_weights, map_location=torch.device('cpu')))
    
    #device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    #net.load_state_dict(torch.load(net_weights,map_location='cuda:0'))
    device = 'cpu'
    #net.cuda()
    net.load_state_dict(torch.load(net_weights, map_location=torch.device('cpu')))
    #net.load_state_dict(torch.load(net_weights))
 
    
    net.eval()
 
    fileName = os.path.basename(image_file_path)
    print('Analysing image: '+fileName)
    fileName = ('').join(fileName.split('.')[0:-1])

    image, img_original = loadImage(image_file_path)
    img_original = img_original/np.max(img_original)
    image_torch = torch.from_numpy(image).float().to(device)
    masks_pred = net(image_torch)
    bw = np.argmax(masks_pred.cpu().data.numpy(), axis=1).astype(np.uint8)
    bw = cv2.resize(bw[0,:,:,:], (img_original.shape[0], img_original.shape[1]),interpolation = cv2.INTER_NEAREST)
    if img_original.shape[2]>bw.shape[2]:
        bw = np.pad(bw,( (0, 0), (0, 0), ( img_original.shape[2]-bw.shape[2],0)))
    elif img_original.shape[2]<bw.shape[2]:
        bw = bw[:,:,bw.shape[2]-img_original.shape[2]:]

    #Irgb = np.zeros((img_original.shape[0],img_original.shape[1],3,img_original.shape[2])) 
    #for i in range(0,img_original.shape[2]):
    #           Irgb[:,:,0,i] = cv2.max(img_original[:,:,i].astype(float),(bw[:,:,i]==1).astype(float))
     #           Irgb[:,:,1,i] = img_original[:,:,i]
     #           Irgb[:,:,2,i] = img_original[:,:,i]
    #with tifffile.TiffWriter(output_path+fileName+'_seg.tif') as stack: 
     #   for idx in range(0,Irgb.shape[3]):
     #       stack.save((Irgb[:,:,:,idx]*255).astype('uint8'),compress='DEFLATE')

    print("DEBUG")
    print("path exists:", os.path.exists(output_path))
    print(os.path.join(output_path, fileName+'_bw.tif'))

    # with open(os.path.join(output_path, fileName+'_bw.txt'), "w") as f:  # Removed txt file creation
    #     f.write("test")  # Removed txt file creation
    
    with tifffile.TiffWriter(os.path.join(output_path, fileName+'_bw.tif')) as stack: 
        for idx in range(0,bw.shape[2]):
            stack.save((bw[:,:,idx]*255).astype('uint8'),compression='DEFLATE')
    #2024-04-10:We changed in line 99 "compress" to "compression" because the package has been updated
    #image_flat = img_original.max(2)
    #bw_flat_edge = ndi.morphology.binary_dilation(bw.max(2))-bw.max(2)
    #Irgb_flat = np.zeros((img_original.shape[0],img_original.shape[1],3)) 

    #Irgb_flat[:,:,2] = cv2.max(image_flat,bw_flat_edge.astype(float))
    #Irgb_flat[:,:,1] = image_flat
    #Irgb_flat[:,:,0] = image_flat
    #cv2.imwrite(output_path+'maxproj'+fileName+'.png',(Irgb_flat*255).astype(int))
  #  image_flat = img_original.max(2)
  #  bw_flat_edge = ndi.morphology.binary_dilation(bw.max(2))-bw.max(2)
  #  Irgb_flat = np.zeros((img_original.shape[0],img_original.shape[1],3)) 

 #   Irgb_flat[:,:,2] = cv2.max(image_flat,bw_flat_edge.astype(float))
 #   Irgb_flat[:,:,1] = image_flat
 #   Irgb_flat[:,:,0] = image_flat
 #   cv2.imwrite(output_path+'maxproj'+fileName+'.png',(Irgb_flat*255).astype(int))
    
# cv2.imwrite(output_path+fileName+'.tif', image_contour)

    # cv2.imwrite(output_path+fileName+'_seg.tif', segment_result*255) 

 

    

if __name__ == "__main__":
    #output_path = 'C:/Users/amial914/work/src/DL_seg_Novo/test_data/'
    #image_file_path = r"C:\Users\amial914\work\src\DL_seg_Novo\test_data\InsH2Bfabp_CD15_10dpf_VAST02_plate01_ACVR1C_201019_E04.tif"
    image_file_path = sys.argv[1]
    output_path = sys.argv[2]
    
    # Clean the output path - remove any quotes or problematic characters
    output_path = output_path.strip().strip('"').strip("'")
    
    print("DEBUG: Received arguments:")
    print(f"  image_file_path: '{image_file_path}'")
    print(f"  output_path: '{output_path}'")
    print(f"  output_path exists: {os.path.exists(output_path)}")
    
    analyzeFolder(image_file_path, output_path)
 
  