import torch
from torchvision import transforms
import cv2
import numpy as np
import types
from numpy import random
import gryds as gryds
import utils.implot  as imp
from scipy import ndimage
import scipy.misc
import scipy.ndimage
from random import randint

def subtractMean(image):

    mean = np.mean(image, axis=(0, 1))
    image[:,:,0]-=mean[0]
    image[:, :, 1] -= mean[1]
    image[:, :, 2] -= mean[2]
    return image.astype(float)
def randomFlip(image, mask):
    flip =  flip = randint(-1, 2)
    image = cv2.flip(image,flip )
    mask = cv2.flip(mask,flip )
    return image, mask
def shrinkExp(image, mask,rx,ry):
    # rx = random.uniform(0.7, 1.1)
    # ry = random.uniform(0.7, 1.1)
    imsize = image.shape
    M = np.float32([[rx, 0, 0], [0, ry, 0]])
    image = cv2.warpAffine(image, M, (imsize[1],imsize[0]))
    mask = (np.round(cv2.warpAffine(mask.astype(float), M, (imsize[1], imsize[0]),flags=cv2.INTER_NEAREST))).astype(int)
    return image,mask

def rotate(image, mask,maxminrot=(-90,90)):
    rot = random.uniform(maxminrot[0], maxminrot[1])
    #img = scipy.misc.imrotate(img, rot, interp='bilinear').astype(float)/255
    image = scipy.ndimage.rotate(image,rot,reshape=False,mode='reflect').astype(float)
    mask =scipy.ndimage.rotate(mask,rot,reshape=False,order=0,mode='constant',cval=0) 

    return image, mask
def bspline(image,mask):
    random_grid = np.random.rand(2, 3, 3)
    random_grid -= 0.5
    random_grid /= random.uniform(10, 20)
    bspline = gryds.BSplineTransformation(random_grid)
    if len(image.shape)>2:

        for ch in range(0,image.shape[2]):
            interpolator = gryds.Interpolator(image[:,:,ch])
            image[:,:,ch] = interpolator.transform(bspline)
    else:
        interpolator = gryds.Interpolator(image)
        image = interpolator.transform(bspline)

    mask2 = np.zeros((mask.shape))
    for label in range(1,mask.max()+1):
        mask_interpolator = gryds.Interpolator((mask==label).astype(int))
        mask2[np.round(mask_interpolator.transform(bspline,order=0))==1]=label
    mask = mask2
    return image, mask


def normalize_meanstd_gray(image):
    # axis param denotes axes along which mean & std reductions are to be performed

    mean = np.mean(image, axis=(0, 1))
    std = np.std(image, axis=(0, 1))

    image[:, :] = (image[:, :] - mean) / std
    #image[:, :, 1] = (image[:, :, 1] - mean[1]) / std[1]
    #image[:, :, 2] = (image[:, :, 2] - mean[2]) / std[2]
    return image.astype(float)

def normalize_meanstd(image):
    # axis param denotes axes along which mean & std reductions are to be performed

    mean = np.mean(image, axis=(0, 1))
    std = np.std(image, axis=(0, 1))

    for ch in range(0,image.shape[2]):
        image[:, :, ch] = (image[:, :, ch] - mean[ch]) / std[ch]
    #image[:, :, 1] = (image[:, :, 1] - mean[1]) / std[1]
    #image[:, :, 2] = (image[:, :, 2] - mean[2]) / std[2]
    return image.astype(float)

def normalize_meanstd_3D(image):
    # axis param denotes axes along which mean & std reductions are to be performed


    mean =  np.mean(image)
    std =  np.std(image)
    image  = (image - mean) / std
 
    #image[:, :, 1] = (image[:, :, 1] - mean[1]) / std[1]
    #image[:, :, 2] = (image[:, :, 2] - mean[2]) / std[2]
    return image.astype(float)

def random_RGBshift(image,scale=0.1):
    r_shift = random.uniform(-0.2*scale, 0.2*scale)
    g_shift = random.uniform(-0.2 * scale, 0.2 * scale)
    b_shift = random.uniform(-0.2 * scale, 0.2 * scale)

    image[..., 0] += r_shift
    image[..., 1] += g_shift
    image[..., 2] += b_shift
    return image

def normalize_float(image):
    image = (image - image.min()) / (image.max() - image.min())
    return image


def random_gamma(image):
    # build a lookup table mapping the pixel values [0, 255] to
    # their adjusted gamma values
    gamma = random.uniform(0.5,2)
    invGamma = 1.0 / gamma
    table = np.array([((i / 255.0) ** invGamma) * 255
                      for i in np.arange(0, 256)]).astype("uint8")
    image = cv2.LUT(image, table).reshape(image.shape)
    # apply gamma correction using the lookup table
    return image