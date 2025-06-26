import torch
from torchvision import transforms
import cv2
import numpy as np
import types
from numpy import random

def subtractMean(image):

    mean = np.mean(image, axis=(0, 1))
    image[:,:,0]-=mean[0]
    image[:, :, 1] -= mean[1]
    image[:, :, 2] -= mean[2]
    return image.astype(float)

def shrinkExp(image, mask,rx,ry):
    # rx = random.uniform(0.7, 1.1)
    # ry = random.uniform(0.7, 1.1)
    imsize = image.shape
    M = np.float32([[rx, 0, 0], [0, ry, 0]])
    image = cv2.warpAffine(image, M, (imsize[1],imsize[0]))
    mask = cv2.warpAffine(mask, M, (imsize[1], imsize[0]))
    return image,mask

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