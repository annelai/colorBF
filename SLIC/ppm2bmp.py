import cv2
import os

ppmdir = './color_ppm/'
bmpdir = './color_bmp/'

# create a folder to save bmp files
os.mkdir(bmpdir)

# get all ppm files
for files in os.listdir(ppmdir):
    if files.endswith(".ppm"):
        ppm_name = ppmdir + files
        img = cv2.imread(ppm_name)
        bmp_name = bmpdir + files.split('.')[0] + '.bmp'
        cv2.imwrite(bmp_name, img)
