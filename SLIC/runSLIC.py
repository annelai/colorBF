import os
import subprocess


spnumber = 1000
bmpdir = './color_bmp/'

# create a folder for dat files
spdir = './superpixel_%d/' % spnumber
os.mkdir(spdir)

# get all bmp files
for files in os.listdir(bmpdir):
    if files.endswith(".bmp"):
        bmp_name = bmpdir + files
        cmd = 'SLICSuperpixelSegmentation %s 20 %s %s' % (bmp_name, spnumber, spdir)
        subprocess.call(cmd, shell=True)
