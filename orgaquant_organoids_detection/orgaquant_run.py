###################################### 
### The code is modified from orgaquant.py from orgaquant package
### Set up this orgaquant script
### 1. Follow the instruction to install orgaquant in a virtual environment: https://github.com/TKassis/OrgaQuant
### 2. Create directory Orgaquant/plugin/, put statistics.py and visualization.py in it
### 3. activate orgaquant virtual environment
### 4. Run orgaquant_run.py instead of orgaquant.py
######################################

import tensorflow as tf
from tensorflow import keras
import seaborn as sns
from keras_retinanet import models
from keras_retinanet.utils.image import preprocess_image, resize_image, adjust_contrast
from plugin.statistics import area_ellipse
from plugin.visualization import nice_annotation, my_histogram
import cv2
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
import time 


folder_path = 'PATH_TO_DIR_WITH_INPUT_BRIGHT_FIELD_PICTURES'
output_path = 'PATH_TO_OUTPUT_DIR'
model_path = "PATH_TO_OrgaQuant/trained_models/orgaquant_intestinal_v3.h5"

# Larger "Image Size" allows you to detect smaller orgaoids at the cost of computational demand
min_side = 2048 
max_side = 2048
# 'Larger "Contrast" can improve detection sometimes
contrast = 2.5
# Use larger "Threshold" to eliminate false positives
threshold = 0.9
# How many pixels represents 1 micrometer
pixel_to_micrometer_rate = 0.6153



def load_orga_model():
    model = models.load_model(model_path, backbone_name='resnet50')
    return model

model = load_orga_model()

imagelist = []
prefixlist = []
outlist = []
outroot, _, _ = next(os.walk(output_path), [None, None, None]) # return Nones if not found
if not outroot: # if the output file not exist
    os.mkdir(output_path)
    outroot, _, _ = next(os.walk(output_path))
else:
    pass 

for root, directories, filenames in os.walk(folder_path):
    prefixlist = prefixlist + [x.rstrip('.jpg').rstrip('.tif').rstrip('.TIF').rstrip('.png').rstrip('.jpeg').rstrip('.tiff') for x in filenames if x.endswith(
        ('.jpg', '.tif', '.TIF', '.png', '.jpeg', '.tiff'))]
    imagelist = imagelist + [os.path.join(root,x) for x in filenames if x.endswith(('.jpg','.tif','.TIF', '.png', '.jpeg', '.tiff'))]
    outlist = outlist + [os.path.join(outroot,x) for x in filenames if x.endswith(('.jpg','.tif','.TIF', '.png', '.jpeg', '.tiff'))]

final_output = pd.DataFrame(np.empty((0,5)), 
                            columns=['x1', 'y1', 'x2', 'y2', 'filename'])
final_output = final_output.astype(dtype={'x1': int, 'y1': int, 'x2': int, 'y2': int, 'filename': str})


for i, filename in enumerate(imagelist):
    try:
        #IMAGE_PATH = os.path.join(root,filename)
        IMAGE_PATH = filename
        OUTPUT_PATH = outlist[i]
        PREFIX = prefixlist[i]
        # load image
        # image = read_image_bgr(IMAGE_PATH)
        image = cv2.imread(IMAGE_PATH)

        # copy to draw on
        draw = image.copy()
        draw = cv2.cvtColor(draw, cv2.COLOR_BGR2RGB)

        # preprocess image for network
        image = adjust_contrast(image,contrast)
        image = preprocess_image(image)
        image, scale = resize_image(image, min_side=min_side, max_side=max_side)

        # process image
        boxes, scores, labels = model.predict_on_batch(np.expand_dims(image, axis=0))

        # correct for image scale
        boxes /= scale

        out = np.empty((0,4), dtype=np.float32)

        # visualize detections
        for box, score, label in zip(boxes[0], scores[0], labels[0]):
            # exlude low score detection
            if score < threshold:
                continue
            
            # exclude near-border detection
            if np.mean([box[0],box[2]]) < min_side * 0.05 or \
                np.mean([box[0],box[2]]) > min_side * 0.95 or \
                np.mean([box[1],box[3]]) > min_side * 0.95 or \
                np.mean([box[1],box[3]]) < min_side * 0.05 :
                annot_color = (128, 128, 0)
                close_to_border = True
                continue
                
            else:
                annot_color = (247, 37, 133)
                close_to_border = False
                

            b = box.astype(int)
            box_micrometre = box/pixel_to_micrometer_rate

            nice_annotation(draw, b, annot_color,  str(area_ellipse(diameter1=box_micrometre[2] - box_micrometre[0],
                                                                        diameter2=box_micrometre[3] - box_micrometre[1])))
            out = np.append(out, box_micrometre.reshape(1,4), axis=0)
            final_output = final_output.append({'x1': box_micrometre[0], 'y1': box_micrometre[1],
                                            'x2': box_micrometre[2], 'y2': box_micrometre[3], 
                                            'filename': IMAGE_PATH, 'threshold': threshold,
                                            'close_to_border': close_to_border}, 
                                            ignore_index=True)

        output = pd.DataFrame(out, columns=['x1', 'y1', 'x2', 'y2'], dtype=np.int16)
        output['Diameter 1 (μm)'] = output['x2'] - output['x1']
        output['Diameter 2 (μm)'] = output['y2'] - output['y1']
        output['Ellipse_Area'] = area_ellipse(diameter1=output['Diameter 1 (μm)'],
                                                diameter2=output['Diameter 2 (μm)'])
        output['threshold'] = round(threshold,2)
        output.to_csv('{}_{}.csv'.format(OUTPUT_PATH, round(threshold,2)), index=False)
        plt.imsave('{}_{}_detected.pdf'.format(OUTPUT_PATH, round(threshold,2)), draw)
        fig = my_histogram(df=output)
        fig.savefig('{}_{}_histogram.pdf'.format(OUTPUT_PATH, round(threshold,2)), bbox_inches="tight")
    except Exception as e:
        print('Error: {}'.format(e))
        logging.error('Error: {}'.format(e))


final_output['Diameter 1 (μm)'] = final_output['x2'] - final_output['x1']
final_output['Diameter 2 (μm)'] = final_output['y2'] - final_output['y1']
final_output['Ellipse_Area (μm^2)'] = area_ellipse(diameter1=final_output['Diameter 1 (μm)'],
                                        diameter2=final_output['Diameter 2 (μm)'])
final_output.to_csv(os.path.join(outroot, 
            'summary_{}.csv'.format(time.strftime("%Y-%m-%d %H-%M-%S", time.localtime()))), index=False)

print('Analysis complete!')
logging.info('Analysis complete!')