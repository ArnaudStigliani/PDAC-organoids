import cv2
import matplotlib.pyplot as plt
import seaborn as sns

def nice_annotation(img, bbox, color, label):
    '''
    Adapted from https://stackoverflow.com/questions/56108183/python-opencv-cv2-drawing-rectangle-with-text
    '''
    x1, y1, x2, y2  = bbox
    img = cv2.rectangle(img, (x1, y1), (x2, y2), color, 2)
    (w, h), _ = cv2.getTextSize(
        label, cv2.FONT_HERSHEY_SIMPLEX, 0.6, 1)

    # Prints the text.
    img = cv2.rectangle(img, (x1, y1 - 20), (x1 + w, y1), color, -1)
    img = cv2.putText(img, label, (x1, y1 - 5),
                      cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 255), 2) # white text color

def my_histogram(df, log_scale = True):
    f = plt.figure(figsize=(6, 3))
    gs = f.add_gridspec(1, 1)
    with sns.axes_style("ticks"):
        ax = f.add_subplot(gs[0, 0])
        ax = sns.histplot(df,
                    x="Ellipse_Area",
                bins=30,
                ax = ax,
                kde=True,
                log_scale=log_scale)
        ax.set(xlabel='Ellipse_Area (μm^2)', ylabel='Count')
        ax.axvline(x=df.Ellipse_Area.mean(),
                color='red', linestyle = '--')
        ax.axvline(x=df.Ellipse_Area.median(),
                color='red', linestyle = '-')
    return f
