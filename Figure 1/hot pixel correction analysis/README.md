Code to correct new dead pixels using the PixelFly Camera. While the camera still works well, the model is too old to update to use the companies provided software for pixel correction.

This code uses coordinates of dead pixels (chosen specifically for the camera based on thresholding), then uses the average of the neighboring pixels to fill the dead pixel with an estimated value.
Improves the performance of movement correction software.
