import numpy as np

# Define grid dimensions
nx = 600
ny = 267

# Define grid spacing
dx = 0.2 / nx
dy = 0.089 / ny

# Initialize volume fractions
material1_fraction = np.zeros((nx, ny))
material2_fraction = np.zeros((nx, ny))

# Define analytical interface equation (circle)
def circle(x, y, xc, yc, r):
    return (x - xc)**2 + (y - yc)**2 - r**2

# Define circle parameters
circle_radius = 0.025
circle_center_x = 0.130
circle_center_y = 0.045

# Loop over grid points
for i in range(nx):
    for j in range(ny):
        x0 = i * dx  # Cell corner coordinates
        y0 = j * dy
        x1 = (i + 1) * dx
        y1 = (j + 1) * dy
        # Check if cell intersects the circle
        if circle(x0, y0, circle_center_x, circle_center_y, circle_radius) <= 0 and \
           circle(x1, y0, circle_center_x, circle_center_y, circle_radius) <= 0 and \
           circle(x0, y1, circle_center_x, circle_center_y, circle_radius) <= 0 and \
           circle(x1, y1, circle_center_x, circle_center_y, circle_radius) <= 0:
            # Cell is fully inside the circle
            material1_fraction[i, j] = 1.0
            material2_fraction[i, j] = 0.0
        elif circle(x0, y0, circle_center_x, circle_center_y, circle_radius) >= 0 and \
           circle(x1, y0, circle_center_x, circle_center_y, circle_radius) >= 0 and \
           circle(x0, y1, circle_center_x, circle_center_y, circle_radius) >= 0 and \
           circle(x1, y1, circle_center_x, circle_center_y, circle_radius) >= 0:
            # Cell is fully outside the circle
            material1_fraction[i, j] = 0.0
            material2_fraction[i, j] = 1.0            
        else:
            # Cell intersects the circle
            # Calculate volume fraction based on intersection area
            intersection_area = 0.0
            nsubcells = 1000
            for sub_i in range(nsubcells):
                for sub_j in range(nsubcells):
                    sub_x = x0 + (sub_i + 0.5) * (dx / nsubcells)
                    sub_y = y0 + (sub_j + 0.5) * (dy / nsubcells)
                    if circle(sub_x, sub_y, circle_center_x, circle_center_y, circle_radius) <= 0:
                        intersection_area += (dx / nsubcells) * (dy / nsubcells)
            cell_area = dx * dy
            volume_fraction = intersection_area / cell_area
            material1_fraction[i, j] = volume_fraction
            material2_fraction[i, j] = 1 - volume_fraction
            
        # Write volume fractions to output file
with open('volume_fractions.dat', 'w') as f:
    for j in range(ny):
        for i in range(nx):
            f.write(f'{i} {j} {material1_fraction[i, j]} {material2_fraction[i, j]}\n')
