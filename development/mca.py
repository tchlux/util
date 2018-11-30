

import numpy as np
from util.plot import Plot


square_test = False
random_points = True

if square_test:
    points = np.array([
        [1, 1],
        [0, 1],
        [1, 0],
        [0, 0],
        [1.1, 1.1],
        [-.1, 1.1],
        [1.1, -.1],
        [-.1, -.1],
        [.5, .5]
    ], dtype=float)
elif not random_points:
    points = np.array([
        [1, 1.5],
        [1.3, 1.3],
        [1.6, .8],
        [.2, -1.2],
        [.3, -1]
    ], dtype=float)
else:
    # Two normals
    size = 200
    points = np.random.normal(size=(size, 2))
    points[:,0] *= 1/3
    # points += np.array([1.,-2.])
    points = np.concatenate((points, np.random.normal(0, 2, size=(size, 2))))
    points[-size:,1]  *= 1/3


dim = 2
vecs = np.random.random(size=(dim, points.shape[1]))
# Normalize the vectors.
for i in range(vecs.shape[0]):
    vecs[i] /= np.linalg.norm(vecs[i], ord=2)

p = Plot()

points = points - np.mean(points, axis=0)

# Given a vector, return all points flipped onto the same side of the vector.
def flipped_points(vec):
    # Get the signs
    signs = np.matmul(points, vec)
    signs[signs == 0] = 1.
    signs /= abs(signs)
    return (points.T * signs).T

# Update each of the vectors to be the average of flipped points.
def update_vecs():
    # vecs[0,:] = np.sum(np.abs(points), axis=0)
    # vecs[0,:] /= sum(vecs[0,:])
    # vecs[1,:] = 0.
    for i in range(vecs.shape[0]):
        # Take the average of the appropriately flipped points
        adjusted_points = flipped_points(vecs[i])
        vecs[i] = np.sum(adjusted_points, axis=0)
        vecs[i] /= np.linalg.norm(vecs[i], ord=2)

# Get the value of each of the vectors.
def value_vecs():
    values = []
    for v in vecs:
        new_pts = flipped_points(v)
        projections = np.matmul(new_pts, v)
        abs_sum = np.sum(abs(projections))
        values.append( abs_sum )
    return np.array(values)

# Reorder the vectors by magnitude of of 1-norm sum (their values).
def reorder_vecs():
    values = value_vecs()
    order = np.argsort(values)[::-1]
    vecs[:,:] = vecs[order][:,:]

# Orthogonalize the vectors.
def orthogonalize_vecs(flip=True):
    q,r = np.linalg.qr(vecs.T)
    vecs[:,:] = q.T[:,:]
    if flip:
        # Flip the vectors to be in an expected direction. (positive)
        signs = np.matmul(vecs, np.ones(vecs.shape[1]))
        signs[signs == 0] = 1.
        signs /= abs(signs)
        vecs[:,:] = (vecs.T * signs).T[:,:]

# Update the plot with the current points and vecs.
def update_p(frame):
    # Add the points to the frame first (so they are on the bottom).
    p.add("points", *(points.T), frame=frame, color=p.color(vecs.shape[0]))
    # Plot the vectors ordered by their value.
    values = value_vecs()
    values /= len(points)
    for i,(v,(x,y)) in enumerate(zip(values, vecs)):
        x, y = x*v, y*v
        p.add("Vec "+str(v+1), [0,x], [0,y], mode="lines",
              marker_size=0, frame=frame, color=p.color(i))
        p.add_annotation(str(round(v,3)), x, y, frame=frame)

previous_vecs = vecs.copy()
steps = 50
update_p("Begin")
for step in range(steps):
    update_vecs()
    update_p(str(step) + " updated")
    reorder_vecs()
    update_p(str(step) + " reordered")
    orthogonalize_vecs()
    update_p(str(step) + " orthogonalized")
    if np.linalg.norm(previous_vecs - vecs) <= 0: break
    else: previous_vecs = vecs.copy()
update_p("End")

p.show(x_range=(-4.5,4.5), y_range=(-2,2), width=1000, height=600,
       show_legend=False, file_name="test.html")

print()
print(vecs)
print(value_vecs())
print()
vecs[:,:] = np.diag(np.ones(2))[:,:]
print(vecs)
print(value_vecs())

print()
from util.stats import ca
a, b = ca(points)
vecs[:,:] = a[:,:]
print(vecs)
print(value_vecs())
print(value_vecs() / sum(value_vecs()))
print(b)
