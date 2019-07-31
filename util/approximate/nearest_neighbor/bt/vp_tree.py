# This code seems to be complete.
# There might be a bug, but there probably isn't.
# The VP tree method tends to have nearly no pruning in high dimension (random) data.
# This makes sense, because all points are nearly equally spaced in high dimension.


import numpy as np

class VPTree:
    # A node in this tree is either a leaf or a branch. Branches are
    # defined by a center and a radius.
    class VPNode:
        leaf = True          # True -> This is a leaf node.
        parent = None        # Node preceding this node in tree.
        points = None        # List of all points.
        sq_lens = None       # List of squared lengths of all points.
        indices = None       # List of original indices of all points.
        center = None        # Point from original array defining center.
        index = None         # Index of "center" point in original array.
        radius = None        # Median distance for non-leaf nodes.
        shell = float('inf') # Furthest distance from center.
        contains = 0         # Number points contained (including center).

        # Handle the construction of this node.
        def __init__(self, points, indices, sq_lens, leaf_size, parent=None):
            self.parent = parent
            # Compute the number of contained points.
            self.contains = len(points)
            # Compute a center (a point on the boundary of the space), 
            # swap it to the front of the set of points.
            min_corner = np.min(points, axis=0)
            sq_dist_to_min = sq_lens + np.sum(min_corner**2) - 2 * np.dot(points,min_corner)
            to_swap = np.argmin(sq_dist_to_min)
            temp = points[0].copy()
            points[0,:] = points[to_swap,:]
            points[to_swap,:] = temp[:]
            # Store the information about the center point.
            self.index = 0
            self.center = points[0]
            self.index = indices[0]
            center_sq = sq_lens[0]
            # Store contained data.
            self.points = points[1:]
            self.indices = indices[1:]
            self.sq_lens = sq_lens[1:]
            # Compute the distances to all contained points (from center).
            if len(self.points > 0):
                distances_sq = (center_sq + self.sq_lens - 2 * np.dot(self.points, self.center))
            # Partition the set of points further if we are not done.
            if (len(self.points) > leaf_size):
                max_val = np.max(distances_sq)
                self.shell = np.sqrt(max_val) if max_val > 0 else 0
                self.leaf = False
                radius_sq = np.median(distances_sq)
                inside = distances_sq <= radius_sq
                outside = np.logical_not(inside)
                self.radius = np.sqrt(radius_sq) if radius_sq > 0 else 0
                self.inside = type(self)(self.points[inside], self.indices[inside], 
                                         self.sq_lens[inside], leaf_size, self)
                self.outside = type(self)(self.points[outside], self.indices[outside],
                                          self.sq_lens[outside], leaf_size, self)
            # else:
            #     print("shell:", self.shell)

        # Method for displaying information about this Node.
        def __str__(self):
            s =  f"Node:\n"
            s += f"  leaf     - {self.leaf}\n"
            s += f"  center   - {self.center}\n"
            s += f"  contains - {self.contains}\n"
            s += f"  radius   - {self.radius}"
            return s
            
    # Build this tree.
    def __init__(self, points, leaf_size=10):
        self.points = points
        self.leaf_size = leaf_size
        indices = np.arange(len(points))
        sq_lens = np.sum(points**2, axis=1)
        self.root = self.VPNode(points, indices, sq_lens, leaf_size)


    # # Print out the points (with edges showing ownership on a plot).
    # def plot(self):
    #     from util.plot import Plot
    #     p = Plot()
    #     nodes = {}
    #     queue = [(0, self.root)]
    #     while len(queue) > 0:
    #         level, node = queue.pop(0)
    #         nodes[level] = nodes.get(level,[]) + [node]
    #         if not node.leaf:
    #             queue += [(level+1,node.inside), (level+1,node.outside)]
    #     # Plot all levels of the tree.
    #     for level in sorted(nodes):
    #         for n in nodes[level]:
    #             c = p.color_num
    #             name = f"{n.index} (level {level})"
    #             p.add(name, [
    #             if n.leaf:


    # Print out the hierarchical structure of this tree.
    def graph(self):
        from util.plot import Plot
        g = Plot()
        y = 0
        queue = [(0,self.root)]
        # Add the root node.
        g.add_node(f"{self.root.index} {{{self.root.contains}}}", 0, 0,
                   label=True, label_y_offset=0)
        while (len(queue) > 0):
            y -= 1
            new_queue = []
            for (x, node) in queue:
                if not node.leaf:
                    left_x = x - (2**y) + 1/(2*(y-1))
                    right_x = x + (2**y) - 1/(2*(y-1))
                    parent_name = f"{node.index} {{{node.contains}}}"
                    left_name  = f"{node.inside.index} {{{node.inside.contains}}}"
                    right_name = f"{node.outside.index} {{{node.outside.contains}}}"
                    # Make the names inlude contents for leaves.
                    if (node.inside.leaf):  left_name  += f" {node.inside.indices}"
                    if (node.outside.leaf): right_name += f" {node.outside.indices}"
                    # Add the nodes to the graph.
                    g.add_node(left_name, left_x, y, label=True, label_y_offset=0)
                    g.add_node(right_name, right_x, y, label=True, label_y_offset=0)
                    g.add_edge([left_name, parent_name])
                    g.add_edge([right_name, parent_name])
                    new_queue += [(left_x, node.inside), (right_x, node.outside)]
            queue = new_queue
        g.graph(file_name="tree.html")
        

    # Compute the nearest "k" neighbors to each query point.
    def query(self, pts, k=1, return_distance=True):
        d_calcs = 0
        nodes = 1
        # Make sure "pts" is 2D
        try:    len(pts[0])
        except: pts = [pts]
        k = min(k, len(self.points))
        all_distances = []
        all_indices = []
        # Cycle through all provided points.
        for pt in pts:
            # Perform VP tree search.
            # Initialize storage for nearest points, num found, explore stack.
            indices = np.ones(k+2 + self.leaf_size, dtype=int) + len(self.points)+1
            distances = np.ones(k+2 + self.leaf_size, dtype=float) * float('inf')
            pt_sq = np.sum(pt**2)
            count = 0
            stack = [self.root]
            # Loop until the stack is empty (all appropriate leaf nodes explored).
            while (len(stack) > 0):
                node = stack.pop() # Pop out last element (depth first).
                # Measure distance to vantage point.
                d = np.linalg.norm(node.center - pt)
                d_calcs += 1
                indices[count] = node.index
                distances[count] = d
                count += 1
                # If the furthest point being kept is closer than the
                # nearest possible point in the node, then skip.
                if (distances[k-1] <= (d - node.shell)): continue
                # If this is a leaf node, measure distance to all points.
                elif node.leaf and (node.contains > 1):
                    end = count + len(node.indices)
                    # Compute distances to all points.
                    indices[count:end] = node.indices
                    d = pt_sq + node.sq_lens - 2*np.dot(node.points, pt)
                    gt_zero = d > 0
                    d[gt_zero] = np.sqrt(d[gt_zero])
                    distances[count:end] = d
                    d_calcs += len(d)
                    # Update the "count" of stored indices and distances.
                    count = end
                # Otherwise this must be a branch node, determine exploration.
                elif (not node.leaf):
                    max_dist = max(distances[:k-1])
                    is_inside = (d <= node.radius)
                    # If there could be a closer point inside or there are
                    # not enough points outside to fill "k" requirement.
                    if (is_inside or ((d - node.radius) < max_dist) or
                        (node.outside.contains < (k - count))):
                        # print(f"Adding inside node {node.inside.index}..   ",
                        #       f"{d:.2f} {max_dist:.2f} {node.radius:.2f}    ",
                        #       (is_inside, ((d - node.radius) < max_dist),
                        #        (node.outside.contains < (k - count))))
                        stack.append( node.inside )
                        nodes += 1
                    # If there could be a closer point outside or there 
                    # aren't enough points inside to fill "k" requirement.
                    if (not is_inside or ((node.radius - d) < max_dist) or
                        (node.inside.contains < (k - count))):
                        # print(f"Adding outside node {node.outside.index}..   ",
                        #       f"{d:.2f} {max_dist:.2f} {node.radius:.2f}    ",
                        #       (not is_inside, ((node.radius - d) < max_dist),
                        #        (node.inside.contains < (k - count))))
                        stack.append( node.outside )
                        nodes += 1
                    # If the point is inside, then make sure that the
                    # "inside" node is the first to be explored (by
                    # being last in the stack).
                    if (is_inside and (stack[-1] == node.outside)):
                        stack[-2], stack[-1] = stack[-1], stack[-2]
                # If there are too many, move all lower distances
                # into the first "k" spots, update count to reflect.
                if (count > k):
                    to_keep = np.argpartition(distances[:count], k)[:k]
                    indices[:k] = indices[to_keep]
                    distances[:k] = distances[to_keep]
                    count = k

            assert(count >= k)
            sorted_idxs = np.argsort(distances[:k])
            all_indices.append( indices[sorted_idxs] )
            all_distances.append( distances[sorted_idxs] )
        # Convert all results into a NumPy array.
        indices = np.array(all_indices)
        print("Nodes explored:      ", nodes)
        print("Distances calculated:", d_calcs)
        # Return the appropriate info depending on usage.
        if return_distance:
            distances = np.array(all_distances)
            return (distances, indices)
        else:
            return indices
            

