from util.approximate import Approximator, test_plot
import numpy as np

# Inverse distance weighted values
class Shepard(Approximator):
    def fit(self, x, y):
        self.x = x
        self.y = y

    def predict(self, x, exp=5):
        response = []
        for pt in x:
            distances = np.sqrt(np.sum((self.x - pt)**2, axis=1))
            # Handle the special case of predicting a previously seen point
            if min(distances) == 0:
                response.append( self.y[np.argmin(distances)] )
                continue
            distances /= max(distances)
            weights = distances**(-1)
            weights = weights**exp
            weights /= sum(weights)
            response.append( sum(weights*self.y) )
        return np.array(response)

# Linear interpolation between nearest containing points. If the
# chosen dimension 
class NearbySimplex(Approximator):
    def fit(self, x, y):
        self.x = x
        self.y = y

    def predict(self, x, num_points=None):
        # Default number of points is a full simplex
        if (type(num_points) == type(None)):
            num_points = x.shape[1] + 1
        # Initial storage for response values
        response = []
        # Cycle making predictions for each point independently
        for pt in x:
            # Get the distances between all known points and the current point
            distances = np.sqrt(np.sum((self.x - pt)**2, axis=1))
            # Handle the special case of predicting a previously seen point
            if min(distances) == 0:
                response.append( self.y[np.argmin(distances)] )
                continue
            # Identify the simplex (of desired dimension) that roughly
            # contains this interpolation point.
            simp = [np.argmin(distances)]
            for _ in range(1,num_points):
                # Look for a point on the other side of the
                # interpolation point from all current points
                found = False
                others = float('inf') * distances
                for i in range(len(self.x)):
                    if (i in simp): continue
                    # If this point is on the other side, use it
                    if all(np.dot(self.x[s] - pt, self.x[i] - pt) < 0
                           for s in simp):
                        others[i] = distances[i]
                        found = True
                # If no point exists on other side, exit
                if not found: break
                # Otherwise, add the closest point to the simplex
                simp += [np.argmin(others)]
            # Use the voronoi mesh to do interpolation with this
            # (potentially underdefined) simplex
            from util.approximate import VoronoiMesh
            sub_model = VoronoiMesh()
            sub_model.fit(self.x[simp], self.y[simp])
            response.append( sub_model(pt) )
        # Return all response values as an array
        return np.array(response)

# The minimimum feasible Lipschitz approximator
class LipschitzMedian(Approximator):
    def fit(self, x, y, verbose=True):
        self.x = x
        self.y = y
        self.l = 0.0
        for p1,fp1 in zip(x,y):
            for p2,fp2 in zip(x,y):
                dist = np.sqrt(sum((p1 - p2)**2))
                if (fp1 == fp2): continue
                pair_l = abs(fp2 - fp1) / dist
                self.l = max(self.l, pair_l)
        print(" Lipschitz Constant:", self.l)

    def predict(self, x):
        response = []
        for pt in x:
            distances = np.sqrt(np.sum((self.x - pt)**2, axis=1))
            # Handle the special case of predicting a previously seen point
            if min(distances) == 0:
                response.append( self.y[np.argmin(distances)] )
                continue
            # Calculate the lipschitz feasible region
            lower = max(self.y - distances * self.l)
            upper = min(self.y + distances * self.l)
            response.append( (upper + lower) / 2 )
        return np.array(response)


# The minimimum feasible Lipschitz approximator, uses the value
class MinimumLipschitz(Approximator):
    def fit(self, x, y):
        self.x = x
        self.y = y

    def predict(self, x):
        response = []
        for pt in x:
            distances = np.sqrt(np.sum((self.x - pt)**2, axis=1))
            # Handle the special case of predicting a previously seen point
            if min(distances) == 0:
                response.append( self.y[np.argmin(distances)] )
                continue
            # Find the value such that 
            num_steps = 0
            min_l = 1.0
            change = 1.0
            while (change > 2**(-8)):
                num_steps += 1
                # Calculate the lipschitz feasible region
                lower = max(self.y - distances * min_l)
                upper = min(self.y + distances * min_l)
                # Adjust the estimated lipschitz constant
                if (upper > lower):
                    change /= 2
                    min_l -= change
                else:
                    min_l += change

            # Need a way to remove 'peaks', add second derivative
            # continuity around known interpolation points.

            # Use overtaking to smooth the function
            lower = self.y - distances * min_l
            upper = self.y + distances * min_l
            max_influence = np.argmax(lower)
            min_influence = np.argmin(upper)
            # Get the next closest max
            lower[max_influence] = 0
            next_max = np.argmax(lower)
            # Get the next closest min
            upper[min_influence] = float('inf')
            next_min = np.argmin(upper)
            # Recompute the lower and upper
            lower = self.y - distances * min_l
            upper = self.y + distances * min_l
            # Check for overtaking
            od = min(abs(lower[max_influence] - lower[next_max]),
                     abs(upper[min_influence] - upper[next_min]))
            min_l += max(0,(od+.2) if (od < -.1) else (-od))
            min_l += max(0,(od) if (od < .1) else (.2-od))
            

            # # When the influencers are on the same side, use the closer.
            # same_side = 0 < np.dot(self.x[max_influence] - pt,
            #                        self.x[min_influence] - pt)
            # if same_side:
            #     if distances[max_influence] < distances[min_influence]:
            #         response.append( self.y[max_influence] )
            #     else:
            #         response.append( self.y[min_influence] )
            #     continue

            # print(max_influence, min_influence)
            # print(distances[max_influence])
            # print(distances[min_influence])
            # print(min_l)
            # min_l /= min(distances)

            # If the min and max bounds are provided from the same
            # side, tend towards the closer one with increased
            # distance away from them both.
            
            # Consider what points are about to overtake in terms of
            # the defining bound (by closeness to max / min).
            # If an overtake is about to occur,
            #  raise the lipschitz constant estimate
            # If the distance to the closest bound is dominating ,
            #  raise the lipschitz constant estimate

            response.append( (max(self.y - distances * min_l) + 
                              min(self.y + distances * min_l)) / 2 )
        return np.array(response)


if __name__ == "__main__":
    from util.plotly import Plot
    f = lambda od: (max(0,(od+.2) if (od < -.1) else (-od)) +
                    max(0,(od) if (od < .1) else (.2-od)))
    p = Plot()
    p.add_func("", f, [-1, 1])
    _ = p.show()
    exit()

    # from util.approximate import Delaunay as model
    # model = LipschitzMedian
    model = MinimumLipschitz

    p,_,_ = test_plot(model, N=3, D=1, random=False)


    # Generate a test plot showing this algorithm
    # p = test_plot(model, N=90, D=2, plot_points=2000, random=False,
    #               low=-.5, upp=1.5)

    p.show(file_name="test_plot.html")

    # - define 'neighbors'
    # - calculate distance to neighbors
    # - modify distances to have desirable properties
    # - produce output


    # - only consider interpolation, no extrapolation
    # - when moving through the space, the change should be nearly linear
    #   with respect to the points which have the largest distance change
    # - the point that is moved away from the most and is closest should
    #   be combined with the point that is moved towards most and is the
    #   closest
    # - consider closeness and alignment

    # - find nearest neighbor, find nearest point on other side
    # - repeat until simplex is formed or done
    # - solve for linear combination of vertices that is closest to point
