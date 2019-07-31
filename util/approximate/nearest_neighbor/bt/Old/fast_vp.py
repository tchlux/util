
# This function rearranges a set of provided points into a vantage
# point tree and returns the split radius and shell size of all nodes
# in the tree on output.
# 
# The tree is stored implicitly with each 
# 
# Arguments:
#   points    -- matrix of floats (D x N)
#   leaf_size -- positive integer determining points per leaf node.
#   random    -- method for determining center, True -> random selection,
#                False -> point furthest from center of mass.
# 
# Returns:
#   points -- matrix of floats (D x N)
#   splits -- vector of floats for non-leaf node split radii (N)
#   shells -- vector of floats for node outer shell distance (N)
def build_vp_tree(points, leaf_size=10, random=True):
    
    if not random:
        # Compute the center of mass and make it the first point.
        pass
    return points, splits, shells


 
   
