import os, time, pickle  # System libraries

# Item tree building constants
MAX_ITEMS_PER_LAYER = 10000   # Max number of items kept per layer
SUPPORT_INCREMENT = 0.025     # The amount to incrment support per trim
                              # if there are still too many items
MIN_SUPPORT = 0.1
DEFAULT_ORDER = False

def DEFAULT_OUT(*args, **kwargs):
    kwargs.update({"flush":True}) # Default output forces a "flush"
    print(*args, **kwargs)
 
def _check_consistency(data):
    if (type(data) != list or type(data[0]) != list):
        raise(TypeError("Data given must be a 2D python list."))


# This class will generate an associative item list, given minimum
# support and confidence in a item. Currently only implements support
class AprioriTree:
    # "output" should be a function that takes a string as an argument
    # and has an argument named "end" that declares the line-end behavior
    def __init__(self, support=MIN_SUPPORT, output=DEFAULT_OUT, ordered=DEFAULT_ORDER):
        self.data = []
        self.support = support
        self.ordered = ordered
        if ordered: 
            raise(Exception("No ordered implementation available yet. Consider making all values tuples (<column index>, value)."))
        self.current_items = {}
        self.base_items = {}
        self.base_file_name = "ap_tree_items[%i].pkl"
        self.file_number = 0
        self.print = output # Function used for output
        self.item_sets = [] # Where the final results of mining are stored


    ######## BEGIN ITEM MINING CODE #########

    # Pre:  "self.data" has been set and checked for consistency
    # Post: self.item_sets is established with length 1 items
    def _initialize_items(self):
        self.base_items = {}
        for row in self.data:
            for val in row:
                if (val,) not in self.base_items: 
                    self.base_items[(val,)] = 0
        self.current_items = self.base_items
        

    # Pre:  "new_items" is a dictionary of items
    # Post: "self.current_items" are written to file and the new items
    #       are transitioned into their place
    def _transition_items(self, new_items):
        self.file_number += 1
        with open(self.base_file_name%(self.file_number), "wb") as f:
            pickle.dump(self.current_items, f)
        self.current_items = new_items


    # Pre:  "self.current_items" and "self.base_items" are populated
    #       dictionaries 
    # Post: All the combinatorial possible combinations of the most
    #       recently built item set and the first one are created
    def _generate_items(self):
        new_items = {}
        len_curr_items = len(self.current_items)
        # Cycle the current item sets
        for i,item in enumerate(self.current_items):
            # Update the user
            if not i%(len_curr_items//100 + 1):
                self.print("Generating items (%0.1f%%)"%
                           (100.0*i/len_curr_items),end="\r")
            # Create entries for all possible combinations of the
            # current item set with single items from "self.base_items" 
            for single in self.base_items:
                if (single[0] in item): continue # Skip repeats
                # Sort to ensure there's no redundancy
                item_set = list(item + single); item_set.sort()
                new_items[tuple(item_set)] = 0
        self.print(" \nGenerated %i size %i items."%
                   (len(new_items),self.file_number+2))
        self._transition_items(new_items)
        return len(new_items) > 0


    # Pre:  "self.item_sets" has a new dictionary on the end where
    #       values of keys are defined by the amount of support in data
    # Post: All items that do not have enough support are removed
    def _trim_items(self, support=None):
        min_support = self.support if (support == None) else support
        to_remove = []
        # Create a generator so we can pop while iterating
        curr_items_gen = (item for item in list(self.current_items.keys()))
        len_curr_items = len(self.current_items)
        for i,item in enumerate(curr_items_gen):
            if not i%(len_curr_items//100 + 1):
                self.print("Removing non-supported items (%0.1f%%)"%
                           (100.0*i/len_curr_items),end="\r")
            # Pop out all items that do not meet minimum support
            if self.current_items[item] < min_support:
                self.current_items.pop(item)
        # If there are still too many items, trim off more
        if len(self.current_items) > MAX_ITEMS_PER_LAYER:
            self._trim_items(min_support + SUPPORT_INCREMENT)
        else: # Otherwise, done trimming
            self.print("%i items at %0.1f%% support"%
                       (len(self.current_items), 100.0*min_support))


    # Pre:  "self.current_items" and "self.data" have been
    #       appropriately populated
    # Post: Support for all items in "self.current_items" is found
    #       through "self.data"
    def _find_item_support(self):
        support_per_instance = 1.0 / len(self.data)
        start = time.time()
        etr = lambda per: (time.time()-start) / per * (1-per)
        for i,row in enumerate(self.data):
            if not i % 10:
                self.print("(%0.1f%%) with [%is] left"%
                           (100.0*i/len(self.data), 
                            etr(float(i+1)/len(self.data))), end="\r")

            for item in self.current_items:
                if all(value in row for value in item):
                    self.current_items[item] += support_per_instance


    # Pre:  "data" is a python 2D list, "longest_item" integer > 0
    # Post: All items are mined from data using the apriori tree
    #       pruning method
    def mine_items(self, data, longest_item=None):
        # Reset the class
        self.file_number = 0; self.base_items = {};
        self.item_sets = [];  self.current_items = {}
        # Check and store data
        _check_consistency(data)
        self.data = data
        # Automatically set longest item if it is not given
        if longest_item == None:
            longest_item = max([len(row) for row in self.data])
        # Create "base_items" dictionary
        self._initialize_items()
        # Loop through and mine all potential item sets
        for size in range(longest_item):
            # Provide update
            self.print("Potential items of size %i (%0.1f%%): %i"%
                       ( size + 1, 100.0 * (size+0.5) / longest_item, 
                         len(self.current_items) ))
            # Find support for all of the items in self.current_items
            self._find_item_support()
            # Remove items without desired support
            self._trim_items()
            # Break if there are no items to build from
            if (len(self.current_items) == 0): break 
            elif (size+1 < longest_item): # generate next group
                # Break if no new items could be generated (this check
                # actually performs the item generation process)
                if not self._generate_items(): break 
            else: # If there are current items and there are no longer
                  # item sets, then save current items to file
                self._transition_items({})                

        # Load all of the item sets into self.item_sets from temporary files
        self.print("\nLoading all item dictionaries from files..")
        for i in range(self.file_number):
            self.print("Loading file %s - (%0.1f%%)"%(
                self.base_file_name%(i+1),100.0*i/self.file_number),end="\r")
            with open(self.base_file_name%(i+1), "rb") as f:
                self.item_sets.append(pickle.load(f))
            os.remove(self.base_file_name%(i+1))
        self.print("Done loading files.\n")
        # Return the list of dictionaries of ( item set : support )
        return self.item_sets


# Test code for this file
if __name__ == "__main__":
    a = AprioriTree()#output=lambda *args, **kwargs:None)
    d = [[0,1,2,6],
         [2,1,5,8],
         [8,4,3,2],
         [7,2,5,3],
         [1,2,0,5],
         [4,1,7,9],
         [3,6,8,1],
         [1,8,9,0],
         [0,2,1,4],
         [1,9,4,7]]
    item_sets = a.mine_items(d)
    for item_set in item_sets: 
        items = list(item_set.items())
        items.sort(key=lambda i: -i[1])
        for i in items: print("%s: %0.1f%%"%(i[0], 100.0*i[1]))
        print()
