from util.data import *

# Define a class for reading a video file without actually storing it in memory.
class IndexedVideo:
    import cv2
    import numpy as np
    class FailedSet(Exception): pass
    class FailedRead(Exception): pass
    class BadIndex(Exception): pass
    class IndexTooLarge(Exception): pass
    # Initialize by getting the information about this video from file.
    def __init__(self, video_file_name, flatten=False, use_float=False):
        # Store the 'flatten' attirbute for flattening images.
        self.flatten = flatten
        self.float = use_float
        # Store the video using an OpenCV video capture object (with its attributes)
        self.vid = self.cv2.VideoCapture(video_file_name)
        self.frames = int(self.vid.get(self.cv2.CAP_PROP_FRAME_COUNT))
        self.width  = int(self.vid.get(self.cv2.CAP_PROP_FRAME_WIDTH))
        self.height = int(self.vid.get(self.cv2.CAP_PROP_FRAME_HEIGHT))
        self.fps    = int(self.vid.get(self.cv2.CAP_PROP_FPS))
        # Find the last available index in the video using the "set"
        # method. This must be done because of a bug in OpenCV.
        self.length = self.frames-1
        for i in range(-1, -self.frames, -1):
            # If the index is successfully accessed, break.
            try:   
                img = self[i]
                self.length = self.frames + i
                break
            # Otherwise, ignore the failed read operation.
            except self.FailedRead: pass
        # Use the first image to get the full shape of this video.
        self.shape = (self.length,) + self[0].shape
        self.file_name = video_file_name

    # Define the "length" of this video.
    def __len__(self): return self.length

    # Define the index operator for this video.
    def __getitem__(self, index):
        # Cover two possible (common) usage errors.
        if type(index) != int: raise(self.BadIndex("Only integer indices allowed."))
        if index > self.length: raise(self.IndexTooLarge(f"Index {index} greater than length of self, {self.length}."))
        # Convert negative indices into positive indices.
        if (index < 0):
            while (index < 0): index += self.length
        # Move the capture to the correct frame of the video.
        success = self.vid.set(self.cv2.CAP_PROP_POS_FRAMES, index)
        if not success: raise(self.FailedSet(f"Failed to set to video frame {index}."))
        # Get the image at that frame.
        success, image = self.vid.read()
        if not success: raise(self.FailedRead(f"Failed to read video frame {index}."))
        # Flatten / float the image if that setting is enabled.
        if self.flatten: image = image.flatten()
        if self.float:   image = self.np.asarray(image, dtype=float)
        # Return the image at that frame.
        return image
