from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF, QObject
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF, QObject
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import math

# Some global color constants that might be useful
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)

# Global variable that controls the speed of the recursion automation, in seconds
#
PAUSE = 0.25

# nlogn time and space. logn for recursion and n for merge calls
def convex_hull_solver(points):
    # Base case- if length of array is 4 we return in
    if len(points) < 4:
        return points

    # get new left and right array by calling function to split in half
    left, right = get_halves(points)
    # recursive calls
    left = convex_hull_solver(left)
    right = convex_hull_solver(right)

    # if right or left array has 3 points, we need to order clockwise based on slope
    if len(right) == 3 and compute_slope(right[0], right[2]) > compute_slope(right[0], right[1]):
        # time: O(1) space: O(1)
        temp = right[1]
        right[1] = right[2]
        right[2] = temp
    if len(left) == 3 and compute_slope(left[0], left[2]) > compute_slope(left[0], left[1]):
        # time: O(1) space: O(1) because hull size is max of 3
        temp = left[1]
        left[1] = left[2]
        left[2] = temp
    hull = merge_hull(left, right)

    # return final array
    return hull

# time: O(1) # space O(n) because we have two n/2 arrays
def get_halves(n):
    # takes in array n as parameter and returns two arrays splitting n from the middle
    half = len(n) // 2
    return n[:half], n[half:]


# time: O(n) because we call top_tan_indexes and bottom_tan_indexes which are both O(n)
# space: O(n) because we are creating 2 arrays of size n/2 worst case
def merge_hull(left, right):
    # calls functions to find top and bottom tangents
    top_tan_indexes = find_top_tan(left, right)
    bottom_tan_indexes = find_bottom_tan(left, right)
    # merge the lists of top and bottom tangents together
    new_hull_indexes = top_tan_indexes + bottom_tan_indexes
    sorted_hull = merge_clockwise(left, right, new_hull_indexes)

    return sorted_hull

# time: would be O(2n) but simplifies to O(n) because worst case is you keep all points
# space: O(n)- create new list of size n worst case
def merge_clockwise(left, right, tan_indexes):
    # this function merges the two hulls and sorts them clockwise
    clockwise_list = []

    # start with left most point in array and append up until upper left tangent
    for point in left[:tan_indexes[0] + 1]:
        clockwise_list.append(point)

    # if upper right tangent index is greater than bottom right tangent index loop to end of right list
    # and then all the way back to the bottom right tangent point (including bottom right point)
    if tan_indexes[1] > tan_indexes[2]:
        for point in right[tan_indexes[1]: len(right)]:
            clockwise_list.append(point)
        for point in right[:tan_indexes[2] + 1]:
            clockwise_list.append(point)
    else:
        for point in right[tan_indexes[1]: tan_indexes[2] + 1]:
            clockwise_list.append(point)

    # if the bottom right tangent index is not 0, append all points until end of right list
    if tan_indexes[3] != 0:
        for point in left[tan_indexes[3]:len(left)]:
            clockwise_list.append(point)

    return clockwise_list

# time: would be O(2n) but simplifies to O(n) because worst case is you keep all points
# space: O(1)
def find_top_tan(left, right):
    new_hull_indexes = []

    # O(n) being added to total complexity but dropped because we are adding
    curr_left = find_rightmost_ind(left)
    curr_right = 0
    is_upper_tangent = False
    curr_slope = compute_slope(left[curr_left], right[curr_right])
    while not is_upper_tangent:
        left_updated = False
        right_updated = False
        is_tangent_left = False
        # start outer loop that runs until no changes in the left or right hulls have been made
        while not is_tangent_left:
            # account for the chance that we need to move clockwise past index 0 or left array
            if curr_left == 0:
                next_left = len(left) - 1
            else:
                next_left = curr_left - 1
            new_slope = compute_slope(left[next_left], right[curr_right])
            if new_slope < curr_slope:
                curr_slope = new_slope
                curr_left = next_left
                left_updated = True
            else:
                is_tangent_left = True
        is_tangent_right = False
        curr_slope = compute_slope(left[curr_left], right[curr_right])
        # start outer loop that runs until no changes in the left or right hulls have been made
        while not is_tangent_right:
            next_right = curr_right + 1
            if curr_right == len(right) - 1:
                next_right = 0
            new_slope = compute_slope(left[curr_left], right[next_right])
            if new_slope > curr_slope:
                curr_slope = new_slope
                curr_right = next_right
                right_updated = True
            else:
                is_tangent_right = True
        # If we updated a left or right point we need to check again, otherwise we are done
        if left_updated or right_updated:
            is_upper_tangent = False
        else:
            is_upper_tangent = True

    # append top tangent points in clockwise order
    new_hull_indexes.append(curr_left)
    new_hull_indexes.append(curr_right)
    return new_hull_indexes

# time: would be O(2n) but simplifies to O(n) because worst case is you keep all points
# space: O(1)
def find_bottom_tan(left, right):
    new_hull_indexes = []
    curr_left = find_rightmost_ind(left)
    curr_right = 0
    is_upper_tangent = False

    curr_slope = compute_slope(left[curr_left], right[curr_right])

    while not is_upper_tangent:
        left_updated = False
        right_updated = False
        is_tangent_left = False

        while not is_tangent_left:

            # if current right index is last point of array our next point will be the beginning of array
            if curr_left == len(left) - 1:
                next_left = 0
            else:
                next_left = curr_left + 1

            new_slope = compute_slope(left[next_left], right[curr_right])

            if new_slope > curr_slope:
                curr_slope = new_slope
                curr_left = next_left
                left_updated = True
            else:
                is_tangent_left = True

        is_tangent_right = False
        curr_slope = compute_slope(left[curr_left], right[curr_right])

        while not is_tangent_right:

            # if the current right index is zero our next point will be starting at end or array
            if curr_right == 0:
                next_right = len(right) - 1
            else:
                next_right = curr_right - 1

            new_slope = compute_slope(left[curr_left], right[next_right])

            if new_slope < curr_slope:
                curr_slope = new_slope
                curr_right = next_right
                right_updated = True
            else:
                is_tangent_right = True

        # If we updated a left or right point we need to check again, otherwise we are done
        if left_updated or right_updated:
            is_upper_tangent = False
        else:
            is_upper_tangent = True

    # append top tangent points in clockwise order
    new_hull_indexes.append(curr_right)
    new_hull_indexes.append(curr_left)
    return new_hull_indexes

# time: O(n) space:
def find_rightmost_ind(arr):
    right_ind = 0
    for i in range(len(arr)):
        if arr[i].x() > arr[right_ind].x():
            right_ind = i
    return right_ind

# time: O(1) space:
def compute_slope(p1, p2):
    return (p2.y() - p1.y()) / (p2.x() - p1.x())


#
# This is the class you have to complete.
#
class ConvexHullSolver(QObject):

    # Class constructor
    def __init__(self):
        super().__init__()
        self.pause = False

    # Some helper methods that make calls to the GUI, allowing us to send updates
    # to be displayed.

    def showTangent(self, line, color):
        self.view.addLines(line, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseTangent(self, line):
        self.view.clearLines(line)

    def blinkTangent(self, line, color):
        self.showTangent(line, color)
        self.eraseTangent(line)

    def showHull(self, polygon, color):
        self.view.addLines(polygon, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseHull(self, polygon):
        self.view.clearLines(polygon)

    def showText(self, text):
        self.view.displayStatusText(text)

    # This is the method that gets called by the GUI and actually executes
    # the finding of the hull
    def compute_hull(self, points, pause, view):
        self.pause = pause
        self.view = view
        assert (type(points) == list and type(points[0]) == QPointF)

        t1 = time.time()

        # time complexity O(nlogn)
        # space complexity 0(n)
        ordered_points = sorted(points, key=lambda QPointF: QPointF.x())

        t2 = time.time()

        print('Time Elapsed (Sorting): {:3.3f} sec'.format(t2 - t1))

        t3 = time.time()
        hull = convex_hull_solver(ordered_points)

        t4 = time.time()
        polygon = [QLineF(hull[i], hull[(i + 1) % len(hull)]) for i in range(len(hull))]
        # when passing lines to the display, pass a list of QLineF objects.  Each QLineF
        # object can be created with two QPointF objects corresponding to the endpoints
        self.showHull(polygon, RED)
        self.showText('Time Elapsed (Convex Hull): {:3.3f} sec'.format(t4 - t3))



# Empirical Data: Convex Hull Elapsed Times                 Finding k value                 k = mean measured time/nlogn
# When n = 10 [0.000, 0.000, 0.000, 0.000, 0.000]
# When n = 100 [0.001, 0.000, 0.000, 0.000, 0.000]          100*log(100) = 200              k = 0.0002/200 = 0.000001
# When n = 1000 [0.008, 0.006, 0.006, 0.006, 0.007]         1000*log(1000) = 3000           k = 0.0066/3000 = 0.0000022
# When n = 10000 [ 0.058, 0.059, 0.058, 0.059, 0.058]       10000*log(10000) = 40000        k = 0.0584/40000 = 0.00000146
# When n = 100000 [0.000, 0.000, 0.000, 0.000, 0.000]       100000*log(100000) = 500000     k = 0.5688/500000 = 0.00000114
# When n = 500000 [0.000, 0.000, 0.000, 0.000, 0.000]       500000*log(500000) = 2849485    k = 3.531/2849485 = 0.00000124
# When n = 1000000 [0.000, 0.000, 0.000, 0.000, 0.000]      1000000*log(1000000) = 6000000  k = 6.4116/6000000 = 0.00000107



# Empirical Data for Sorting Elapsed Times
# When n = 10 [0.000, 0.000, 0.000, 0.000, 0.000]
# When n = 100 [0.000, 0.000, 0.000, 0.000, 0.000]
# When n = 1000 [0.000, 0.000, 0.000, 0.000, 0.000]
# When n = 10000 [0.000, 0.000, 0.000, 0.000, 0.000]
# When n = 100000 [0.000, 0.000, 0.000, 0.000, 0.000]
# When n = 500000 [0.000, 0.000, 0.000, 0.000, 0.000]
# When n = 1000000 [0.000, 0.000, 0.000, 0.000, 0.000]



