import math
import copy


class State:

    def __init__(self, parent_state, city_ind, cities):
        # properties for states other than initial
        if parent_state != None:
            self.parent_state = parent_state
            # copy the parent states matrix
            self.matrix = copy.deepcopy(parent_state.matrix)

            # increment depth
            self.depth = self.parent_state.depth + 1

            # initialize lower bound to zero, will be updated in reduce_matrix()
            self.state_lower_bound = 0

            # set states parent lower bound to the parent's state lower bound
            self.parent_lower_bound = self.parent_state.state_lower_bound

            # copy rows and columns visited, and also the path from the parent
            self.set_rows = copy.deepcopy(parent_state.set_rows)
            self.set_cols = copy.deepcopy(parent_state.set_cols)
            self.path = copy.deepcopy(parent_state.path)
            self.path_set = copy.deepcopy(parent_state.path_set)

            # set to and from city indexes
            self.from_city_ind = parent_state.to_city_ind
            self.to_city_ind = city_ind

            self.reduce_matrix(cities)

    # Time: O(n^2)
    # Space: O(n) copy parent's matrix
    def reduce_matrix(self, cities):
        # function to reduce states other than initial state

        # add cost of from city to city we want to go to
        cost_of_path = self.matrix[self.from_city_ind][self.to_city_ind]

        # block off column of city you are potentially going to
        # block off row of city you are at
        # set inverse to infinity
        # set rows and columns to infinity

        # infinite out row from_city_index
        for i in range(len(self.matrix)):
            self.matrix[self.from_city_ind][i] = math.inf

        # infinite out column to_city_index
        for i in range(len(self.matrix)):
            self.matrix[i][self.to_city_ind] = math.inf

        self.matrix[self.to_city_ind][self.from_city_ind] = math.inf

        # add from city to set rows and to city to set columns so we don't reduce on these values
        self.set_rows.add(self.from_city_ind)
        self.set_cols.add(self.to_city_ind)

        # get cost of reducing rows and columns
        cost_reduce = self.reduce_rows()
        cost_reduce += self.reduce_cols(self.to_city_ind)

        # add cost of row and col to state's lower bound
        self.state_lower_bound = cost_reduce + cost_of_path + self.parent_lower_bound

        # add current city to path set and array
        self.path.append(cities[self.to_city_ind])
        self.path_set.add(self.to_city_ind)

    def reduce_first_state(self):
        cost_row = self.reduce_rows()
        cost_col = self.reduce_cols(self.to_city_ind)
        self.state_lower_bound = cost_row + cost_col

    def set_first_state(self, matrix, cities, start_ind):
        self.matrix = copy.deepcopy(matrix)
        self.parent_lower_bound = 0
        self.depth = 1

        # keep track of the cities that we have visited
        self.path_set = set()
        self.path = []

        # when we travel between two cities, we will infinite-out the column and row
        # we will keep track of which columns and rows are infinited-out with these sets
        self.set_cols = set()
        self.set_rows = set()

        # select arbitrary start city

        # save the start city as the to_index
        self.to_city_ind = start_ind

        self.path_set.add(start_ind)
        self.path.append(cities[start_ind])

        self.reduce_first_state()

    def reduce_rows(self):
        # row reduce
        cost_reduce = 0
        for row in range(len(self.matrix)):
            # if row is not in our set of rows get min value
            if not self.set_rows.__contains__(row):
                # get min value for a row, if it's less than infinity subtract that value from rows
                min = self.get_min_row(row)
                if min > 0 and min != math.inf:
                    for col in range(len(self.matrix)):
                        cur_val = self.matrix[row][col]
                        self.matrix[row][col] = cur_val - min
                    cost_reduce += min
        return cost_reduce

    def reduce_cols(self, colTo=-1):
        # column reduce
        cost_reduce = 0
        can_get_min = True
        for row in range(len(self.matrix)):
            # if colTo = -1 we are reducing the initial state matrix
            if colTo == -1:
                can_get_min = True
            # if we have visited col before or our row if column of city we want to go to, don't get min
            elif self.set_cols.__contains__(row) or row is colTo:
                can_get_min = False
            if can_get_min:
                # get min value for a column, if it's less than infinity subtract that value from column
                min = self.get_min_col(row)
                if min > 0 and min != math.inf:
                    for col in range(len(self.matrix)):
                        cur_val = self.matrix[col][row]
                        self.matrix[col][row] = cur_val - min
                        # increment cost to reduce
                    cost_reduce += min
            # reset can get min to true
            can_get_min = True
        # return the cost to reduce
        return cost_reduce

    def get_min_row(self, row):
        row = self.matrix[row]
        min = math.inf
        # loop through rows
        for i in range(len(row)):
            # get num and update min if it is less
            num = row[i]
            if num == 0:
                return num
            elif num < min:
                min = num
        return min

    def get_min_col(self, col):
        min = math.inf
        # loop through columns
        for i in range(len(self.matrix)):
            # get num and update min if it is less
            num = self.matrix[i][col]
            if num == 0:
                return num
            elif num < min:
                min = num
        return min

    def get_key(self):
        # key for the queue is states lower bound * 2 / 2. This helps incorporate depth and lower bound to key
        return ((self.state_lower_bound * 2) / self.depth)

    # when a state with same lower bound and depth as another state on the queue, heapq.heappush
    # will compare the state objects using this function
    # __lt__ will select one to be placed 'ahead' of the other in the queue
    def __lt__(self, other):
        return True
