#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np
from TSPClasses import *
from State import *
import heapq
import itertools



class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution, 
		time spent to find solution, number of permutations tried during search, the 
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''
	
	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation(ncities)
			route = []
			# Now build the route using the random permutation
			for i in range(ncities):
				route.append(cities[perm[i]])
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results


	''' <summary>
		This is the entry point for the greedy solver, which you must implement for 
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

	# Time Complexity: O(n^3) may need to pick up to n start cities
	# Space Complexity: O(n) containing path
	def greedy(self,time_allowance=60.0):

		# create results
		results = {}

		# path array to keep track of path in order
		path = []

		# get cities
		cities = self._scenario.getCities()

		# create visited set to keep track of cities we have already seen
		visited = set()

		# random stat index
		start_index = random.randint(0, len(cities) - 1)

		index = start_index

		start_time = time.time()

		found_path = False
		valid_path = False
		shortest_edge = 0
		while not found_path and time.time()-start_time < time_allowance:

			# if we found a shortest edge add it to our path, visited set, and get the next city
			if shortest_edge is not np.inf:
				curr_city = cities[index]
				path.append(cities[index])
				visited.add(cities[index])
			else:
				# break from loop, no out going edges from city
				found_path = True
				break

			# set shortest edge to infinity, see if we can update it to something smaller
			shortest_edge = np.inf

			# if we have visited all of our cities
			if len(cities) == len(visited):
				first_city = cities[start_index]
				final_edge = curr_city.costTo(first_city)
				# see if path from final edge to start edge is possible
				if final_edge is not np.inf:
					valid_path = True
				found_path = True
				break
			for city in cities:
				# if we are not on our current city and we havent visited city, get edge
				if city._name is not curr_city._name and not visited.__contains__(city):
					edge = curr_city.costTo(city)
					if edge is not np.inf:
						# base case
						if shortest_edge == np.inf:
							shortest_edge = edge
							index = city._index

						elif edge < shortest_edge:
							shortest_edge = edge
							index = city._index

		end_time = time.time()
		# if we ran out of time
		if not found_path:
			results['cost'] = math.inf
			results['time'] = end_time - start_time
			results['count'] = 0
			results['soln'] = None
			results['max'] = None
			results['total'] = None
			results['pruned'] = None
			return results
		# if we found a path and it was valid
		elif found_path and valid_path:
			bssf = TSPSolution(path)
			results['cost'] = bssf.cost
			results['time'] = end_time - start_time
			results['count'] = 1
			results['soln'] = bssf
			results['max'] = None
			results['total'] = None
			results['pruned'] = None
			return results
		# start over with new city
		elif not valid_path and found_path:
			return self.greedy()


	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints:
		max queue size, total number of states created, and number of pruned states.</returns>
	'''

	def branchAndBound( self, time_allowance=60.0 ):
		start_time = time.time()
		# values to return
		self.solution_found = 0
		self.total_pruned = 0
		self.states_created = 0
		self.max_queue = 0
		results = {}
		# get cities
		cities = self._scenario.getCities()

		# initialize queue and heapify
		self.queue = []
		heapq.heapify(self.queue)

		# get bssf from greedy
		self.bssf = self.greedy()['soln']

		# initialize and populate matrix
		matrix = [[math.inf for x in range(len(cities))] for y in range(len(cities))]
		for i in range(len(cities)):
			for j in range(len(cities)):
				matrix[i][j] = cities[i].costTo(cities[j])

		# pick random starting index
		self.start_index = random.randint(0, len(cities) - 1)

		# create initial state and reduce it
		#initial_state = State(matrix, None, self.start_index, cities)
		initial_state = State(None, None, None)
		initial_state.set_first_state(matrix, cities, self.start_index)

		# push initial state to queue with key of 1
		heapq.heappush(self.queue, (initial_state.get_key(), initial_state))


		while len(self.queue) != 0 and time.time()-start_time <= time_allowance:

			# update max size of queue
			if len(self.queue) > self.max_queue:
				self.max_queue = len(self.queue)

			# get parent state from the top of the queue
			key, parent_state = heapq.heappop(self.queue)

			# look at states to visit
			if parent_state.state_lower_bound < self.bssf.cost:
				self.create_states(parent_state)

		end_time = time.time()

		# return results
		results['cost'] = self.bssf.cost
		results['time'] = end_time - start_time
		results['count'] = self.solution_found
		results['soln'] = self.bssf
		results['max'] = self.max_queue
		results['total'] = self.states_created
		results['pruned'] = self.total_pruned
		return results

	# Time:
	# Space
	def create_states(self, parent_state):
		# if we have visited all cities, check if we can get back to the first
		if len(parent_state.path_set) == len(self._scenario.getCities()):
			last_city = parent_state.path[-1]
			start_city = parent_state.path[0]
			if last_city.costTo(start_city) is not np.inf:
				# update solution if we have found a better one
				solution = TSPSolution(parent_state.path)
				if solution.cost < self.bssf.cost:
					self.bssf = solution
					self.prune()
					# increment number of solutions found
					self.solution_found += 1
			# if path back to first city is not possible we are done
			else:
				return
		else:
			for i in range(len(self._scenario.getCities())):
				# if current city hasn't been visited
				if i not in parent_state.path_set:
					new_state = State(parent_state, i, self._scenario.getCities())
					# increment number of states created
					self.states_created += 1
					# if new created state has cost less than bssf add to queue with key
					if new_state.state_lower_bound < self.bssf.cost and new_state.state_lower_bound != math.inf:
						key = new_state.get_key()
						heapq.heappush(self.queue, (key, new_state))
					else:
						self.total_pruned += 1

	def prune(self):
		for tuple in self.queue:
			key, state = tuple
			# remove state from queue if its lower bound is greater than bssf
			if state.state_lower_bound >= self.bssf.cost:
				self.queue.remove(tuple)
				self.total_pruned += 1

	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''
		
	def fancy( self,time_allowance=60.0 ):
		pass
		



