# lab_7_examples.py
#
# Author: Quinn
#
# These classes were written as part of Physics 102: Computational Laboratory
# in Physics, Fall 2014.

# import required libraries
import numpy as np
import matplotlib as mpl
# Needed for making the animation work on Mac
import platform
if platform.system() != 'Windows' and platform.system() != 'Linux':
    mpl.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from random import random
from random import gauss



# Example 1
class flocking_sim_1:
    # This function is run when a class is instantiated. It sets up the class and
    # initializes the positions and velocities of the fish.
    def __init__(self, size_of_school=5, total_time=10):
        # Create a class element to store the size of the school
        self.size_of_school = size_of_school
        
        # Create class elements for the postions and velocities of the fish
        self.positions = np.zeros((total_time, size_of_school, 2))        
        self.velocities = np.zeros((total_time, size_of_school, 2))
        
        # Create a class element for the times
        self.times = np.arange(total_time)
        
        # Initialize the positions of the fish in cartesian cordinates (x,y)
        for i in xrange(size_of_school):
            self.positions[0,i,0] = 20.*random()
            self.positions[0,i,1] = 20.*random()

        # Initialize the velocities of the fish in polar cordinates (v, theta)
        for i in xrange(size_of_school):
            self.velocities[0,i,0] = gauss(0.5, 0.025)
            self.velocities[0,i,1] = 2.* np.pi*random()

    # This function obtains the new postions of the fish using their velocities
    def step_fish_pos_to(self, time):
        # Calculate the change in the postions using the velocities 
        #       delta_x = v_x \delta_t (with delta_t = 1)
        #       delta_y = v_y \delta_t (with delta_t = 1)
        # We also need to covert the polar velocities to cartesian cordinates
        delta_x = self.velocities[time-1,:,0]*np.cos(self.velocities[time-1,:,1])
        delta_y = self.velocities[time-1,:,0]*np.sin(self.velocities[time-1,:,1])
        
        # The new position coordinates are the old cordinates plus the changes
        # in position. The %20 makes the positions wrap. If the fish swims
        # off the right, it shows up on left. If the fisth swims off the top
        # it shows up on the bottom.
        self.positions[time,:,0] = (self.positions[time-1,:,0] + delta_x)%20
        self.positions[time,:,1] = (self.positions[time-1,:,1] + delta_y)%20
    
    # This function updates the fish velocities. This function is what
    # determines the simulation.
    def step_fish_vel_to(self, time):
        # The magnitude of the fish velocity remains unchanged
        self.velocities[time,:,0] = self.velocities[time-1,:,0]
        
        # The new velocity angle is the old angle plus a small change
        self.velocities[time,:,1] = self.velocities[time-1,:,1] + 0.07*np.pi
    
    # This function carries out the simulation
    def run_sim(self):
        # For every time in the time array
        for time in self.times[1:]:
            # Get the new fish positions
            self.step_fish_pos_to(time)
            
            # Get the new fish velocities
            self.step_fish_vel_to(time)

    # Function to create an animation of the simulation
    def animate_sim(self):
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_axes((0.05, 0.05, 0.9, 0.9), xlim=(0, 20), ylim=(0,20),
                          autoscale_on=False, aspect='equal', xticks=[], yticks=[])

        self.dots = ax.scatter([], [])

        def init():
            self.dots.set_offsets([[]])
            return self.dots,
        
        def animate(i):
            self.dots.set_offsets(self.positions[i])
            return self.dots,

        ani = animation.FuncAnimation(fig, animate, self.times,
                                      interval=50, blit=True, init_func=init)

        plt.show()

# Example 2
class flocking_sim_2:
    # This function is run when the class is instantiated. It sets up the class and
    # initializes the positions and velocities of the fish.
    def __init__(self, size_of_school=5, total_time=10):
        # Create a class element to store the size of the school
        self.size_of_school = size_of_school
        
        # Create class elements for the postions and velocities of the fish
        self.positions = np.zeros((total_time, size_of_school, 2))        
        self.velocities = np.zeros((total_time, size_of_school, 2))
        
        # Create a class element for the times
        self.times = np.arange(total_time)
        
        # Initialize the positions of the fish in cartesian cordinates (x,y)
        for i in xrange(size_of_school):
            self.positions[0,i,0] = 20.*random()
            self.positions[0,i,1] = 20.*random()

        # Initialize the velocities of the fish in polar cordinates (v, theta)
        for i in xrange(size_of_school):
            self.velocities[0,i,0] = gauss(0.5, 0.025)
            self.velocities[0,i,1] = 2.* np.pi*random()

    # This function obtains the new postions of the fish using their velocities
    def step_fish_pos_to(self, time):
        # Calculate the change in the postions using the velocities 
        #       delta_x = v_x \delta_t (with delta_t = 1)
        #       delta_y = v_y \delta_t (with delta_t = 1)
        # We also need to covert the polar velocities to cartesian cordinates
        delta_x = self.velocities[time-1,:,0]*np.cos(self.velocities[time-1,:,1])
        delta_y = self.velocities[time-1,:,0]*np.sin(self.velocities[time-1,:,1])
        
        # The new position coordinates are the old cordinates plus the changes
        # in position. The %20 makes the positions wrap. If the fish swims
        # off the right, it shows up on left. If the fisth swims off the top
        # it shows up on the bottom.
        self.positions[time,:,0] = (self.positions[time-1,:,0] + delta_x)%20
        self.positions[time,:,1] = (self.positions[time-1,:,1] + delta_y)%20
    
    # This function updates the fish velocities. This function is what
    # determines the simulation.
    def step_fish_vel_to(self, time):
        # The magnitude of the fish velocity remains unchanged
        self.velocities[time,:,0] = self.velocities[time-1,:,0]
        
        # Loop over fish indices
        for i in xrange(self.size_of_school):
            # The new velocity angle is the weighted (0.05 to 0.95) average
            # of the fishes current velocity angle and the current velocity
            # angle of a 'neighbor' fish.
            self.velocities[time,i,1] = (0.1*self.velocities[time-1,(i+1)%self.size_of_school,1]
                                         + 0.9*self.velocities[time-1,i,1])
    
    # This function carries out the simulation
    def run_sim(self):
        # For every time in the time array
        for time in self.times[1:]:
            # Get the new fish positions
            self.step_fish_pos_to(time)
            
            # Get the new fish velocities
            self.step_fish_vel_to(time)

    # Function to create an animation of the simulation
    def animate_sim(self):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_axes((0.05, 0.05, 0.9, 0.9), xlim=(0, 20), ylim=(0,20),
                          autoscale_on=False, aspect='equal', xticks=[], yticks=[])

        self.dots = ax.scatter([], [])

        def init():
            self.dots.set_offsets([[]])
            return self.dots,
        
        def animate(i):
            self.dots.set_offsets(self.positions[i])
            return self.dots,

        ani = animation.FuncAnimation(fig, animate, self.times,
                                      interval=50, blit=True, init_func=init)

        plt.show()


# Example 3
class flocking_sim_3:
    # This function is run when a class is instantiated. It sets up the class and
    # initializes the positions and velocities of the fish and a shark.
    def __init__(self, size_of_school=5, total_time=10):
        # Create a class element to store the size of the school
        self.size_of_school = size_of_school
        
        # Create class elements for the postions and velocities of the fish
        self.positions = np.zeros((total_time, size_of_school, 2))        
        self.velocities = np.zeros((total_time, size_of_school, 2))
        
        # Create class elements for the positions and velocities of a shark
        self.shark_pos = np.zeros((total_time, 2))        
        self.shark_vel = np.zeros((total_time, 2))
        
        # Create a class element for the times
        self.times = np.arange(total_time)
        
        # Initialize the positions of the fish in cartesian cordinates (x,y)
        for i in xrange(size_of_school):
            self.positions[0,i,0] = 20.*random()
            self.positions[0,i,1] = 20.*random()

        # Initialize the velocities of the fish in polar cordinates (v, theta)
        for i in xrange(size_of_school):
            self.velocities[0,i,0] = gauss(0.5, 0.025)
            self.velocities[0,i,1] = 2.* np.pi*random()
        
        # Initializes the positions and velocity of a shark
        self.shark_pos[0,0] = 20.*random()
        self.shark_pos[0,1] = 20.*random()
        
        self.shark_vel[0,0] = 0.15
        self.shark_vel[0,1] = 2.* np.pi*random()

    # This function obtains the new postions of the fish using their velocities
    def step_fish_pos_to(self, time):
        # Calculate the change in the postions using the velocities 
        #       delta_x = v_x \delta_t (with delta_t = 1)
        #       delta_y = v_y \delta_t (with delta_t = 1)
        # We also need to covert the polar velocities to cartesian cordinates
        delta_x = self.velocities[time-1,:,0]*np.cos(self.velocities[time-1,:,1])
        delta_y = self.velocities[time-1,:,0]*np.sin(self.velocities[time-1,:,1])
        
        # The new position coordinates are the old cordinates plus the changes
        # in position. The %20 makes the positions wrap. If the fish swims
        # off the right, it shows up on left. If the fish swims off the top
        # it shows up on the bottom.
        self.positions[time,:,0] = (self.positions[time-1,:,0] + delta_x)%20
        self.positions[time,:,1] = (self.positions[time-1,:,1] + delta_y)%20
    
    # This function updates the fish velocities. This function is what
    # determines the simulation.
    def step_fish_vel_to(self, time):
        # The magnitude of the fish velocity remains unchanged
        self.velocities[time,:,0] = self.velocities[time-1,:,0]
        
        # Loop over fish indices
        for i in xrange(self.size_of_school):
            # Because the edges wrap, we need to consider whether
            # the fish is 'near' the shark on the opposite edge.
            # Define an array of translation vectors.
            vectors = np.array([[-20, -20], [-20, 0], [-20, 20],
                         [0, -20], [0, 0], [0, 20],
                         [20, -20], [20, 0], [20, 20]])
            
            # Define a variable assuming the fish is not close to the shark
            close = False
            
            for vec in vectors:
                # Calculate the current distance from the shark
                delta = self.positions[time, i] - self.shark_pos[time-1]-vec
            
                # If the fish is near the shark
                if np.linalg.norm(delta) < 4:
                    # Set close to true
                    close = True
                    
                    # Move away from the shark
                    self.velocities[time,i,1] = np.arctan2(delta[1], delta[0])
                    
                    # Break out of for loop
                    break
                
            if not close:
                # Move in circles
                self.velocities[time,i,1] = self.velocities[time-1,i,1] + 0.07*np.pi     
    
    # This function updtates the position of the shark
    def step_shark_pos_to(self, time):
        # Calculate the change in the postions using the velocities 
        #       delta_x = v_x \delta_t (with delta_t = 1)
        #       delta_y = v_y \delta_t (with delta_t = 1)
        # We also need to covert the polar velocities to cartesian cordinates
        delta_x = self.shark_vel[time-1, 0]*np.cos(self.shark_vel[time-1, 1])
        delta_y = self.shark_vel[time-1, 0]*np.sin(self.shark_vel[time-1, 1])
        
        # The new position coordinates are the old cordinates plus the changes
        # in position. The %20 makes the positions wrap. If the shark swims
        # off the right, it shows up on left. If the shark swims off the top
        # it shows up on the bottom.
        self.shark_pos[time,0] = (self.shark_pos[time-1,0] + delta_x)%20
        self.shark_pos[time,1] = (self.shark_pos[time-1,1] + delta_y)%20
        
    # This function updates the velocity of the shark
    def step_shark_vel_to(self, time):
        # The magnitude of the shark velocity remains unchanged
        self.shark_vel[time, 0] = self.shark_vel[time-1, 0]
        
        # The shark will randomly change directions
        if random() > 0.95:
            self.shark_vel[time, 1] = self.shark_vel[time-1, 1] + 0.25*np.pi*(random() - 1)
        else:
            self.shark_vel[time, 1] = self.shark_vel[time-1, 1]
            
    # This function carries out the simulation
    def run_sim(self):
        # For every time in the time array
        for time in self.times[1:]:
            # Get the new fish positions
            self.step_fish_pos_to(time)
            
            # Get the new fish velocities
            self.step_fish_vel_to(time)
            
            # Get the new shark positions
            self.step_shark_pos_to(time)
            
            # Get the new shark velocities
            self.step_shark_vel_to(time)

    # Function to create an animation of the simulation
    def animate_sim(self):
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_axes((0.05, 0.05, 0.9, 0.9), xlim=(0, 20), ylim=(0,20),
                          autoscale_on=False, aspect='equal', xticks=[], yticks=[])

        self.fish_dots = ax.scatter([], [])
        self.shark_dot = ax.scatter([], [], c='red', s=150)

        def init():
            self.fish_dots.set_offsets([[]])
            self.shark_dot.set_offsets([[]])
            return self.fish_dots, self.shark_dot
        
        def animate(i):
            self.fish_dots.set_offsets(self.positions[i])
            self.shark_dot.set_offsets(self.shark_pos[i])
            return self.fish_dots, self.shark_dot

        ani = animation.FuncAnimation(fig, animate, self.times,
                                      interval=50, blit=True, init_func=init)

        plt.show()

