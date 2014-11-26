# lab_7_examples.py
#
# Author:
# Student ID:
#
# This library was written as part of Physics 102: Computational Laboratory
# in Physics, Fall 2014.

# import required libraries
import math
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


# Exercise
def to_polar(x, y):
    r = math.sqrt(x**2 + y**2)
    theta = math.atan2(y, x)
    return [r, theta]

def to_rect(r, theta):
    x = r * math.cos(theta)
    y = r * math.sin(theta)
    return [x, y]

class my_flocking_sim:
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
            self.positions[0,i,0] = (19./np.sqrt(2.))*random() - 9.5/np.sqrt(2.)
            self.positions[0,i,1] = (19./np.sqrt(2.))*random() - 9.5/np.sqrt(2.)

        # Initialize the velocities of the fish in polar cordinates (v, theta)
        for i in xrange(size_of_school):
            self.velocities[0,i,0] = 0.2
            self.velocities[0,i,1] = 2.* np.pi*random()

        # Initializes the position of the shark
        self.shark_pos[0,0] = (19./np.sqrt(2.))*random() - 9.5/np.sqrt(2.)
        self.shark_pos[0,1] = (19./np.sqrt(2.))*random() - 9.5/np.sqrt(2.)

        # Initializes the velocity of the shark
        self.shark_vel[0,0] = 0.05
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
        # in position.
        self.positions[time,:,0] = self.positions[time-1,:,0] + delta_x
        self.positions[time,:,1] = self.positions[time-1,:,1] + delta_y

        # Reset the updated array
        self.already_updated = [False]*self.size_of_school

    # This function updates the fish velocities. This function is what
    # determines the simulation.
    def step_fish_vel_to(self, time):
        # The magnitude of the fish velocity remains unchanged
        self.velocities[time,:,0] = self.velocities[time-1,:,0]

        # Loop over fish indices
        for i in xrange(self.size_of_school):

            # This code turns the fish around at the edge of the tank
            if np.linalg.norm(self.positions[time,i]) > 9.5:
                position_angle = np.arctan2(self.positions[time,i,1], self.positions[time,i,0])
                self.velocities[time,i,1] = 2.0*position_angle - self.velocities[time-1,i,1] + np.pi
                self.already_updated[i] = True

            # If the fish has not had a colision with the edge of the tank
            else:
                # Because the edges wrap, we need to consider whether
                # the fish is 'near' the shark on the opposite edge.
                # Define an array of translation vectors.
                vectors = np.array([[-20, -20], [-20, 0], [-20, 20],
                         [0, -20], [0, 0], [0, 20],
                         [20, -20], [20, 0], [20, 20]])

                # Define a variable assuming the fish is not close to the shark
                close = False

                # Loop over vectors
                for vec in vectors:
                    # Calculate the current distance from the shark
                    delta = self.positions[time, i] - self.shark_pos[time-1]-vec

                    # If the fish is near the shark
                    if np.linalg.norm(delta) < 3:
                        # Move away from the shark
                        self.velocities[time,i,1] = np.arctan2(delta[1], delta[0])
                        close = True

                        # Break out of for loop
                        break

                if not close:
                    #       THIS IS THE BLOCK OF CODE YOU NEED TO CHANGE
                    #
                    # When not near a shark the fish just swim in straigth lines.
                    # This is not flocking behavoir. The goal of this week's
                    # exercise is for you to come up with a rule for how the
                    # fish change their velocities over time. Once you have an
                    # idea for a rule, code it up and see what happens. If it
                    # looks like flocking, you are done. If not, try a new rule
                    # until you get something that looks like flocking. You know
                    # you have flocking behavoir if the fish behave like a single
                    # organism.
                    # if self.turning_progress[i] < 200 and self.turning_progress[i] != 0:
                    #     #self.velocities[time,i,0] = self.velocities[time-1,i,0]
                    #     #self.velocities[time,i,1] = self.velocities[time-1,i,1]
                    #     self.turning_progress[i] += 1
                    #     return
                    # self.turning_progress[i] = 0
                    alignment = self.aligmnent(time, i)
                    cohesion = self.cohesion(time, i)
                    separation = self.separation(time, i)
                    r_theta_pair = to_polar(alignment[0] + cohesion[0] + separation[0], alignment[1] + cohesion[1] + separation[1])
                    self.velocities[time,i,0] = r_theta_pair[0]
                    self.velocities[time,i,1] = r_theta_pair[1]
                    if self.velocities[time,i,0] < 0.0001:
                        print "calculated speed is not enough"
                        self.velocities[time,i,0] = self.velocities[time-1,i,0]

    # This function updtates the position of the shark
    def step_shark_pos_to(self, time):
        # Calculate the change in the postions using the velocities
        #       delta_x = v_x \delta_t (with delta_t = 1)
        #       delta_y = v_y \delta_t (with delta_t = 1)
        # We also need to covert the polar velocities to cartesian cordinates
        delta_x = self.shark_vel[time-1, 0]*np.cos(self.shark_vel[time-1, 1])
        delta_y = self.shark_vel[time-1, 0]*np.sin(self.shark_vel[time-1, 1])

        # The new position coordinates are the old cordinates plus the changes
        # in position.
        self.shark_pos[time,0] = self.shark_pos[time-1,0] + delta_x
        self.shark_pos[time,1] = self.shark_pos[time-1,1] + delta_y

    # This function updates the velocity of the shark
    def step_shark_vel_to(self, time):
        # The magnitude of the shark velocity remains unchanged
        self.shark_vel[time, 0] = self.shark_vel[time-1, 0]

        # This code turns the shark around at the edge of the tank
        if np.linalg.norm(self.shark_pos[time]) > 9.5:
            position_angle = np.arctan2(self.shark_pos[time,1], self.shark_pos[time,0])
            self.shark_vel[time,1] = 2.0*position_angle - self.shark_vel[time-1,1] + np.pi

        # The shark will randomly change directions
        elif random() > 0.95:
            self.shark_vel[time, 1] = self.shark_vel[time-1, 1] + 0.25*np.pi*(random() - 1)
        else:
            self.shark_vel[time, 1] = self.shark_vel[time-1, 1]

    def aligmnent(self, time, fishIndex):
        x = 0.0
        y = 0.0
        neibourCount = 0
        for i in xrange(self.size_of_school):
            if fishIndex != i:
                dist = math.sqrt((self.positions[time,i,0] - self.positions[time,fishIndex,0])**2 + (self.positions[time,i,1] - self.positions[time,fishIndex,1])**2)
                print "distance = {0}".format(dist)
                if (dist < 1):
                    x_y_pair = to_rect(self.velocities[time, i, 0], self.velocities[time, i, 1])
                    x += x_y_pair[0]
                    y += x_y_pair[1]
                    neibourCount += 1
        if (neibourCount != 0):
            x /= neibourCount
            y /= neibourCount
            # x = normalize(x)
            # y = normalize(y)
            print "alignment: {0}, {1}".format(x, y)
        return (x, y)

    def cohesion(self, time, fishIndex):
        x = 0.0
        y = 0.0
        neibourCount = 0
        for i in xrange(self.size_of_school):
            if fishIndex != i:
                dist = math.sqrt((self.positions[time,i,0] - self.positions[time,fishIndex,0])**2 + (self.positions[time,i,1] - self.positions[time,fishIndex,1])**2)
                print "distance = {0}".format(dist)
                if (dist < 1):
                    x += self.positions[time, i, 0]
                    y += self.positions[time, i, 1]
                    neibourCount += 1
        if (neibourCount != 0):
            x /= neibourCount
            y /= neibourCount
            x = x - self.positions[time,fishIndex,0]
            y = y - self.positions[time,fishIndex,1]
            # x = normalize(x)
            # y = normalize(y)
            print "cohesion: {0}, {1}".format(x, y)
        return (x, y)

    def separation(self, time, fishIndex):
        x = 0.0
        y = 0.0
        neibourCount = 0
        for i in xrange(self.size_of_school):
            if fishIndex != i:
                dist = math.sqrt((self.positions[time,i,0] - self.positions[time,fishIndex,0])**2 + (self.positions[time,i,1] - self.positions[time,fishIndex,1])**2)
                print "distance = {0}".format(dist)
                if (dist < 1):
                    x += self.positions[time, i, 0] - self.positions[time, fishIndex, 0]
                    y += self.positions[time, i, 1] - self.positions[time, fishIndex, 1]
                    neibourCount += 1
        if (neibourCount != 0):
            x /= neibourCount
            y /= neibourCount
            x = -x
            y = -y
            print "separation: {0}, {1}".format(x, y)
        return (x, y)

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
        ax = fig.add_axes((0.05, 0.05, 0.9, 0.9), xlim=(-10.1, 10.1), ylim=(-10.1, 10.1),
                          autoscale_on=False, aspect='equal', xticks=[], yticks=[])

        circ = plt.Circle((0,0), radius=10, color='black', fill=False)
        fig.gca().add_artist(circ)

        self.fish_dots = ax.scatter([], [], c='cyan')
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
