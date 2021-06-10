#!/usr/bin/env python

import numpy as np

class Quaternion:
    def __init__(self, q):
        self.q = q
    def __repr__(self):
        return '%6.3f %6.3f %6.3f %6.3f'%tuple(self.q)
    @classmethod
    def from_axis_and_angle(cls, axis, angle):
        v = np.array(axis)
        v /= np.sqrt(np.dot(v,v))
        q = np.zeros(4)
        q[0] = np.cos(angle/2.0)
        q[1:] = v*np.sin(angle/2.0)
        return cls(q)
    def rotate(self):
        if 'R' in dir(self):
            return self.R
        #
        self.R = np.zeros((3,3))
        #
        self.R[0][0] = self.q[0]**2 + self.q[1]**2 - self.q[2]**2 - self.q[3]**2
        self.R[0][1] = 2.0*(self.q[1]*self.q[2] - self.q[0]*self.q[3])
        self.R[0][2] = 2.0*(self.q[1]*self.q[3] + self.q[0]*self.q[2])
        #
        self.R[1][0] = 2.0*(self.q[1]*self.q[2] + self.q[0]*self.q[3])
        self.R[1][1] = self.q[0]**2 - self.q[1]**2 + self.q[2]**2 - self.q[3]**2
        self.R[1][2] = 2.0*(self.q[2]*self.q[3] - self.q[0]*self.q[1])
        #
        self.R[2][0] = 2.0*(self.q[1]*self.q[3] - self.q[0]*self.q[2])
        self.R[2][1] = 2.0*(self.q[2]*self.q[3] + self.q[0]*self.q[1])
        self.R[2][2] = self.q[0]**2 - self.q[1]**2 - self.q[2]**2 + self.q[3]**2
        return self.R
    def __mul__(self, othr):
        s = np.zeros(4, dtype=float)
        s[0] = self.q[0]*othr.q[0] - np.dot(self.q[1:], othr.q[1:])
        tmp = np.cross(self.q[1:], othr.q[1:])
        s[1:] = self.q[0]*othr.q[1:] + othr.q[0]*self.q[1:] + tmp
        return Quaternion(s)
    @property
    def euler(self):
        q = self.q
        phi = np.arctan2( 2.0*(q[0]*q[1] + q[2]*q[3]), 1.0-2.0*(q[1]**2+q[2]**2) )
        theta = np.arcsin( 2.0*(q[0]*q[2] - q[3]*q[1]) )
        psi = np.arctan2( 2.0*(q[0]*q[3] + q[1]+q[2]), 1.0-2.0*(q[2]**2+q[3]**2) )
        return np.array([phi, theta, psi], dtype=float)

def sample_fibonacci_sphere(n_sample=1, randomize=True):
    if randomize:
        rnd = np.random.random() * n_sample
    else:
        rnd = 1.0
    #
    point_s = []
    offset = 2.0 / n_sample
    increment = np.pi * (3.0 - np.sqrt(5.0))
    #
    for i in range(n_sample):
        z = (i*offset - 1.0) + (offset / 2.0)
        r = np.sqrt(1 - z**2)
        phi = ((i + rnd)%n_sample) * increment
        x = np.cos(phi) * r
        y = np.sin(phi) * r
        point_s.append((x,y,z))
    #
    point_s = np.array(point_s, dtype=float)
    #
    return point_s
