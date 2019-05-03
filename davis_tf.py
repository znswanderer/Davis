import numpy as np
import tensorflow as tf
import scipy.spatial.distance as dist
import random
import math
import time
import fibonacci_sphere

FLOAT_TYPE = np.float_

class World:
    
    def __init__(self, num_particles, dt=0.1, gamma=0.0, cutoff=0.0, binning=None):
        self.N = num_particles
        self.dt = dt
        self.gamma = gamma
        self.cutoff = cutoff
        self.binning = binning

    def build_graph(self):
        self.graph = tf.Graph()
        with self.graph.as_default():
            self.initial_conds()
            self.init_op = tf.global_variables_initializer()
            self.step = self.timestep_graph()

        self.sess = tf.Session(graph=self.graph)
        self.sess.run(self.init_op)

        #tf.summary.FileWriter(logdir='./logs/', graph=self.graph)

    def initial_conds(self):
        self.tf_r = tf.Variable(FLOAT_TYPE(fibonacci_sphere.create_points(self.N)))
        self.tf_v = tf.Variable(np.zeros(shape=(self.N, 3), dtype=FLOAT_TYPE))
        self.tf_a = tf.Variable(np.zeros(shape=(self.N, 3), dtype=FLOAT_TYPE))

    def get_3dPositions(self):
        self.update_arrays()
        return self.r

    def timestep(self):
        self.sess.run(self.step)

    def update_arrays(self):
        self.r = self.tf_r.eval(self.sess)
        self.v = self.tf_v.eval(self.sess)
        self.a = self.tf_a.eval(self.sess)

    def timestep_graph(self):
        r_, v_ = self.advance()
        a_ = self.calc_forces(r_, v_)
        r_, v_ = self.correct(r_, v_, a_)
        step = tf.group(
            self.tf_r.assign(r_),
            self.tf_v.assign(v_),
            self.tf_a.assign(a_))
        return step

    def advance(self):
        v_ = self.tf_v + (self.dt / 2.0) * self.tf_a
        r_ = self.tf_r + self.dt * v_

        # RATTLE
        r0_dot_r = tf.einsum('ij, ij->i', self.tf_r, r_)
        r0_sq = tf.einsum('ij, ij->i', self.tf_r, self.tf_r)

        # TODO beware of negative sqrt arguments!
        # (But in this case the timestep dt is way too large anyway)        
        lambda_ = tf.sqrt(1.0 - r0_sq + r0_dot_r * r0_dot_r) - r0_dot_r
        lambda_ = tf.reshape(lambda_, (-1, 1))
        dr = lambda_ * self.tf_r

        r_ = r_ + dr
        v_ = v_ + dr / self.dt

        return r_, v_

    def correct(self, r_, v_, a_):
        v_ = v_ + (self.dt / 2.0) * a_
        # RATTLE_v
        lambda_ = tf.einsum('ij, ij->i', v_, r_)
        lambda_ = tf.reshape(lambda_, (-1, 1))
        v_ = v_ - lambda_ * r_
        return r_, v_

    def _calc_forces(self, this_r, this_v, other_r, other_v, fill_diagonal=False):
        dr = tf.subtract(this_r[:, tf.newaxis], other_r)

        dr2 = tf.reduce_sum(dr**2, axis=2)
        if fill_diagonal:
            dr2 = tf.matrix_set_diag(dr2, tf.ones(tf.shape(dr2)[0:-1], dtype=dr2.dtype) * (self.cutoff**2 + 1))
        mask = dr2 < (self.cutoff**2)
        masked_dr2 = tf.boolean_mask(dr2, mask)
        masked_dr = tf.boolean_mask(dr, mask)

        d = tf.sqrt(masked_dr2)
        dr_normed = masked_dr / tf.expand_dims(d, -1)

        # Coulomb force
        force_mag = (1.0 / masked_dr2) - 1.0 / self.cutoff**2
        masked_forces_pairs = dr_normed * tf.expand_dims(force_mag, -1)

        # Damping
        # Note: the damping is done with v(t + dt/2).  To be correct, we
        # might should use v(t+dt) (which is not available at this
        # point). But as we just model this dissipative force in a
        # handwaving sort anyway this should not matter. The damping has
        # not real physical meaning and is just introduced to cool the
        # system down.
        dv = tf.subtract(this_v[:, tf.newaxis], other_v)
        masked_dv = tf.boolean_mask(dv, mask)
        masked_damping_pairs = masked_dv * self.gamma

        forces_pairs = tf.scatter_nd(tf.where(mask), masked_forces_pairs + masked_damping_pairs, dr.shape)
        forces = tf.reduce_sum(forces_pairs, axis=1)
        if not fill_diagonal:
            other_forces = tf.reduce_sum(forces_pairs, axis=0)
            return forces, other_forces
        else:
            return forces

    def calc_forces(self, r_, v_):
        forces = self._calc_forces(r_, v_, r_, v_, fill_diagonal=True)
        return forces



if __name__ == '__main__':
    #tf.InteractiveSession()
    world = World(10, cutoff=0.5)
    world.build_graph()
    for i in range(10000):
        world.timestep()


