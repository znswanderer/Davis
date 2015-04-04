# from http://enja.org/2011/03/22/adventures-in-pyopencl-part-2-particles-with-pyopengl/

#basic glut setup learned from here:
#http://www.java2s.com/Open-Source/Python/Game-2D-3D/PyOpenGL/PyOpenGL-Demo-3.0.1b1/PyOpenGL-Demo/NeHe/lesson2.py.htm

from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import sys

#helper modules
from vector import Vec


class Window(object):

    def __init__(self, simulation, name, *args, **kwargs):
        #mouse handling for transforming scene
        self.mouse_down = False
        self.mouse_old = Vec([0., 0.])
        self.rotate = Vec([0., 0., 0.])
        self.translate = Vec([0., 0., 0.])
        self.initrans = Vec([0., 0., -2.])

        self.width = 640
        self.height = 480

        glutInit(sys.argv)
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)
        glutInitWindowSize(self.width, self.height)
        glutInitWindowPosition(0, 0)
        self.win = glutCreateWindow(name)

        #gets called by GLUT every frame
        glutDisplayFunc(self.draw)
        glutReshapeFunc(self.reshape);

        #handle user input
        glutKeyboardFunc(self.on_key)
        glutMouseFunc(self.on_click)
        glutMotionFunc(self.on_mouse_motion)
        
        #this will call draw every 30 ms
        glutTimerFunc(30, self.timer, 30)

        #setup OpenGL scene
        self.glinit()
        self.light()

        self.simulation = simulation
        simulation.start()
        glutMainLoop()

    def light(self):
        # from https://bazaar.launchpad.net/~mcfletch/pyopengl-demo/trunk/view/head:/PyOpenGL-Demo/redbook/scene.py
        light_ambient =  [0.5, 0.5, 0.5, 1.0]
        light_diffuse =  [.5, .5, .5, 1.0]
        light_specular =  [.5, .5, .5, 1.0]
        #  light_position is NOT default value
        light_position =  [1.0, 1.0, 1.0, 0.0]
        
        glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse)
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular)
        glLightfv(GL_LIGHT0, GL_POSITION, light_position)

    def glinit(self):
        glViewport(0, 0, self.width, self.height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(60., self.width / float(self.height), .1, 1000.)
        glMatrixMode(GL_MODELVIEW)
        glEnable(GL_DEPTH_TEST)

    ###GL CALLBACKS
    def timer(self, t):
        glutTimerFunc(t, self.timer, t)
        glutPostRedisplay()

    def on_key(self, *args):
        ESCAPE = '\033'
        if args[0] == ESCAPE or args[0] == 'q':
            sys.exit()
        elif args[0] == 's':
            self.simulation.save()
        elif args[0] == 'g':
            self.simulation.start_run()
        else:
            self.simulation.handle_key(args[0])

    def on_click(self, button, state, x, y):
        if state == GLUT_DOWN:
            self.mouse_down = True
            self.button = button
        else:
            self.mouse_down = False
        self.mouse_old.x = x
        self.mouse_old.y = y

    def on_mouse_motion(self, x, y):
        dx = x - self.mouse_old.x
        dy = y - self.mouse_old.y
        if self.mouse_down and self.button == 0: #left button
            self.rotate.x += dy * .2
            self.rotate.y += dx * .2
        elif self.mouse_down and self.button == 2: #right button
            self.translate.z -= dy * .01 
        self.mouse_old.x = x
        self.mouse_old.y = y
    ###END GL CALLBACKS

    def reshape(self, w, h):
        size = 1.5
        # from http://cs.lmu.edu/~ray/notes/openglexamples/
        glViewport(0, 0, w, h)
        glMatrixMode(GL_PROJECTION)
        aspect = float(w) / float(h)
        glLoadIdentity();
        if (w <= h):
            # width is smaller, so stretch out the height
            glOrtho(-size, size, -size/aspect, size/aspect, -10.0, 10.0)
        else:
            # height is smaller, so stretch out the width
            glOrtho(-size*aspect, size*aspect, -size, size, -10.0, 10.0)

    def draw(self):
        """Render the particles"""        
        glFlush()

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

        #handle mouse transformations
        glTranslatef(self.initrans.x, self.initrans.y, self.initrans.z)
        glRotatef(self.rotate.x, 1, 0, 0)
        glRotatef(self.rotate.y, 0, 1, 0) #we switched around the axis so make this rotate_z
        glTranslatef(self.translate.x, self.translate.y, self.translate.z)
        
        #render the particles
        self.simulation.render()

        glutSwapBuffers()




