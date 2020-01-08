"""
Written by dsw7@sfu.ca
Here I animate the workings of the Met-aromatic algorithm

------------------
To run the program:
$ python3 /path/to/OpenGLMetAromatic.py

==============
Install notes:
==============

I think there may be a version incompatibility somewhere.
The program appears to run fine but I get a few warnings
at compilation:

OpenGL_accelerate seems to be installed, but unable to import error checking entry point!
Unable to load ArrayDatatype accelerator from OpenGL_accelerate
Unable to load converters accelerators (wrapper, arraydatatype) from OpenGL_accelerate
Unable to load arrayhelpers accelerator from OpenGL_accelerate
OpenGL_accelerate seems to be installed, but unable to import expected wrapper entry points!
Unable to load VBO accelerator from OpenGL_accelerate

Path to this Python interpreter & libs:
/Library/Frameworks/Python.framework/Versions/3.7/bin/python3
/Library/Frameworks/Python.framework/Versions/3.7/lib/python3.7/site-packages/
"""

import pygame
from pygame.locals   import *
from OpenGL.GL       import *
from OpenGL.GLU      import *
from OpenGL.GLUT     import *
from retrieve_1rcy   import import_coords
from numpy           import array, linalg
from re              import search
from itertools       import groupby, chain              
from operator        import itemgetter 
from collections     import deque    
from time            import time             

CUTOFF            = 6.0  # Angstroms
RES               = r'MET|PHE|TYR|TRP'
ATOMS_MET         = r'CE|SD|CG' 
ATOMS_TYR         = r'CD1|CE1|CZ|CG|CD2|CE2'   
ATOMS_TRP         = r'CD2|CE3|CZ2|CH2|CZ3|CE2'  
ATOMS_PHE         = r'CD1|CE1|CZ|CG|CD2|CE2'

DISPLAY_DIMENSION = 800
SLICES_STACKS     = (20, 20)
LP_WIDTHS         = 5.0  # non default GLfloat that specifies lone pair width

# control timing
# --------------
CLOCK_TIME = 20     # frame rate period in ms -> corresponds to (1 / ms * 1/1000) -> frequency in Hz)
SCAL       = 0.10   # scale event speed - first game loop
SCAL2      = 1.00   # scale event speed - second game loop

"""
=============
Timing notes:
=============

@ CLOCK_TIME -> 20, SCAL -> 0.08, SCAL2 -> 1.00
-----------------------------------------------

First frame: 0.021381855010986328 s
* the time when protein first appears on window
* Say "Rusticyanin obtained from PDB"

Remove residues: 3.7628419399261475 s
* the time at which residues are removed
* Say "Only banking methionine and the aromatics"

Remove atoms: 6.114755868911743 s
* the time at which unnecessary atoms are removed
* Say "Strip off unnecessary atoms"

Web with v: 8.206990718841553 s
* the time when the entire protein webbed with vectors v
* Say "Project vectors from methionine sulfur to aromatics"

Bank closely spaced residues: 10.735826015472412 s
* the time when closely spaced residues are removed
* Say "Bank only closely spaced methionine aromatic pairs"

Zoom into SD: 14.840545892715454 s
* here I zoom into the SD
* I also put in an eigenvector
Say "Zoom in a bit"

Start rotation of LPs: 21.787418842315674 s
* start of rotation of lone pairs
Say "Lastly rotate the lone pairs into their natural position"


@ CLOCK_TIME -> 20, SCAL -> 0.10, SCAL2 -> 1.00
-----------------------------------------------
First frame: 0.021639108657836914 s
Remove residues: 5.1682679653167725 s
Remove atoms: 8.098443031311035 s
Web with v: 10.719419002532959 s
Bank closely spaced residues: 13.92305874824524 s
Zoom into SD: 18.514067888259888 s
Start rotation of LPs: 25.40074586868286 s
"""

# a dictionary of atom colors
DICT_COLORS = {
    'N': (0.0, 0.0, 1.0),
    'O': (1.0, 0.0, 0.0),
    'C': (0.75, 0.75, 0.75),
    'S': (1.0, 1.0, 0)
}

# a dictionary of atom radii
# carbon radius normalized to 0.01
DICT_RADII = {
    'N': 0.009,
    'O': 0.008,
    'C': 0.010,
    'S': 0.0145
}

def drawVector(start, end, width=1.0):
    # a function for drawing a line from one point to another
    # does not include an arrowhead
    glLineWidth(width)
    glBegin(GL_LINES)
    glVertex3f(*start)
    glVertex3f(*end)
    glEnd()

def drawAtoms(data):
    # populate the current matrix stack
    # stack is populated with customized spheres for atom identities
    glPushMatrix()
    glTranslatef(*data[5]) 
    glColor3f(*DICT_COLORS.get(data[4]))  # color different atoms
    glutSolidSphere(DICT_RADII.get(data[4]), *SLICES_STACKS)  # custom radii
    glPopMatrix()

DICT_ATOMS_PHE = {
'CG':'A', 'CD2':'B', 'CE2':'C', 
'CZ':'D', 'CE1':'E', 'CD1':'F'
}

DICT_ATOMS_TYR = {
'CG':'A', 'CD2':'B', 'CE2':'C', 
'CZ':'D', 'CE1':'E', 'CD1':'F'
}

DICT_ATOMS_TRP = {
'CD2':'A', 'CE3':'B', 'CZ3':'C', 
'CH2':'D', 'CZ2':'E', 'CE2':'F'
}

# first game loop
# ---------------
INC_1  =           (0,  SCAL * 1000)               # do nothing
INC_2  = (SCAL * 1000,  SCAL * 2000)               # rotate the entire molecule for a second
INC_3  = (SCAL * 2000,  SCAL * 3000)               # pause
INC_4  = (SCAL * 3000,  SCAL * 4000)               # remove unneeded residues
INC_5  = (SCAL * 4000,  SCAL * 5000)               # remove unneeded atoms

# second game loop
# ----------------
INC_6  = (0,            100)                       # zoom in is fixed at 100 frames
INC_7  = (100,          200)                       # rotation rate can also be held constant
INC_8  = (300,  SCAL2 * 400)                       # draw in antiparallel CG-SD / CE-SD vectors
INC_9  = (INC_8[1], INC_8[1] + 90)                 # fixed 90 frame rotation

# --------------------------------------------------
# preprocessing prior to rendering

# get 1rcy data
RAW_INPUT_A = import_coords()  

# append in key data, coordinate data and the Euclidean norm of coordinate data
# new data of form: ['N', 'THR', 'A', '5', 'N', array([-13.081,   4.669,  24.745]), 28.37652457578271]
RAW_INPUT_B = []
for data in RAW_INPUT_A:
    coord = array(data[6:9]).astype(float)
    RAW_INPUT_B.append([*data[2:6], data[11], coord, linalg.norm(coord)])

# get the most "far out" coordinate such that we can normalize vectors to frustrum
N_C = 1 / max(i[6] for i in RAW_INPUT_B)
VEC_V_NORM_NORMALIZED = CUTOFF * N_C

# normalize data
for i in range(0, len(RAW_INPUT_B)):
    RAW_INPUT_B[i][5] = N_C * RAW_INPUT_B[i][5]

# drop Euclidean norm data - no need to overload memory
ALL_ATOMS = [i[0:6] for i in RAW_INPUT_B]

# strip down to residues of interest - need all atoms for animation purposes
STRIPPED_TO_RES = [i for i in ALL_ATOMS if search(RES, i[1]) != None]

# strip down to specific atoms using regex
# need separated data for Rodrigues rotation, midpoints and pooled for animation
DATA_MET = [i for i in STRIPPED_TO_RES if i[1] == 'MET' and search(ATOMS_MET, i[0]) != None]
DATA_PHE = [i for i in STRIPPED_TO_RES if i[1] == 'PHE' and search(ATOMS_PHE, i[0]) != None]
DATA_TYR = [i for i in STRIPPED_TO_RES if i[1] == 'TYR' and search(ATOMS_TYR, i[0]) != None]
DATA_TRP = [i for i in STRIPPED_TO_RES if i[1] == 'TRP' and search(ATOMS_TRP, i[0]) != None]
STRIPPED_TO_ATOMS = DATA_MET + DATA_PHE + DATA_TYR + DATA_TRP

# precompute midpoints
# sort data prior to applying groupby operations
DATA_ARO = DATA_PHE + DATA_TYR + DATA_TRP
DATA_ARO = sorted(DATA_ARO, key=itemgetter(3))  # note lexicographic ordering

# apply groupby operator
DATA_ARO = [list(group) for _, group in groupby(DATA_ARO, lambda x: x[3])]  

# get midpoints
MIDPOINTS = []
for grouped in DATA_ARO:
    # map unique values to atomic label keys
    for row in grouped:
        if row[1] == 'PHE':
            row[0] = DICT_ATOMS_PHE.get(row[0])
        elif row[1] == 'TYR':
            row[0] = DICT_ATOMS_TYR.get(row[0])
        else:
            row[0] = DICT_ATOMS_TRP.get(row[0])
    
    # then sort based on these values which are just A, B, C, D, E, F
    ordered = sorted(grouped, key=itemgetter(0))

    # rotate data
    xyz_A = [c[5] for c in ordered]
    deque_xyz_A = deque(xyz_A)
    deque_xyz_A.rotate(1)
    xyz_B = list(deque_xyz_A)

    # get midpoints
    for A, B in zip(xyz_A, xyz_B):
        MIDPOINTS.append([ordered[0][3], 0.5 * (A + B)])
    
# get all vectors v
SD = [i for i in DATA_MET if i[0] == 'SD']
VECS_V = []
for sds in SD:
    for mps in MIDPOINTS:
        VECS_V.append([sds[5], mps[1], sds[3], mps[0]])

# apply distance condition to vectors v
VECS_V_STRIPPED  = [i for i in VECS_V if linalg.norm(i[1] - i[0]) <= VEC_V_NORM_NORMALIZED]
LABELS           = list(set(chain(*[i[2:4] for i in VECS_V_STRIPPED])))
ATOMS_V_STRIPPED = [i for i in STRIPPED_TO_ATOMS if i[3] in LABELS]

# hard code for zoom in, eigenvector about CG - SD - CE
Z_TRANS, INC = 0.00, 0.05
CG = ATOMS_V_STRIPPED[0][5]  # rename for clarity
SD = ATOMS_V_STRIPPED[1][5]
CE = ATOMS_V_STRIPPED[2][5]
MP_CG_CE = 0.5 * (CG + CE)
EIG_VEC = SD - MP_CG_CE      # the eigenvector of a tetrahedral rotation

# get CE-SD, CG-SD antiparallel vectors mapped to origin
AP_CG_SD = -1 * (CG - SD)
AP_CE_SD = -1 * (CE - SD)

# --------------------------------------------------
# game loop

pygame.init() 
pygame.display.set_caption('MA-Animation')
display = (DISPLAY_DIMENSION, DISPLAY_DIMENSION)
pygame.display.set_mode(display, DOUBLEBUF|OPENGL)

# setup the frustrum (viewing space)
angle_fov = 45.0                                        # field of view angle
aspect_ratio = display[0] / display[1]                  # aspect ratio
z_start = 0.1                                           # the "start" of the viewing plane
z_end = 50.0                                            # the "end" of the viewing plane
gluPerspective(angle_fov, aspect_ratio, z_start, z_end)

"""
# starting camera position
position_camera = (0.00, 0.00,  3.00)
position_object = (0.00, 0.00,  0.00)
rotation_camera = (0.00, 1.00,  0.00)
gluLookAt(*position_camera, *position_object, *rotation_camera)
"""
glTranslatef(0.00, 0.00, -3.00)  # I'm using glTranslatef() for simplicity
                                 # gluLookAt() wraps glTranslatef(), glRotatef() code anyways


# timing
# ------
t_start = time()
TICK_1 = True
TICK_2 = True
TICK_3 = True
TICK_4 = True
TICK_5 = True
TICK_6 = True
TICK_7 = True


FRAMES = 0
while True:
    pygame.time.wait(CLOCK_TIME) # set frame rate
    FRAMES += 1                  # our counter is the number of frames that have been rendered
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)


    # poll events
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            quit()


    # render raw protein structure
    # ----------------------------
    
    if FRAMES >= INC_1[0]   and FRAMES < INC_1[1]: 
        if TICK_1: 
            print('First frame: {} s'.format(time() - t_start))
            TICK_1 = False
        for data in ALL_ATOMS:
            drawAtoms(data)


    # remove unneeded residues
    # ------------------------
    elif FRAMES >= INC_2[0] and FRAMES < INC_2[1]:
        if TICK_2: 
            print('Remove residues: {} s'.format(time() - t_start))
            TICK_2 = False
        for data in STRIPPED_TO_RES:
            drawAtoms(data)


    # remove unneeded atoms
    # ---------------------
    elif FRAMES >= INC_3[0] and FRAMES < INC_3[1]:
        if TICK_3: 
            print('Remove atoms: {} s'.format(time() - t_start))
            TICK_3 = False
        for data in STRIPPED_TO_ATOMS:
            drawAtoms(data)


    # web the feature space with vectors v
    # ------------------------------------
    elif FRAMES >= INC_4[0] and FRAMES < INC_4[1]:
        if TICK_4: 
            print('Web with v: {} s'.format(time() - t_start))
            TICK_4 = False
        for data in STRIPPED_TO_ATOMS:
            drawAtoms(data)
        for pairs in VECS_V:  # draw in vectors v
            drawVector(*pairs[0:2])


    # strip atoms + vectors not meeting distance condition
    # ----------------------------------------------------
    elif FRAMES >= INC_5[0] and FRAMES < INC_5[1]:
        if TICK_5: 
            print('Bank closely spaced residues: {} s'.format(time() - t_start))
            TICK_5 = False
        for data in ATOMS_V_STRIPPED:
            drawAtoms(data)
        for pairs in VECS_V_STRIPPED:
            drawVector(*pairs[0:2])


    # GOTO next loop
    # where I zoom in to a single met aromatic interaction
    else:
        break

    glFlush()
    glutSwapBuffers()
    glutPostRedisplay()
    pygame.display.flip()
    

# scaled SD x, y coords by 1/100 () and translated towards origin
# we move from (-0.12433176, -0.26009507, -3.00) to (0.00, 0.00, -1.00)
TRANSLATE_TO_XY_ORIGIN = (1 / INC_6[1] * -SD[0], 1 / INC_6[1] * -SD[1], 0.02)

FRAMES = 0 # reset frame count to zero
while True:
    pygame.time.wait(CLOCK_TIME)
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
    FRAMES += 1

    # poll events
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            quit()


    # here we create a new origin - the SD coordinate
    # -----------------------------------------------
    if FRAMES < INC_6[1]:
        glTranslatef(*TRANSLATE_TO_XY_ORIGIN)

    
    # rotate origin for a better view
    # -------------------------------
    if FRAMES >= INC_7[0] and FRAMES < INC_7[1]:
        if TICK_6: 
            print('Zoom into SD: {} s'.format(time() - t_start))
            TICK_6 = False
        glTranslatef(*SD)                                          # map to origin
        glRotatef( 90 / (INC_7[1] - INC_7[0]), 1.00, 0.00, 0.00)   # rotate about x axis, do mult in place to preserve precision
        glTranslatef(*-SD)                                         # map back to original coordinates


    # draw in eigenvector
    # fixed pause of 100 frames
    # -------------------
    if FRAMES >= INC_7[1]:
        glColor3f(0.00, 0.00, 1.00)                      # render blue eigenvector
        drawVector(-4 * EIG_VEC + SD, 4 * EIG_VEC + SD)  # here we also extend eigenvector in both directions
    

    # draw in antiparallel CG-SD / CE-SD vectors
    # need custom pause here to describe these vectors
    # ------------------------------------------
    if FRAMES >= INC_8[0] and FRAMES < INC_8[1]:
        glColor3f(0.00, 0.00, 1.00)
        drawVector(SD, AP_CG_SD + SD, width=LP_WIDTHS)
        drawVector(SD, AP_CE_SD + SD, width=LP_WIDTHS)

    
    # do an Euler rotation on antiparallel vectors
    # here we strictly iterate over 90 frames
    # --------------------------------------------
    if FRAMES >= INC_9[0] and FRAMES <= INC_9[1]:
        if TICK_7: 
            print('Start rotation of LPs: {} s'.format(time() - t_start))
            TICK_7 = False
        glPushMatrix()
        glTranslatef(*SD)
        glRotatef((FRAMES - INC_8[1]), *EIG_VEC)
        glTranslatef(*-SD)
        drawVector(SD, AP_CG_SD + SD, width=LP_WIDTHS)
        drawVector(SD, AP_CE_SD + SD, width=LP_WIDTHS)
        glPopMatrix()


    # hold the Euler rotation stationary at 90.00 for 100 frames
    # ----------------------------------------------------------
    if FRAMES > INC_9[1]:
        glPushMatrix()
        glTranslatef(*SD)
        glRotatef(90.00, *EIG_VEC)
        glTranslatef(*-SD)
        drawVector(SD, AP_CG_SD + SD, width=LP_WIDTHS)
        drawVector(SD, AP_CE_SD + SD, width=LP_WIDTHS)
        glPopMatrix()


    # rotate about x axis to better see lone pair positions
    # -----------------------------------------------------
    if FRAMES >= INC_9[1] + 100 and FRAMES < INC_9[1] + 146: # 146 <- rotate by 46 deg, (46 frames)
        glTranslatef(*SD)
        glRotatef(1.00, 1.00, 0.00, 0.00) # note we don't increment by FRAMES - INC_9... here
        glTranslatef(*-SD)                # -> because we're not doing Euler rotation


    for data in ATOMS_V_STRIPPED:
        drawAtoms(data)

    for pairs in VECS_V_STRIPPED:
        drawVector(*pairs[0:2])

    glFlush()
    glutSwapBuffers()
    glutPostRedisplay()
    pygame.display.flip()
