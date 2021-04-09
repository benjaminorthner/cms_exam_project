import pygame
import pygame.gfxdraw # for anti aliased drawing
import pygame.freetype # for rendering text
import numpy as np
import random
import gc # garbage collector for finding sliders

# removes ALSA LIB error
import os
os.environ['SDL_AUDIODRIVER'] = 'dsp'

pygame.init()

display_width = 800
display_height = 600

black = (0, 0, 0)
screenColor = (200, 200, 200)
white = (255, 255, 255)
red = (255, 0, 0)
blue = (30,30,255)
spinUpColor = blue
spinDownColor = red

screen = pygame.display.set_mode((display_width, display_height))
pygame.display.set_caption('QD Simulator')
FONT1 = pygame.freetype.SysFont("Helvetica", 17)

clock = pygame.time.Clock()
TARGET_FPS = 60

attractor_speed = 3
max_boid_speed = 3

# will be rounded down to nearest whole number
total_boids = 400

# init global variables for attractor movement
move_up, move_down, move_left, move_right = False, False, False, False

screen.fill(screenColor)
pygame.display.flip()