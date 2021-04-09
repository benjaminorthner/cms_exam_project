from parameters import *

class Attractor:
    
    instances = []

    def __init__(self, position) -> None:

        # adds itself to list of all attractor instances
        self.__class__.instances.append(self)
        # list of all boids assigned to this attractor instance
        # not in use yet
        self.myboids = [] 

        # set initial position and velocity
        self.position = position
        self.velocity = pygame.Vector2(0,0)

    def update(self, velocity, dt) -> None:
        self.velocity = velocity
        self.position += velocity * dt

    def draw(self) -> None:
        pygame.gfxdraw.filled_circle(screen, int(self.position[0]), int(self.position[1]), 10, black)
        pygame.gfxdraw.aacircle(screen, int(self.position[0]), int(self.position[1]), 10, black)