from parameters import *

# for extra efficiency create flock bins (maye by attractorid)
# and combine cohesion alignment and separation into one function that only loops over all other boids once
class Boid:

    instances = []

    def __init__(self, attractorid=0, color=(255, 0, 0)) -> None:

        # adds itself to list of all boid instances
        self.__class__.instances.append(self)

        # inital position
        self.position = pygame.Vector2(display_width /2, display_height/2)

        # initial velocity
        self.velocity = pygame.Vector2(0.5 - random.random(), 0.5 - random.random())*2
        self.velocity = self.velocity / random.random()
        #self.velocity = pygame.Vector2(0,0)

        # initial acceleration
        self.acceleration = pygame.Vector2(0,0)

        # attractorid and color
        self.attractorid = attractorid
        self.color = color

    def update(self, attractors, dt) -> None:

        # collect forces
        self.acceleration += self.attraction(attractors)
        #self.acceleration += self.alignment()
        #self.acceleration += self.cohesion()
        #self.acceleration += self.separation()
        self.acceleration += self.thermal()

        # actually update
        self.velocity += self.acceleration
        speed = self.velocity.magnitude()
        if speed > max_boid_speed:
           self.velocity = max_boid_speed * self.velocity / speed

        self.position += self.velocity * dt
        
        # set acceleration back to 0
        self.acceleration = self.acceleration*0

    #temp solution only works for exactly 1 attractor
    def attraction(self, attractors: list) -> float:
        attraction_strength = 0.01
        return attraction_strength*(attractors[self.attractorid].position - self.position)

    def alignment(self):
        perceptionRadius = 60
        alignment_strength = 1
        steering = pygame.math.Vector2(0,0)

        for boid in self.__class__.instances:
            diff = self.position - boid.position
            d = diff.magnitude()
            boid_count = 0
            if boid != self and 0 < d <= perceptionRadius and boid.attractorid == self.attractorid:
                steering += boid.velocity
                boid_count += 1
        
        if boid_count >0:
            return alignment_strength*steering / boid_count
        else:
            return pygame.math.Vector2(0,0)
        
    def cohesion(self):
        perceptionRadius = 100
        cohesion_strength=0.00001
        COM = pygame.math.Vector2(0,0)
        boid_count = 0

        for boid in self.__class__.instances:
            d = self.position.distance_to(boid.position)
            if boid != self and 0 < d <= perceptionRadius and boid.attractorid == self.attractorid:
                COM += boid.position
                boid_count += 1
        
        if boid_count > 0:
            return cohesion_strength*(COM-self.position)
        else:
            return pygame.math.Vector2(0,0)


    def separation(self):
        perceptionRadius = 30
        separation_strength = 0.1
        steering = pygame.math.Vector2(0,0)

        for boid in self.__class__.instances:
            diff = self.position - boid.position
            d = diff.magnitude()
            if boid != self and 0 < d <= perceptionRadius and boid.attractorid == self.attractorid:
                steering += separation_strength*(1.0/d)*diff
            else: # small random jitter to induce initial motion
                steering += pygame.math.Vector2(0.5-random.random(),0.5-random.random())*0.01
            
        return steering

    def thermal(self):
        thermal_strength = max_boid_speed / 2
        return pygame.math.Vector2(0.5-random.random(),0.5-random.random())*thermal_strength


    def edges(self) -> None:
        if self.position[0] > display_width:
            self.position[0] = 0
        elif self.position[0] < 0:
            self.position[0] = display_width

        if self.position[1] > display_height:
            self.position[1] = 0
        elif self.position[1] < 0:
            self.position[1] = display_height

    def draw(self) -> None:
        pygame.gfxdraw.filled_circle(screen, int(self.position.x), int(self.position.y), 3, self.color)