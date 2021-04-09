#from os import F_LOCK
from parameters import *
from boid import *
from attractor import *
from slider import *
import time



def quit_game():
    pygame.quit()
    quit()

def handle_inputs(dt, occ_list):
    global move_up, move_down, move_left, move_right
    
    if pygame.mouse.get_pressed()[0] == True:
            reassign_attractors(occ_list)

    # checks and processes User Inputs
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            quit_game()

        if event.type == pygame.KEYDOWN:

            if event.key == pygame.K_UP:
                move_up = True
            if event.key == pygame.K_DOWN:
                move_down = True
            if event.key == pygame.K_LEFT:
                move_left = True
            if event.key == pygame.K_RIGHT:
                move_right = True


            if event.key == pygame.K_ESCAPE:
                quit_game()


        if event.type == pygame.KEYUP:
            if event.key == pygame.K_UP:
                move_up = False
            if event.key == pygame.K_DOWN:
                move_down = False
            if event.key == pygame.K_LEFT:
                move_left = False
            if event.key == pygame.K_RIGHT:
                move_right = False

    # perform movements
    if move_up:
        update_background(pygame.Vector2(0,-attractor_speed),dt)
    if move_down:
        update_background(pygame.Vector2(0,attractor_speed),dt)
    if move_left:
        update_background(pygame.Vector2(-attractor_speed,0),dt)
    if move_right:
        update_background(pygame.Vector2(attractor_speed,0),dt)

##############################
#         BACKGROUND         #
##############################

def init_background(occ_list):

    # make QD
    Attractor(pygame.Vector2(display_width / 2, display_height / 4))

    # make bath sites
    NBATH = (len(occ_list) - 2) // 2
    for i in range(NBATH):
        Attractor(pygame.Vector2((i+1)*display_width/(NBATH+1), 3*display_height/4))

def draw_background():
    screen.fill(screenColor)
    for attractor in Attractor.instances:
        attractor.draw()

def update_background(velocity, dt):
    for attractor in Attractor.instances:
        attractor.update(velocity, dt)

##############################
#           Boids            #
##############################

def boid_number_per_site(occ_list):
    # the following code assigns attractors to boids according to distribution (as close as possible) given by occ_list
    # it does so in such a way to preserve the given total_boids number, and spin

    total_weights = sum(occ_list)
    NSITES = len(occ_list) // 2

    # "ideal boids number per spin-x site"
    iUp = [occ_list[i] * total_boids/total_weights for i in range(NSITES)]
    iDown = [occ_list[i+NSITES] * total_boids/total_weights  for i in range(NSITES)]

    # "floor of boid number per spin-x site"
    fUp = [int(i) for i in iUp]
    fDown = [int(i) for i in iDown]

    # "remainder boid number per spin-x site"
    rUp = [iUp[i] - fUp[i] for i in range(NSITES)]
    rDown = [iUp[i] - fUp[i] for i in range(NSITES)]

    # IF SPINS SHOULD BE INDEPENDENT THEN CHANGE THE FF CODE
    # (BECAUSE OF THE FF CODE ONLY CASES WITH EQUAL UP AND DOWN SPIN OCCUPATIONS WILL WORK -> need even total_boids)
    # assign boids leftover due to rounding down
    leftoverUp = total_boids // 2 - sum(fUp)
    leftoverDown = total_boids // 2 - sum(fDown)

    for i in range(leftoverUp):
        # find which site is missing the most boids
        shortest = rUp.index(max(rUp))
        fUp[shortest] += 1
        rUp[shortest] = 0

    for i in range(leftoverDown):
        shortest = rDown.index(max(rDown))
        fDown[shortest] += 1
        rDown[shortest] = 0

    return fUp, fDown

def init_boids(occ_list):
    
    fUp, fDown = boid_number_per_site(occ_list)

    # actually initialise boids
    for i, n in enumerate(fUp):
        for _ in range(n):
            Boid(attractorid=i, color=spinUpColor)
    
    for i, n in enumerate(fDown):
        for _ in range(n):
            Boid(attractorid=i, color=spinDownColor)

def reassign_attractors(occ_list):
    
    NSITES = len(occ_list) // 2
    # get new boid numbers per site
    up, down = boid_number_per_site(occ_list)

    # get current number of Up and Down boids per site
    currentUp = [0]*NSITES
    currentDown = [0]*NSITES
    for boid in Boid.instances:
        if boid.color == spinUpColor:
            currentUp[boid.attractorid] += 1
        elif boid.color == spinDownColor:
            currentDown[boid.attractorid] += 1

    # generate list showing how many boids each site needs to gain (or loose if negative)
    needUp = [up[i] - currentUp[i] for i in range(NSITES)]
    needDown = [down[i] - currentDown[i] for i in range(NSITES)]

    # reassign boids until need lists are empty. Shuffle boids for maybe more organic movement
    # random.shuffle(boids)
    for boid in Boid.instances:
        if boid.color == spinUpColor:
            # check if boid belongs to overpopulated attractor
            if needUp[boid.attractorid] < 0:
                # remove 1 boid from that attractors needs to give away
                needUp[boid.attractorid] +=1 

                # find new attractor for boid (first index with positive value in needUp)
                newAttractorId = [i for i, x in enumerate(needUp) if x > 0][0]
                
                # reassign boid
                boid.attractorid = newAttractorId

                # remove 1 boid to new attractors need to accept boids
                needUp[newAttractorId] -= 1
        
        # same for spin-down
        if boid.color == spinDownColor:
            if needDown[boid.attractorid] < 0:

                needDown[boid.attractorid] +=1 

                newAttractorId = [i for i, x in enumerate(needDown) if x > 0][0]
                boid.attractorid = newAttractorId

                needDown[newAttractorId] -= 1

    
def update_boids(dt):
    random.shuffle(Boid.instances)
    for boid in Boid.instances:
        boid.edges()
        boid.update(Attractor.instances, dt)
        boid.draw()

##############################
#          Sliders           #
##############################

def init_sliders(occ_list):
    global Vg, X

    Vg = Slider("Vg", sx=70, sy=50, width=200, height=30, initialValue=0, valueMin=0, valueMax=10)
    X = Slider("X", sx=70, sy=130, width=200, height=30, initialValue=0, valueMin=0, valueMax=10)

    update_sliders(occ_list)
    reassign_attractors(occ_list)

def update_sliders(occ_list):
    occ_list[0] = Vg.value
    occ_list[5] = Vg.value
    occ_list[1] = X.value
    occ_list[6] = X.value

    Slider.updateAll()
    Slider.drawAll()

def game_loop(occ_list) -> None:

    init_background(occ_list)
    init_boids(occ_list)
    init_sliders(occ_list)

    while True:

        dt = clock.tick(60) *.001 * TARGET_FPS

        handle_inputs(dt, occ_list)       
        draw_background()
        update_sliders(occ_list)
        update_boids(dt)
    
        pygame.display.flip()

game_loop(occ_list=np.array([0.2, 1, 1, 2, 2, 0.2, 1, 1, 2, 2]))
pygame.quit()
quit()