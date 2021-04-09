from parameters import *

class Slider:
    instances = []

    def __init__(self, name:str, sx:int, sy:int, width:int, height:int, initialValue:float, valueMin:float, valueMax:float) -> None:
        
        # adds itself to list of all slider instances
        self.__class__.instances.append(self)

        # argument read in
        self.name = name
        self.sx = sx
        self.sy = sy
        self.width = width
        self.height = height

        # button properties
        self.buttonColor = (255, 255, 255)
        self.buttonBorderThickness = 1
        self.buttonWidth = 25

        # slider properties
        self.sliderBarColor = (220, 220, 220)
        self.sliderBarBorderThickness = 0
        self.sliderBarHeight = 3
        self.sliderEndPadding = 30
        self.sliderHeightPadding = 20

        # slider steps
        self.sliderStepColor = (30, 30, 30)
        self.sliderStepBorderThickness = 0
        self.sliderStepWidth = 2
        self.sliderStepHeight = 20

        # slider values
        self.value = initialValue
        self.valueMax = valueMax
        self.valueMin = valueMin

        # slider coordinate [0,1]
        self.pos = 0
        self.set_pos()

    def draw(self) -> None:
        
        # draw slider bar
        pygame.draw.rect(screen,self.sliderBarColor,(self.sx,self.sy+self.height/2-self.sliderBarHeight/2,self.width,self.sliderBarHeight)) # draw line
        
        # draw slider steps
        for i in range(0,11):
            pygame.draw.rect(screen,self.sliderStepColor,
                            (round(self.sx+i*10*((self.width)/100)),
                            self.sy+(self.height-self.sliderStepHeight)/2,
                            self.sliderStepWidth,
                            self.sliderStepHeight))

            if self.sliderStepBorderThickness != 0:
                pygame.draw.rect(screen,(0,0,0),
                                (round(self.sx+i*10*((self.width)/100)),
                                self.sy+(self.height-self.sliderStepHeight)/2,
                                self.sliderStepWidth,self.sliderStepHeight),
                                self.sliderStepBorderThickness)
        
        # draw slider button
        Bx = self.sx+self.pos*self.width-self.buttonWidth/2 + self.sliderStepWidth/2
        By = self.sy
        pygame.draw.rect(screen,self.buttonColor,(Bx, By,self.buttonWidth,self.height)) # draw button fill

        # draw value above slider button
        value_surface, value_rect = FONT1.render( "{:.2f}".format(self.value), black)
        screen.blit(value_surface, (Bx - value_rect[2]/2 + self.buttonWidth/2, By - value_rect[3]-10))


        # draw button & slider borders
        if self.buttonBorderThickness != 0:
            pygame.draw.rect(screen,(0,0,0), (Bx, By,self.buttonWidth,self.height), self.buttonBorderThickness)
        if self.sliderBarBorderThickness != 0:
            pygame.draw.rect(screen,(0,0,0),
                            (self.sx,self.sy+self.height/2-self.sliderBarHeight/2,self.width,self.sliderBarHeight),
                            self.sliderBarBorderThickness)
        
        # draw name next to slider
        name_surface, name_rect = FONT1.render(self.name, black)
        screen.blit(name_surface, (self.sx - name_rect[2] - self.buttonWidth - 5, self.sy + self.height/2 - name_rect[3]/2))

        pygame.display.update(pygame.Rect(self.sx,self.sy,self.width,self.height))
    
    def update(self):
        mouse = pygame.mouse.get_pos()
        click = pygame.mouse.get_pressed()
        
        # check if mouse on slider and clicked
        if mouse[0] > self.sx and mouse[0] < self.sx+self.width and mouse[1] > self.sy - self.sliderHeightPadding and mouse[1] < self.sy+self.height + self.sliderHeightPadding:
            if click[0] == 1:
                
                mouseX = mouse[0]
                self.pos = (mouseX- self.sx) / abs(self.sx + self.width - self.sx)
                self.set_value()
        
        # extra region for mouse on extreme slider ends
        if mouse[0] > self.sx - self.sliderEndPadding and mouse[0] < self.sx+self.width + self.sliderEndPadding and mouse[1] > self.sy and mouse[1] < self.sy+self.height:
            if click[0] == 1:
                if mouse[0] < self.sx:
                    self.pos = 0
                    self.set_value()
                elif mouse[0] > self.sx+self.width:
                    self.pos = 1
                    self.set_value()

    @classmethod
    def drawAll(cls) -> None:
        for instance in cls.instances:
            instance.draw()

    @classmethod
    def updateAll(cls) -> None:
        for instance in cls.instances:
            instance.update()

    def set_value(self) -> None:
        self.value = self.pos * abs(self.valueMax - self.valueMin) + self.valueMin

    def set_pos(self) -> None:
        self.pos = (self.value- self.valueMin) / abs(self.valueMax - self.valueMin)
