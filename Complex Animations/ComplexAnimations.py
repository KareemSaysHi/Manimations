from manim import *
import numpy as np

class ComplexStuff(Scene):
    def construct(self):
        complexplane = ComplexPlane(x_range = (-25, 25, .25), y_range = (-25, 25, .25))
        movingcomplexplane = complexplane.copy()
        movingcomplexplane.prepare_for_nonlinear_transform()
        complexplane.set_stroke(BLUE_E, 1)
        
        self.play(
            Write(complexplane, run_time = 3),
            FadeIn(movingcomplexplane),
        )
        self.wait(1)
        self.play(
            movingcomplexplane.animate.apply_complex_function(lambda z: np.exp(z)), #change this function to animate whatever you want
            run_time = 5
        )
        self.wait(2)
        