from manim import *
import numpy as np
import cmath

hbar = 2 #constants
eme =  1 #more constants
points = []

class TimeIndependentHamiltonian(ThreeDScene):
    
    def sphereFunc(self, u, v):
        return np.array([np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v)])

    def hamilFunc(self, t, eplus, eminus, ogpsiquantum, evectorplus, evectorminus):
        thing = self.quantumToSpherical((cmath.exp(complex(0, -1*eplus*t/hbar))*ogpsiquantum[0][0])*evectorplus + (cmath.exp(complex(0, -1*eminus*t/hbar))*ogpsiquantum[1][0])*evectorminus)
        #print ((ogpsiquantum[0][0])*evectorplus)
        #print ("quantum state: " + str(cmath.exp(complex(0, -1*eplus*t/hbar))*ogpsiquantum[0][0])*evectorplus + (cmath.exp(complex(0, -1*eminus*t/hbar))*ogpsiquantum[1][0])*evectorminus)
        #print("spherical: " + str(thing))
        #print ("position: " + str(self.sphericalToRectangular(thing[0], thing[1], thing[2])))
        #print ("")
        return self.sphericalToRectangular(thing[0].real, thing[1].real, thing[2].real)

    def get_param_func(self, axes, t, eplus, eminus, ogpsiquantum, evectorplus, evectorminus):
        tol = 1e-9

        pf = ParametricFunction(
            lambda t: axes.c2p(*self.hamilFunc(t, eplus, eminus, ogpsiquantum, evectorplus, evectorminus)),
            t_range=[0, t+.02, 0.005],
            tolerance_for_point_equality=tol
        )
        
        pf.set_color(color=[RED, YELLOW, BLUE, RED])
        return pf



    def quantumToSpherical(self, quantum):
        #print ("modulus = " + str(abs(quantum[0][0])+abs(quantum[1][0])))
        #print(quantum[0][0])
        #print("old: " + str(quantum))
        overallphase = cmath.phase(quantum[0][0])
        newquantum = np.array([[quantum[0][0] * cmath.exp(complex(0, -1*overallphase))], [quantum[1][0]* cmath.exp(complex(0, -1*overallphase))]], dtype= "complex_")
        #print(newquantum)
        #print("")
        foundtheta = np.arccos(newquantum[0][0])*2
        #print (foundtheta)
        if abs(foundtheta) < 0.000001:
            return [1, 0, 0] #spin up in z basis
        if abs(foundtheta) > PI-0.0000001:
            return [1, PI, 0] #spin down in z basis
        thinginside = (newquantum[1][0])/np.sin(foundtheta/2)
        #print(quantum)
        #print(thinginside)
        foundphi = cmath.log(thinginside)*complex(0, -1)
        return np.array([1, foundtheta, foundphi])

    def sphericalToRectangular(self, rho, theta, phi):
        x = rho*np.sin(theta)*np.cos(phi)
        y = rho*np.sin(theta)*np.sin(phi)
        z = rho*np.cos(theta)
        return np.array([x, y, z])

    def construct(self):       
        axes = ThreeDAxes(x_range=(-2, 2, 1), y_range=(-2, 2, 1), z_range=(-2, 2, 1), x_length=8, y_length=8, z_length=8) #defining axes
        sphere = Surface(
            lambda u, v: axes.c2p(*self.sphereFunc(u, v)), #make a sphere
            u_range=[-PI, PI],
            v_range=[0, TAU],
            fill_opacity=.1
        )
       
        timer = ValueTracker(0) #timer value

        brho = 1
        btheta = PI/4
        bphi = 0
        brectangular = self.sphericalToRectangular(brho, btheta, bphi) #note, might have to add updater
        eplus = eme*hbar/2*brho #defining eigenenergies, these don't change with basis
        eminus = -1*eplus 
        evectorplus = np.array([[np.cos(btheta/2)], [np.sin(btheta/2) * cmath.exp(complex(0, bphi))]]) #defining eigenvectors (written in z basis)
        evectorminus = np.array([[np.sin(btheta/2)], [-np.cos(btheta/2) * cmath.exp(complex(0, bphi))]]) 


        psitheta = PI/2
        psiphi = PI/6
        #ogpsiquantum = np.array([[np.cos(btheta/2)], [np.sin(btheta/2)]])

        ogpsiquantum = np.array([[np.cos(btheta/2)*np.cos(psitheta/2)+np.sin(btheta/2)*np.sin(psitheta/2)*cmath.exp(complex(0, psiphi-bphi))], [np.sin(btheta/2)*np.cos(psitheta/2)-np.cos(btheta/2)*np.sin(psitheta/2)*cmath.exp(complex(0, psiphi-bphi))]])

        ogpsirectangular = self.sphericalToRectangular(1.0, psitheta, psiphi) #note, might have to add updater
        print(ogpsirectangular)
        print(brectangular)

        #changingpsiquantum = ogpsiquantum.add_updater(lambda m: m.become(
        #    cmath.exp(-i*eplus*timer.get_value()/hbar)*ogpsiquantum[0]*evectorplus + cmath.exp(-i*eminus*timer.get_value()/hbar)*ogpsiquantum[1]*evectorminus
        #))

        #psirectangular = self.sphericalToRectangular(self.quantumToSpherical(ogpsiquantum)).add_updater(lambda n: n.become(self.sphericalToRectangular(self.quantumToSpherical(changingpsiquantum))))




        bvector = Arrow((0, 0, 0), axes.c2p(brectangular[0], brectangular[1], brectangular[2]), buff=0, color=BLUE_D) #making b-vector out of 

        psivector = Arrow((0, 0, 0), axes.c2p(ogpsirectangular[0], ogpsirectangular[1], ogpsirectangular[2]), buff=0, color=WHITE).add_updater(lambda m: m.put_start_and_end_on([0,0,0], axes.c2p(self.hamilFunc(timer.get_value(), eplus, eminus, ogpsiquantum, evectorplus, evectorminus)[0], self.hamilFunc(timer.get_value(), eplus, eminus, ogpsiquantum, evectorplus, evectorminus)[1], self.hamilFunc(timer.get_value(), eplus, eminus, ogpsiquantum, evectorplus, evectorminus)[2])))

        #psivector = Arrow((0, 0, 0), axes.c2p(0, 0, 1), buff=0, color=WHITE)

        pf = self.get_param_func(axes, timer.get_value(), eplus, eminus, ogpsiquantum, evectorplus, evectorminus)   
        pf.add_updater(
            lambda mob: mob.become(self.get_param_func(axes, timer.get_value(), eplus, eminus, ogpsiquantum, evectorplus, evectorminus))
        )

        psithetatext = str(np.ceil(1000*psitheta)/1000)
        psiphitext = str(np.ceil(1000*psiphi)/1000)
        bthetatext = str(np.ceil(1000*btheta)/1000)
        bphitext = str(np.ceil(1000*bphi)/1000)

        psitext = Tex("$|\\psi(0)\\rangle = cos(\\frac{" + str(psithetatext) + "}{2})|+\\rangle + sin(\\frac{" + str(psithetatext) + "}{2})e^{i(" + str(psiphitext) + ")}|-\\rangle$").scale(0.5)
        btext = Tex("$|+\\rangle_B$", "$ = cos(\\frac{" + str(bthetatext) + "}{2})|+\\rangle + sin(\\frac{" + str(bthetatext) + "}{2})e^{i(" + str(bphitext) + ")}|-\\rangle$").scale(0.5).next_to(psitext, DOWN)
        btext[0].set_color(BLUE)
        vg = VGroup(psitext, btext)
        vg.to_corner(UR)

      
        self.add_fixed_in_frame_mobjects(vg)

        self.set_camera_orientation(phi=65 * DEGREES, theta=30 * DEGREES)
        self.add(sphere,axes, pf, bvector, psivector)
        self.wait(0.5)
        self.begin_ambient_camera_rotation(rate=0.1)
        #self.wait(.25)
        self.play(timer.animate(run_time=5, rate_func=linear).set_value(2*PI))










