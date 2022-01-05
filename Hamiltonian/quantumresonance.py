from manim import *
import numpy as np
import cmath

hbar = 2 #constants
eme = 1 #more constants
points = []
flipped = False

class VisualizingQunatumResonance(ThreeDScene):
    
    def sphereFunc(self, u, v):
        return np.array([np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v)])

    def hamilFuncButSpherical(self, t, brho, btheta, bphi, psitheta, psiphi, axes):
        eplus = eme*hbar/2*brho
        eminus = -1*eplus 
        evectorplus = np.array([[np.cos(btheta/2)], [np.sin(btheta/2) * cmath.exp(complex(0, bphi))]]) #defining eigenvectors (written in z basis)
        evectorminus = np.array([[np.sin(btheta/2)], [-np.cos(btheta/2) * cmath.exp(complex(0, bphi))]]) 


        psiquantum = np.array([[np.cos(btheta/2)*np.cos(psitheta/2)+np.sin(btheta/2)*np.sin(psitheta/2)*cmath.exp(complex(0, psiphi-bphi))], [np.sin(btheta/2)*np.cos(psitheta/2)-np.cos(btheta/2)*np.sin(psitheta/2)*cmath.exp(complex(0, psiphi-bphi))]])        
        
        thing = self.quantumToSpherical((cmath.exp(complex(0, -1*eplus*t/hbar))*psiquantum[0][0])*evectorplus + (cmath.exp(complex(0, -1*eminus*t/hbar))*psiquantum[1][0])*evectorminus)
        print (cmath.exp(complex(0, -1*eplus*t/hbar)))

        rect = self.sphericalToRectangular(thing[0].real, thing[1].real, thing[2].real)

        nextpoint = list(axes.c2p(rect[0], rect[1], rect[2]))
        if nextpoint[2] < 2.0:
            points.append(nextpoint)

        return np.array([thing[0].real, thing[1].real, thing[2].real])

    def hamilFunc(self, t, brho, btheta, bphi, psitheta, psiphi):
        eplus = eme*hbar/2*brho
        eminus = -1*eplus 
        evectorplus = np.array([[np.cos(btheta/2)], [np.sin(btheta/2) * cmath.exp(complex(0, bphi))]]) #defining eigenvectors (written in z basis)
        evectorminus = np.array([[np.sin(btheta/2)], [-np.cos(btheta/2) * cmath.exp(complex(0, bphi))]]) 


        psiquantum = np.array([[np.cos(btheta/2)*np.cos(psitheta/2)+np.sin(btheta/2)*np.sin(psitheta/2)*cmath.exp(complex(0, psiphi-bphi))], [np.sin(btheta/2)*np.cos(psitheta/2)-np.cos(btheta/2)*np.sin(psitheta/2)*cmath.exp(complex(0, psiphi-bphi))]])        
        
        thing = self.quantumToSpherical((cmath.exp(complex(0, -1*eplus*t/hbar))*psiquantum[0][0])*evectorplus + (cmath.exp(complex(0, -1*eminus*t/hbar))*psiquantum[1][0])*evectorminus)
        #print ((ogpsiquantum[0][0])*evectorplus)
        #print ("quantum state: " + str(cmath.exp(complex(0, -1*eplus*t/hbar))*ogpsiquantum[0][0])*evectorplus + (cmath.exp(complex(0, -1*eminus*t/hbar))*ogpsiquantum[1][0])*evectorminus)
        #print("spherical: " + str(thing))
        #print ("position: " + str(self.sphericalToRectangular(thing[0], thing[1], thing[2])))
        #print ("")
        return self.sphericalToRectangular(thing[0].real, thing[1].real, thing[2].real)

    '''def get_param_func(self, axes, t, brho, btheta, bphi, psitheta, psiphi):
        tol = 1e-9

        pf = ParametricFunction(
            lambda t: axes.c2p(*self.hamilFunc(t, brho, btheta, bphi, psitheta, psiphi)),
            t_range=[0, t, 0.007],
            tolerance_for_point_equality=tol
        )
        
        pf.set_color(color=[RED, YELLOW, BLUE, RED])
        return pf
'''

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

    def setup(self):
        axes = ThreeDAxes(x_range=(-2, 2, 1), y_range=(-2, 2, 1), z_range=(-2, 2, 1), x_length=8, y_length=8, z_length=8)
        points.append(list(axes.c2p(0, 0, 1)))


        
    def construct(self):       
        axes = ThreeDAxes(x_range=(-2, 2, 1), y_range=(-2, 2, 1), z_range=(-2, 2, 1), x_length=8, y_length=8, z_length=8) #defining axes
        sphere = Surface(
            lambda u, v: axes.c2p(*self.sphereFunc(u, v)), #make a sphere
            u_range=[-PI, PI],
            v_range=[0, TAU],
            fill_opacity=.1
        )
       
        timer = ValueTracker(0) #timer value

        drivefrequency = ValueTracker(0)
        bmodfrequency = 0
        bmodstrength = 1

        #bzero = 
        #bone = 
        #btwo = 

        #initialbrho = np.sqrt(bzero**2 + bone**2 + btwo**2)
        initialbrho = 1
        initialbtheta = PI/4
        initialbphi = 0


        brho = ValueTracker(initialbrho).add_updater(lambda m: m.set_value(initialbrho + bmodstrength * np.sin(timer.get_value() * bmodfrequency)))
        btheta = ValueTracker(initialbtheta)
        bphi = ValueTracker(initialbphi).add_updater(lambda m: m.set_value(initialbphi + drivefrequency.get_value()*timer.get_value()))

        ogbrectangular = self.sphericalToRectangular(initialbrho, initialbtheta, initialbphi)


        initialpsitheta = 0
        initialpsiphi = 0

        psiangles = Dot(point = [1, 0, 0], radius = 0).add_updater(lambda m: m.move_to([1, self.hamilFuncButSpherical(0.1, brho.get_value(), btheta.get_value(), bphi.get_value(), psiangles.get_y(), psiangles.get_z(), axes)[1], self.hamilFuncButSpherical(0.1, brho.get_value(), btheta.get_value(), bphi.get_value(), psiangles.get_y(), psiangles.get_z(), axes)[2]]))

        

        ''' psitheta = ValueTracker(initialpsitheta).add_updater(lambda m: m.set_value(
            self.hamilFuncButSpherical(timer.get_value(), brho.get_value(), btheta.get_value(), bphi.get_value(), psitheta.get_value(), psiphi.get_value())[1]
        ))

        psiphi = ValueTracker(initialpsiphi).add_updater(lambda m: m.set_value(
            self.hamilFuncButSpherical(timer.get_value(), brho.get_value(), btheta.get_value(), bphi.get_value(), psitheta.get_value(), psiphi.get_value())[2]
        ))'''
        
        ogpsirectangular = self.sphericalToRectangular(1.0, initialpsitheta, initialpsiphi)

        #changingpsiquantum = ogpsiquantum.add_updater(lambda m: m.become(
        #    cmath.exp(-i*eplus*timer.get_value()/hbar)*ogpsiquantum[0]*evectorplus + cmath.exp(-i*eminus*timer.get_value()/hbar)*ogpsiquantum[1]*evectorminus
        #))

        #psirectangular = self.sphericalToRectangular(self.quantumToSpherical(ogpsiquantum)).add_updater(lambda n: n.become(self.sphericalToRectangular(self.quantumToSpherical(changingpsiquantum))))

        

        bvector = Arrow((0, 0, 0), axes.c2p(ogbrectangular[0], ogbrectangular[1], ogbrectangular[2]), buff=0, color=BLUE_D).add_updater(lambda m: m.put_start_and_end_on([0, 0, 0], axes.c2p(self.sphericalToRectangular(brho.get_value(), btheta.get_value(), bphi.get_value())[0], self.sphericalToRectangular(brho.get_value(), btheta.get_value(), bphi.get_value())[1], self.sphericalToRectangular(brho.get_value(), btheta.get_value(), bphi.get_value())[2])))
            
        psivector = Arrow((0, 0, 0), axes.c2p(ogpsirectangular[0], ogpsirectangular[1], ogpsirectangular[2]), buff=0, color=WHITE).add_updater(lambda m: m.put_start_and_end_on([0,0,0], axes.c2p(self.hamilFunc(0.1, brho.get_value(), btheta.get_value(), bphi.get_value(), psiangles.get_y(), psiangles.get_z())[0], self.hamilFunc(0.1, brho.get_value(), btheta.get_value(), bphi.get_value(), psiangles.get_y(), psiangles.get_z())[1], self.hamilFunc(0.1, brho.get_value(), btheta.get_value(), bphi.get_value(), psiangles.get_y(), psiangles.get_z())[2]))) 

        


        
        
        #psivector = Arrow((0, 0, 0), axes.c2p(ogpsirectangular[0], ogpsirectangular[1], ogpsirectangular[2]), buff=0, color=WHITE).add_updater(lambda m: m.put_start_and_end_on([0,0,0], axes.c2p(self.hamilFunc(timer.get_value(), brho.get_value(), btheta.get_value(), bphi.get_value(), psitheta.get_value(), psiphi.get_value())[0], self.hamilFunc(timer.get_value(), brho.get_value(), btheta.get_value(), bphi.get_value(), psitheta.get_value(), psiphi.get_value())[1], self.hamilFunc(timer.get_value(), brho.get_value(), btheta.get_value(), bphi.get_value(), psitheta.get_value(), psiphi.get_value())[2]))) 

        #pf = self.get_param_func(axes, timer.get_value(), eplus, eminus, ogpsiquantum, evectorplus, evectorminus)   
        #pf.add_updater(lambda mob: mob.become(self.get_param_func(axes, timer.get_value(), eplus, eminus, ogpsiquantum, evectorplus, evectorminus)))

        psithetatext = str(np.ceil(1000*initialpsitheta)/1000)
        psiphitext = str(np.ceil(1000*initialpsiphi)/1000)
        bthetatext = str(np.ceil(1000*initialbtheta)/1000)
        bphitext = str(np.ceil(1000*initialbphi)/1000)

        psitext = Tex("$|\\psi(0)\\rangle = cos(\\frac{" + str(psithetatext) + "}{2})|+\\rangle + sin(\\frac{" + str(psithetatext) + "}{2})e^{i(" + str(psiphitext) + ")}|-\\rangle$").scale(0.5)
        btext = Tex("$|B(0)\\rangle$", "$ = cos(\\frac{" + str(bthetatext) + "}{2})|+\\rangle + sin(\\frac{" + str(bthetatext) + "}{2})e^{i(" + str(bphitext) + ")}|-\\rangle$").scale(0.5).next_to(psitext, DOWN)
        btext[0].set_color(BLUE)
        vg = VGroup(psitext, btext)
        vg.to_corner(UR)

        def get_param_func(group):
            #print(points)
            #print("a")
            newpath = VGroup(color = RED)
            newpath.set_points_smoothly(np.array(points))
            group.become(newpath)
        #return path

        path = VGroup().add_updater(get_param_func)
      
        self.add_fixed_in_frame_mobjects(vg)

        self.set_camera_orientation(phi=75 * DEGREES, theta=30 * DEGREES)
        self.add(sphere, axes, bvector, timer, brho, btheta, bphi, psivector, psiangles, path) #ADD PF BACK
        self.wait(.5)
        self.begin_ambient_camera_rotation(rate=0.3)
        self.wait(.25)
        for i in range (1, 20):
            drivefrequency.set_value(i/101)
            brho.set_value(initialbrho)
            btheta.set_value(initialbtheta)
            bphi.set_value(initialbphi)
            psiangles.set_x(1)
            psiangles.set_y(0)
            psiangles.set_z(0)

            timer.set_value(0)
            points.clear()
            points.append(list(axes.c2p(0, 0, 1)))
            flipped = False          
            self.play(timer.animate(run_time=2.5, rate_func=linear).set_value(25*PI))
            if flipped == True:
                print(i)            

        self.wait(.5)
