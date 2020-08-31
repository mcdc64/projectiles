# Cian McDonnell 2020
# Feel free to use or modify this code as you wish, but please credit me as the creator. Thanks.

import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np

import tkinter as Tk
from tkinter import ttk
from functools import partial
from scipy import integrate


#TO CHANGE PROJECTILE PROPERTIES
#READ THE ARGUMENTS AT THE TOP OF THE PROJECTILE CLASS
#THEN PUT THE ARGUMENTS YOU WANT INTO THE "projectiles" ON LINE 87/88



proj_size=1 #physical size of matplotlib marker to plot object position

def density(h): #make density a function of height for better accuracy (below about 10000m). All formulae taken from a NASA page.
    if(h<=11000):
        T0 = 288 #sea level temperature in K (used in the formula to determine density at a given altitude)
        alpha = 0.0065 #scaling factor determining how cold it gets as altitude goes up
        Th = T0 - (alpha*h)
        p0 = 101.29
        ph = p0* (((Th+273.1)/288)**5.256)
        density0 = 1.25 #density of fluid at sea level in kg/m^3
        densityh = density0*((Th/T0)**4.256)
    if(h>11000)and(h<=25000): #approximate density with an exponential relation here
        T = -56.46
        p = 22.65 * np.exp(1.73-0.000157*h) #pressure
        densityh = p/(0.2869*(T+273.1))
    if(h>25000): #atmosphere is so thin above this height that drag has little effect anyway
        T = -131.21+(0.00299*h)
        p = 2.488 * ((((T+273.1)/216.6))**-11.388)
        densityh = p/(0.2869*(T+273.1))
    return densityh
    
class Projectile:
    # arguments for a projectile: 1: inital x                        |  2: inital y
    #                             3. inital speed                    |  4. launching angle
    #                             5. drag coefficient                |  6. mass
    #                             7. projected area                  |  8. coefficient of restitution with ground
    #                             9. friction coefficient with ground| 10. vacuum? true/false
    def __init__(self,*args):
        self.init_state = [args[0], args[1], args[2] * np.cos(args[3]), args[2] * np.sin(args[3])]
        self.state = self.init_state #the state is the particle's x and y positions and velocities at any time
        self.c_drag = args[4] # drag coefficient
        self.mass = args[5]
        self.area = args[6]
        self.origin = [args[0],args[1]]
        self.time_elapsed = 0
        self.e = args[7]
        self.f = args[8]
        self.vacuum = args[9]
    def position(self):
        return [self.state[0],self.state[1]]

    def dstate_dt(self,state,t):
        s = self.state

        dsdt = np.zeros_like(self.state)
        dsdt[0] = state[2] #derivatives of the position components are merely the velocity components
        dsdt[1] = state[3]
        if(not self.vacuum):
            dsdt[2] = 0-(np.sign(self.state[2])*((self.c_drag*(self.state[2]**2))*(density(self.state[1])*self.area/self.mass))) # x velocity only experiences drag
            dsdt[3] = -g-(np.sign(self.state[3])*((self.c_drag*(self.state[3]**2))*(density(self.state[1])*self.area/self.mass))) #y velocity is changed by both gravity and drag
        else:
            dsdt[2] = 0 # x velocity constant
            dsdt[3] = -g #y velocity changes as t changes
            
        return dsdt
    def step(self,dt):
        self.state = integrate.odeint(self.dstate_dt,self.state,[0,dt])[1] #solves the equation system given in dstate_dt and gives an integrated result
        if(self.state[1]<=0.13*+proj_size and self.state[3]<=0.13*proj_size): #make the particle bounce when near the ground
            self.state[3] = -self.e*self.state[3] #assume the particle bounces back with velocity -ev
            self.state[2]= self.state[2]/self.f
        self.time_elapsed += dt

paused = False

atmo_hidden = False
vac_hidden = False

init_v = 100 #initial velocity in m/s
init_angle = np.pi/3.8 #launch angle in radians

projectile = Projectile(100,0,init_v,init_angle,0.1,3,0.04,0.6,1.4,False)
vac_projectile = Projectile(100,0,init_v,init_angle,0,3,0,0.6,1.4,True)

xdata = []
ydata = []
xdata_vac = []
ydata_vac = []


def pause(): #use a boolean to decide whether to continue the sim or not
    global paused
    paused = not paused

def hide_vac(): #hide the vacuum projectile if we only want to look at the atmospheric one
    global vac_hidden
    vac_hidden = not vac_hidden
def hide_atmo(): #hide the atmospheric projectile if we only want to look at the vacuum one
    global atmo_hidden
    atmo_hidden = not atmo_hidden
def reset(): #give the projectiles new position and velocity when the user enters some.
    global ax
    global xdata,ydata,xdata_vac,ydata_vac
    new_x = float(reset_x_pos.get())
    new_y = float(reset_y_pos.get())
    new_vx = 0
    new_vy = 0
    if(coordinate_var.get()=="Polar"):
        new_vx = float(reset_vel_1.get())*np.cos(float(reset_vel_2.get())*np.pi/180)
        new_vy = float(reset_vel_1.get())*np.sin(float(reset_vel_2.get())*np.pi/180)
    if(coordinate_var.get()=="Cartesian"):
        new_vx = float(reset_vel_1.get())
        new_vy = float(reset_vel_2.get())
        
    xdata = [] #clear position data for the projectiles so we don't get lines jumping from old position to the new one
    ydata = []
    xdata_vac = []
    ydata_vac = []
    
    projectile.time_elapsed = 0
    vac_projectile.time_elapsed = 0
        
    projectile.init_state = [new_x,new_y,new_vx,new_vy]
    projectile.state = [new_x,new_y,new_vx,new_vy]
    vac_projectile.init_state = [new_x,new_y,new_vx,new_vy]
    vac_projectile.state = [new_x,new_y,new_vx,new_vy]

    
    

fig = plt.figure(figsize=[10,10])

root = Tk.Tk()
root.geometry("1200x1000") #set window size
root.title("Projectiles")
frame = Tk.Frame(root)

canvas = FigureCanvasTkAgg(fig, master=root)
toolbar = NavigationToolbar2Tk( canvas, frame )

time_var = Tk.StringVar(root)
time_label = ttk.Label(textvariable=time_var)

pause_button = ttk.Button(text="Pause",command=partial(pause))

coordinate_var = Tk.StringVar(root)
coord_option = ttk.OptionMenu(root,coordinate_var,"Cartesian",*{"Cartesian":"Cartesian","Polar":"Polar"})
coord_option_label = ttk.Label(text="Coordinate System")

vel_comp_1 = Tk.StringVar(root) #velocity components to display. Can be either Cartesian or polar.
vel_comp_2 = Tk.StringVar(root)
vel_comp_1_vac = Tk.StringVar(root)
vel_comp_2_vac = Tk.StringVar(root)

vel_comp_1_label = ttk.Label(textvariable=vel_comp_1)
vel_comp_2_label = ttk.Label(textvariable=vel_comp_2)
vel_comp_1_vac_label = ttk.Label(textvariable=vel_comp_1_vac)
vel_comp_2_vac_label = ttk.Label(textvariable=vel_comp_2_vac)

reset_label = ttk.Label(text="Restart with New Position and Velocity:")

reset_x_pos = ttk.Entry(root) #Entries to reset the projectile position and velocity.
reset_x_pos_label = ttk.Label(text="New Initial X position (Cartesian):")
reset_y_pos = ttk.Entry(root)
reset_y_pos_label = ttk.Label(text="New Initial Y position (Cartesian):")

reset_vel_1_var = Tk.StringVar(root)
reset_vel_2_var = Tk.StringVar(root)
reset_vel_1 = ttk.Entry(root)
reset_vel_1_label = ttk.Label(textvariable=reset_vel_1_var)
reset_vel_2 = ttk.Entry(root)
reset_vel_2_label = ttk.Label(textvariable=reset_vel_2_var)

reset_button = ttk.Button(text="Restart",command=reset)

hide_atmo_button = ttk.Button(text="Show/Hide Atmospheric Projectile",command=hide_atmo)
hide_vac_button = ttk.Button(text="Show/Hide Vacuum Projectile",command=hide_vac)

frame.place(x=250,y=0)
canvas.get_tk_widget().place(x=250,y=30)
pause_button.place(x=10,y=100)
time_label.place(x=10,y=130)

coord_option_label.place(x=145,y=100)
coord_option.place(x=145,y=120)

vel_comp_1_label.place(x=10,y=170)
vel_comp_2_label.place(x=10,y=190)
vel_comp_1_vac_label.place(x=10,y=210)
vel_comp_2_vac_label.place(x=10,y=230)

reset_label.place(x=10,y=270)
reset_x_pos_label.place(x=10,y=290)
reset_x_pos.place(x=10,y=310)
reset_y_pos_label.place(x=10,y=340)
reset_y_pos.place(x=10,y=360)

reset_vel_1_label.place(x=10,y=390)
reset_vel_1.place(x=10,y=410)
reset_vel_2_label.place(x=10,y=440)
reset_vel_2.place(x=10,y=460)

reset_button.place(x=10,y=490)


hide_atmo_button.place(x=10,y=560)
hide_vac_button.place(x=10,y=590)


g = 9.81 #acceleration due to gravity in m/(s^2)
dt = 1.0/100

max_h = (init_v**2 * np.sin(init_angle)**2)/(2*g)
range = (2*(init_v**2)*np.tan(init_angle))/(g*(1+(np.tan(init_angle))**2))




ax = plt.axes()
ax.grid(which="both")
line1, = ax.plot([], [], "b",ms=proj_size)
line2, = ax.plot([], [], "r--",ms=proj_size)
ground,= ax.plot([-1e6,1e6],[0,0],"#000000",ms=0.3)

plt.legend([line1,line2], ['Projectile in Atmosphere', 'Projectile in Vacuum'],loc="upper right")

frame_no = 30

plt.xlim(0, 1000)
plt.ylim(0, 1000)
def init():
    if(not paused):
        line1.set_data([],[])
        line2.set_data([], [])
        ground.set_data([-1e6,1e6],[0,0])
        ax.set_aspect('equal')



    return line1,line2,ground,


def animate(n):
    if(not paused):
        xdata.append(projectile.position()[0])
        ydata.append(projectile.position()[1])
        xdata_vac.append(vac_projectile.position()[0])
        ydata_vac.append(vac_projectile.position()[1])
        projectile.step(dt)
        vac_projectile.step(dt)


        time_var.set("Time: "+str(round(projectile.time_elapsed,3))+" s")
    if(coordinate_var.get()=="Cartesian"):
        vel_comp_1.set("X Velocity: "+str(round(projectile.state[2],2))+" m/s")
        vel_comp_2.set("Y Velocity: "+str(round(projectile.state[3],2))+" m/s")
        vel_comp_1_vac.set("X Velocity (Vacuum): "+str(round(vac_projectile.state[2],2))+" m/s")
        vel_comp_2_vac.set("Y Velocity (Vacuum): "+str(round(vac_projectile.state[3],2))+" m/s")
        reset_vel_1_var.set("New Initial X Velocity (m/s): ")
        reset_vel_2_var.set("New Initial Y Velocity (m/s): ")
    if(coordinate_var.get()=="Polar"):
        vel_comp_1.set("Speed: "+ str(round(np.sqrt(projectile.state[2]**2+projectile.state[3]**2),2))+" m/s")
        vel_comp_2.set("Angle: "+ str(round(np.arctan(projectile.state[3]/projectile.state[2])*(180/np.pi),2))+" deg")
        vel_comp_1_vac.set("Speed (Vacuum): "+ str(round(np.sqrt(vac_projectile.state[2]**2+vac_projectile.state[3]**2),2))+" m/s")
        vel_comp_2_vac.set("Angle (Vacuum): "+ str(round(np.arctan(vac_projectile.state[3]/vac_projectile.state[2])*(180/np.pi),2))+" deg")
        reset_vel_1_var.set("New Initial Speed (m/s): ")
        reset_vel_2_var.set("New Launch Angle (deg): ")
    if(not atmo_hidden):
        line1.set_data(xdata,ydata)
    else:
        line1.set_data([],[])
        
    if(not vac_hidden):
        line2.set_data(xdata_vac,ydata_vac)
    else:
        line2.set_data([],[])
        
    return line1,line2,ground,
time_acceleration = 1
interv = int(dt*1000)
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=frame_no, interval=int((interv)/time_acceleration), blit=True, repeat=True)

canvas.draw() #update the tkinter canvas with the matplotlib animation

Tk.mainloop()