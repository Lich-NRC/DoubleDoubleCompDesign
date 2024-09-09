import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
import re
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from DD_Laminate_Finder_Functions import *

app = tk.Tk()
#app.geometry("500x500")
#app.tk.call('tk', 'scaling', 1.5)
app.title("Double-Double Laminate Finder")
app.option_add("*Font", "Aerial 11")

custom_style = ttk.Style()
custom_style.configure('.', font=('Aerial', 11))

#s = ttk.Style(app)
#s.theme_use("vista")

# Configure the tabs background color
custom_style.configure("TNotebook.Tab", background="light blue", padding=(10, 5))  # Tab background

# Global variables
phi_A = 0
psi_A = 0
phi_D = 0
psi_D = 0
phi_mid = 0
psi_mid = 0
outside_range_A = False
outside_range_D = False
outside_range_mid = False
Tsai = 0
theta_degrees_quad = []
theta_degrees_dd_a = []
theta_degrees_dd_d = []
theta_degrees_dd_mid = []
va1_quad = 0
va2_quad = 0
vd1_quad = 0
vd2_quad = 0
va1_dd_equi_a = 0
va2_dd_equi_a = 0
vd1_dd_equi_a = 0
vd2_dd_equi_a = 0
va1_dd_equi_d = 0
va2_dd_equi_d = 0
vd1_dd_equi_d = 0
vd2_dd_equi_d = 0
va1_dd_equi_mid = 0
va2_dd_equi_mid = 0
vd1_dd_equi_mid = 0
vd2_dd_equi_mid = 0
checkbox_symmetry_variable = tk.BooleanVar()
checkbox_symmetry_variable.set(True)
checkbox_material_variable = tk.BooleanVar()
checkbox_material_variable.set(False)
Ex = 200; Ey = 199; Es = 76.9; vx = 0.3; thickness=0.135;
index = 0

#########################################################
###################  Functions ##########################
#########################################################
# Function to update the theta_degrees list
def quad_input():
    global theta_degrees_quad
    
    # Initialize an empty theta_degrees list
    theta_degrees_quad = []
    theta_degrees_quad_temp  = []
    input_values = textbox_stacking.get("1.0",tk.END)
    try:
        # Split the input by either space or comma
        values = re.split(r'[,\s]+', input_values)
        theta_degrees_quad.extend(map(float, values[0:-1]))
        theta_degrees_quad_size = len(theta_degrees_quad) 
        if checkbox_symmetry_variable.get() == True:
            theta_degrees_quad_temp = theta_degrees_quad
            for i in range(1,theta_degrees_quad_size+1):
                theta_degrees_quad.append(theta_degrees_quad_temp[theta_degrees_quad_size-i])
    except ValueError:
        messagebox.showinfo("Invalid input", "Please enter valid numbers separated by spaces or commas")

def update_material_label(event):
    global Ex, Ey, Es, vx, thickness
    global index
    if checkbox_material_variable.get() == False:
        if listbox_material.curselection():
            index = listbox_material.curselection()[0] # Get the index of the selected item
        if index == 0:
            label_material.config(text="Ex = 164 GPa \nEy = 8.97 GPa\
                                    \nEs = 5.67 GPa \nvx = 0.32\nPly Thickness = 0.135 mm \nTsai's Modulus = Not Calculated")
            Ex = 164; Ey = 8.97; Es = 5.67; vx = 0.32; thickness = 0.135;
        elif index == 1:
            label_material.config(text="Ex = 181 GPa \nEy = 10.3 GPa\
                                    \nEs = 7.17 GPa \nvx = 0.28\nPly Thickness = 0.127 mm \nTsai's Modulus = Not Calculated")
            Ex = 181; Ey = 10.3; Es = 7.17; vx = 0.28; thickness = 0.127;
        elif index == 2:
            label_material.config(text="Ex = 162 GPa \nEy = 9.00 GPa\
                                    \nEs = 5.00 GPa \nvx = 0.40\nPly Thickness = 0.185 mm \nTsai's Modulus = Not Calculated")
            Ex = 162; Ey = 9.00; Es = 5.00; vx = 0.40; thickness = 0.185;
        elif index == 3:
            label_material.config(text="Ex = 159 GPa \nEy = 8.96 GPa\
                                    \nEs = 5.50 GPa \nvx = 0.32\nPly Thickness = 0.125 mm \nTsai's Modulus = Not Calculated")
            Ex = 159; Ey = 8.96; Es = 5.50; vx = 0.32; thickness = 0.125;
        elif index == 4:
            label_material.config(text="Ex = 138 GPa \nEy = 8.96 GPa\
                                    \nEs = 7.10 GPa \nvx = 0.30\nPly Thickness = 0.125 mm \nTsai's Modulus = Not Calculated")
            Ex = 138; Ey = 8.96; Es = 7.10; vx = 0.30; thickness = 0.125;
        elif index == 5:
            label_material.config(text="Ex = 126 GPa \nEy = 8.40 GPa\
                                    \nEs = 4.20 GPa \nvx = 0.31\nPly Thickness = 0.152 mm \nTsai's Modulus = Not Calculated")
            Ex = 126; Ey = 8.40; Es = 4.20; vx = 0.31; thickness = 0.152;


def update_ui(theta_degrees_input):
    global psi_A, phi_A, psi_D, phi_D, psi_mid, phi_mid
    global outside_range_A, outside_range_D, outside_range_mid
    global Tsai
    global theta_degrees_quad, theta_degrees_dd_a, theta_degrees_dd_d, theta_degrees_dd_mid
    global va1_quad ,va2_quad, vd1_quad ,vd2_quad
    global va1_dd_equi_a, va2_dd_equi_a, vd1_dd_equi_a, vd2_dd_equi_a
    global va1_dd_equi_d, va2_dd_equi_d, vd1_dd_equi_d, vd2_dd_equi_d
    global va1_dd_equi_mid, va2_dd_equi_mid, vd1_dd_equi_mid, vd2_dd_equi_mid
    global index
    global Ex, Ey, Es, vx, thickness

    global value_homogenization_AD_A 
    global value_homogenization_AD_D 
    global value_homogenization_AD_mid 
    global value_homogenization_B_A 
    global value_homogenization_B_D 
    global value_homogenization_B_mid 

#    value_homogenization_AD_A = []
#    value_homogenization_AD_D = []
#    value_homogenization_AD_mid = []
#    value_homogenization_B_A = []
#    value_homogenization_B_D = []
#    value_homogenization_B_mid = []

    outside_range_A = False
    outside_range_D = False
    outside_range_mid = False
    v1_star_mid = 0
    v2_star_mid = 0

    A, B, D, A_star, B_star, D_star, va_1_star, va_2_star, vd_1_star, vd_2_star = stiffness_matrix(theta_degrees_input=theta_degrees_quad, Ex=Ex, Ey=Ey, Es=Es, vx=vx, thickness=thickness)
    A_star_quad = A_star
    B_star_quad = B_star
    D_star_quad = D_star
    
    v1_star_mid = (va_1_star * (1 - scrollbar_mid.get()) + vd_1_star * scrollbar_mid.get() )
    v2_star_mid = (va_2_star * (1 - scrollbar_mid.get()) + vd_2_star * scrollbar_mid.get() )

    if all(value == 0 for value in theta_degrees_input) or all(value == 90 for value in theta_degrees_input):
        messagebox.showinfo("Invalid input", "Equivalent DD cannot be found for 0° or 90° unidirectional laminates")
    else:
        # DD Calculations
        Tsai = A_star[0, 0] + A_star[1, 1] + 2 * A_star[2, 2]
        if checkbox_material_variable.get() == False:
            if index == 0:
                label_material.config(text="Ex = 164 GPa \nEy = 8.97 GPa\
                                    \nEs = 5.67 GPa \nvx = 0.32\nPly Thickness = 0.135 mm \nTsai's Modulus = "+str(round(Tsai/10**9*100)/100)+" GPa")
            elif index == 1:
                label_material.config(text="Ex = 181 GPa \nEy = 10.3 GPa\
                                    \nEs = 7.17 GPa \nvx = 0.28\nPly Thickness = 0.127 mm \nTsai's Modulus = "+str(round(Tsai/10**9*100)/100)+" GPa")
            elif index == 2:
                label_material.config(text="Ex = 162 GPa \nEy = 9.00 GPa\
                                    \nEs = 5.00 GPa \nvx = 0.40\nPly Thickness = 0.185 mm \nTsai's Modulus = "+str(round(Tsai/10**9*100)/100)+" GPa")
            elif index == 3:
                label_material.config(text="Ex = 159 GPa \nEy = 8.96 GPa\
                                    \nEs = 5.50 GPa \nvx = 0.32\nPly Thickness = 0.125 mm \nTsai's Modulus = "+str(round(Tsai/10**9*100)/100)+" GPa")
            elif index == 4:
                label_material.config(text="Ex = 138 GPa \nEy = 8.96 GPa\
                                    \nEs = 7.10 GPa \nvx = 0.30\nPly Thickness = 0.125 mm \nTsai's Modulus = "+str(round(Tsai/10**9*100)/100)+" GPa")
            elif index == 5:
                label_material.config(text="Ex = 126 GPa \nEy = 8.40 GPa\
                                    \nEs = 4.20 GPa \nvx = 0.31\nPly Thickness = 0.152 mm \nTsai's Modulus = "+str(round(Tsai/10**9*100)/100)+" GPa")

        # Checking the Double-Double range 
        solve_psi_A = va_1_star + np.sqrt(-1 * va_1_star ** 2 + va_2_star / 2 + 0.5)
        solve_psi_D = vd_1_star + np.sqrt(-1 * vd_1_star ** 2 + vd_2_star / 2 + 0.5)
        solve_psi_mid = v1_star_mid + np.sqrt(-1 * v1_star_mid ** 2 + v2_star_mid / 2 + 0.5)
        if solve_psi_A > 1:
            solve_psi_A = 1
            outside_range_A = True
        if solve_psi_D > 1:
            solve_psi_D = 1
            outside_range_D = True
        if solve_psi_mid > 1:
            solve_psi_mid = 1
            outside_range_mid = True
        psi_A = 180 / np.pi * 0.5 * np.arccos(solve_psi_A)
        psi_D = 180 / np.pi * 0.5 * np.arccos(solve_psi_D)
        psi_mid = 180 / np.pi * 0.5 * np.arccos(solve_psi_mid)
        solve_phi_A = 2 * va_1_star - np.cos(2 * psi_A * np.pi / 180)
        solve_phi_D = 2 * vd_1_star - np.cos(2 * psi_D * np.pi / 180)
        solve_phi_mid = 2 * v1_star_mid - np.cos(2 * psi_mid * np.pi / 180)
        if solve_phi_A < -1:
            solve_phi_A = -1
            outside_range_A = True
        if solve_phi_D < -1:
            solve_phi_D = -1
            outside_range_D = True
        if solve_phi_mid < -1:
            solve_phi_mid = -1
            outside_range_mid = True
        phi_A = 180 / np.pi * 0.5 * np.arccos(solve_phi_A)
        phi_D = 180 / np.pi * 0.5 * np.arccos(solve_phi_D)
        phi_mid = 180 / np.pi * 0.5 * np.arccos(solve_phi_mid)

        label_equivalentdd_equalA.config(text="ψ: "+str(round(2 * psi_A) / 2)+"°\nφ: "+str(round(2 * phi_A) / 2)+"°", font=('Aerial', 15, 'bold'), foreground='blue', bg='#f5e690')
        label_equivalentdd_equalD.config(text="ψ: "+str(round(2 * psi_D) / 2)+"°\nφ: "+str(round(2 * phi_D) / 2)+"°", font=('Aerial', 15, 'bold'), foreground='blue', bg='#bae8e3')
        label_equivalentdd_equalmid.config(text="ψ: "+str(round(2 * psi_mid) / 2)+"°\nφ: "+str(round(2 * phi_mid) / 2)+"°", font=('Aerial', 15, 'bold'), foreground='blue', bg='#bae8b0')
        
        if outside_range_A==True:
            label_equivalentdd_equalA.config(text="ψ: "+str(round(2 * psi_A) / 2)+"°\nφ: "+str(round(2 * phi_A) / 2)+"°\n(Outside DD range)",\
                                            font=('Aerial', 15, 'bold'), foreground='blue', bg='#f5e690')
        if outside_range_D==True:    
            label_equivalentdd_equalD.config(text="ψ: "+str(round(2 * psi_D) / 2)+"°\nφ: "+str(round(2 * phi_D) / 2)+"°\n(Outside DD range)",\
                                            font=('Aerial', 15, 'bold'), foreground='blue', bg='#bae8e3')
        if outside_range_mid==True:    
            label_equivalentdd_equalmid.config(text="ψ: "+str(round(2 * psi_mid) / 2)+"°\nφ: "+str(round(2 * phi_mid) / 2)+"°\n(Outside DD range)",\
                                            font=('Aerial', 15, 'bold'), foreground='blue', bg='#bae8b0')
        
        for i in range(3):
            for j in range(3):
                value = A[i][j]
                formatted_value = f"{value:.2e}"
                if abs(A[i][j])<1e-4:
                    label_stiffness_quad_A = ttk.Label(frame_stiffness_quad_A, text=formatted_value, foreground='red')
                else:
                    label_stiffness_quad_A = ttk.Label(frame_stiffness_quad_A, text=formatted_value, padding=2)
                label_stiffness_quad_A.grid(row=i, column=j)
                value = D[i][j]
                formatted_value = f"{value:.2e}"
                if abs(D[i][j])<1e-4:
                    label_stiffness_quad_D = ttk.Label(frame_stiffness_quad_D, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_quad_D = ttk.Label(frame_stiffness_quad_D, text=formatted_value, padding=2)
                label_stiffness_quad_D.grid(row=i, column=j)
                value = B[i][j]
                formatted_value = f"{value:.2e}"
                if abs(B[i][j])<1e-4:
                    label_stiffness_quad_B = ttk.Label(frame_stiffness_quad_B, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_quad_B = ttk.Label(frame_stiffness_quad_B, text=formatted_value, padding=2)
                label_stiffness_quad_B.grid(row=i, column=j)
                
                # A*, B* and D*
                value = A_star[i][j]
                formatted_value = f"{value:.2e}"
                if abs(A_star[i][j])<1e-4:
                    label_stiffness_quad_A_star = ttk.Label(frame_stiffness_quad_A_star, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_quad_A_star = ttk.Label(frame_stiffness_quad_A_star, text=formatted_value, padding=2)
                label_stiffness_quad_A_star.grid(row=i, column=j)
                value = D_star[i][j]
                formatted_value = f"{value:.2e}"
                if abs(D_star[i][j])<1e-4:
                    label_stiffness_quad_D_star = ttk.Label(frame_stiffness_quad_D_star, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_quad_D_star = ttk.Label(frame_stiffness_quad_D_star, text=formatted_value, padding=2)
                label_stiffness_quad_D_star.grid(row=i, column=j)
                value = B_star[i][j]
                formatted_value = f"{value:.2e}"
                if abs(B_star[i][j])<1e-4:
                    label_stiffness_quad_B_star = ttk.Label(frame_stiffness_quad_B_star, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_quad_B_star = ttk.Label(frame_stiffness_quad_B_star, text=formatted_value, padding=2)
                label_stiffness_quad_B_star.grid(row=i, column=j)


        
        ##############################################################################################
        ############################### Equivalent DD with similar [A*] #############################
        ##############################################################################################

        if (entry_dd_a_repeat.get() == "") or (int(entry_dd_a_repeat.get()) < 1):
            messagebox.showinfo("Invalid input", "Please enter a valid number for the number of repeats of DD building block")

        theta_degrees_dd_a = []
        if int(radio_var_a.get()) == 1:
            for i in range(int(entry_dd_a_repeat.get())):
                theta_degrees_dd_a.append(phi_A)
                theta_degrees_dd_a.append(-psi_A)
                theta_degrees_dd_a.append(-phi_A)
                theta_degrees_dd_a.append(psi_A)
        if int(radio_var_a.get()) == 2:
            for i in range(int(entry_dd_a_repeat.get())):
                theta_degrees_dd_a.append(phi_A)
                theta_degrees_dd_a.append(psi_A)
                theta_degrees_dd_a.append(-phi_A)
                theta_degrees_dd_a.append(-psi_A)
        if int(radio_var_a.get()) == 3:
            for i in range(int(entry_dd_a_repeat.get())):
                theta_degrees_dd_a.append(phi_A)
                theta_degrees_dd_a.append(-psi_A)
                theta_degrees_dd_a.append(psi_A)
                theta_degrees_dd_a.append(-phi_A)
        if int(radio_var_a.get()) == 4:
            for i in range(int(entry_dd_a_repeat.get())):
                theta_degrees_dd_a.append(phi_A)
                theta_degrees_dd_a.append(-phi_A)
                theta_degrees_dd_a.append(psi_A)
                theta_degrees_dd_a.append(-psi_A)

        A, B, D, A_star, B_star, D_star, va_1_star, va_2_star, vd_1_star, vd_2_star = stiffness_matrix(theta_degrees_input=theta_degrees_dd_a, Ex=Ex, Ey=Ey, Es=Es, vx=vx, thickness=thickness)
        for i in range(3):
            for j in range(3):
                value = A_star[i][j]
                formatted_value = f"{value:.2e}"
                if abs(A_star[i][j])<1e-4:
                    label_stiffness_dd_a_A_star = ttk.Label(frame_stiffness_dd_a_A_star, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_dd_a_A_star = ttk.Label(frame_stiffness_dd_a_A_star, text=formatted_value, padding=2)
                label_stiffness_dd_a_A_star.grid(row=i, column=j)
                value = D_star[i][j]
                formatted_value = f"{value:.2e}"
                if abs(D_star[i][j])<1e-4:
                    label_stiffness_dd_a_D_star = ttk.Label(frame_stiffness_dd_a_D_star, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_dd_a_D_star = ttk.Label(frame_stiffness_dd_a_D_star, text=formatted_value, padding=2)
                label_stiffness_dd_a_D_star.grid(row=i, column=j)
                value = B_star[i][j]
                formatted_value = f"{value:.2e}"
                if abs(B_star[i][j])<1e-4:
                    label_stiffness_dd_a_B_star = ttk.Label(frame_stiffness_dd_a_B_star, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_dd_a_B_star = ttk.Label(frame_stiffness_dd_a_B_star, text=formatted_value, padding=2)
                label_stiffness_dd_a_B_star.grid(row=i, column=j)

        # Difference with quad for equivalent [A*]
        for i in range(3):
            for j in range(3):
                value = abs(A_star[i][j] - A_star_quad[i][j] ) / A_star_quad[i][j]  * 100
                if (value < 10000) and (value > -10000):
                    formatted_value = "   {:.{prec}f}".format(value, prec=(2 if value>=10 else 3))
                    label_stiffness_dd_a_A_star_diff = ttk.Label(frame_stiffness_dd_a_A_star_diff, text=formatted_value, padding=2)
                elif (value > 10000):
                    formatted_value = "   > 1e5   "
                    label_stiffness_dd_a_A_star_diff = ttk.Label(frame_stiffness_dd_a_A_star_diff, text=formatted_value, padding=2, foreground='blue')
                else:
                    formatted_value = "   < -1e5  "
                    label_stiffness_dd_a_A_star_diff = ttk.Label(frame_stiffness_dd_a_A_star_diff, text=formatted_value, padding=2, foreground='blue')
                label_stiffness_dd_a_A_star_diff.grid(row=i, column=j)
                value = abs(D_star[i][j] - D_star_quad[i][j] ) / D_star_quad[i][j]  * 100
                if (value < 10000) and (value > -10000):
                    formatted_value = "   {:.{prec}f}".format(value, prec=(2 if value>=10 else 3))
                    label_stiffness_dd_a_D_star_diff = ttk.Label(frame_stiffness_dd_a_D_star_diff, text=formatted_value, padding=2)
                elif (value > 10000):
                    formatted_value = "   > 1e5   "
                    label_stiffness_dd_a_D_star_diff = ttk.Label(frame_stiffness_dd_a_D_star_diff, text=formatted_value, padding=2, foreground='blue')
                else:
                    formatted_value = "   < -1e5  "
                    label_stiffness_dd_a_D_star_diff = ttk.Label(frame_stiffness_dd_a_D_star_diff, text=formatted_value, padding=2, foreground='blue')
                label_stiffness_dd_a_D_star_diff.grid(row=i, column=j)
                value = abs(B_star[i][j] - B_star_quad[i][j] ) / B_star_quad[i][j] * 100
                if (value < 10000) and (value > -10000):
                    formatted_value = "   {:.{prec}f}".format(value, prec=(2 if value>=10 else 3))
                    label_stiffness_dd_a_B_star_diff = ttk.Label(frame_stiffness_dd_a_B_star_diff, text=formatted_value, padding=2)
                elif (value > 10000):
                    formatted_value = "   > 1e5   "
                    label_stiffness_dd_a_B_star_diff = ttk.Label(frame_stiffness_dd_a_B_star_diff, text=formatted_value, padding=2, foreground='blue')
                else:
                    formatted_value = "   < -1e5  "
                    label_stiffness_dd_a_B_star_diff = ttk.Label(frame_stiffness_dd_a_B_star_diff, text=formatted_value, padding=2, foreground='blue')
                label_stiffness_dd_a_B_star_diff.grid(row=i, column=j)
        
        # Homogenization conditions check for equivalent [A*]
        value_homogenization_AD_A = abs(A_star - D_star) / Tsai * 100
        value_homogenization_B_A = abs(B_star) / Tsai * 100
        for i in range(3):
            for j in range(3):
                formatted_value = "{:.{prec}f}".format(value_homogenization_AD_A[i][j], prec=(2 if value_homogenization_AD_A[i][j]>=10 else 3))
                if abs(value_homogenization_AD_A[i][j]) <= 2:
                    label_stiffness_dd_a_homo_ad = ttk.Label(frame_stiffness_dd_a_homo_ad, text=formatted_value, padding=2, foreground='green',background='')
                else:
                    label_stiffness_dd_a_homo_ad = ttk.Label(frame_stiffness_dd_a_homo_ad, text=formatted_value, padding=2, background='#f5a99f')
                label_stiffness_dd_a_homo_ad.grid(row=i, column=j)

                formatted_value = "{:.{prec}f}".format(value_homogenization_B_A[i][j], prec=(2 if value_homogenization_B_A[i][j]>=10 else 3))
                if abs(value_homogenization_B_A[i][j]) <= 2:
                    label_stiffness_dd_a_homo_b = ttk.Label(frame_stiffness_dd_a_homo_b, text=formatted_value, padding=2, foreground='green',background='')
                else:
                    label_stiffness_dd_a_homo_b = ttk.Label(frame_stiffness_dd_a_homo_b, text=formatted_value, padding=2, background='#f5a99f')
                label_stiffness_dd_a_homo_b.grid(row=i, column=j)

        ##############################################################################################
        ############################# Equivalent DD with similar [D*] ###############################
        ##############################################################################################
        if (entry_dd_d_repeat.get() == "") or (int(entry_dd_d_repeat.get()) < 1):
            messagebox.showinfo("Invalid input", "Please enter a valid number for the number of repeats of DD building block")

        theta_degrees_dd_d = []
        if int(radio_var_d.get()) == 1:
            for i in range(int(entry_dd_d_repeat.get())):
                theta_degrees_dd_d.append(phi_D)
                theta_degrees_dd_d.append(-psi_D)
                theta_degrees_dd_d.append(-phi_D)
                theta_degrees_dd_d.append(psi_D)
        if int(radio_var_d.get()) == 2:
            for i in range(int(entry_dd_d_repeat.get())):
                theta_degrees_dd_d.append(phi_D)
                theta_degrees_dd_d.append(psi_D)
                theta_degrees_dd_d.append(-phi_D)
                theta_degrees_dd_d.append(-psi_D)
        if int(radio_var_d.get()) == 3:
            for i in range(int(entry_dd_d_repeat.get())):
                theta_degrees_dd_d.append(phi_D)
                theta_degrees_dd_d.append(-psi_D)
                theta_degrees_dd_d.append(psi_D)
                theta_degrees_dd_d.append(-phi_D)
        if int(radio_var_d.get()) == 4:
            for i in range(int(entry_dd_d_repeat.get())):
                theta_degrees_dd_d.append(phi_D)
                theta_degrees_dd_d.append(-phi_D)
                theta_degrees_dd_d.append(psi_D)
                theta_degrees_dd_d.append(-psi_D)
        A, B, D, A_star, B_star, D_star, va_1_star, va_2_star, vd_1_star, vd_2_star = stiffness_matrix(theta_degrees_input=theta_degrees_dd_d, Ex=Ex, Ey=Ey, Es=Es, vx=vx, thickness=thickness)
        for i in range(3):
            for j in range(3):
                value = A_star[i][j]
                formatted_value = f"{value:.2e}"
                if abs(A_star[i][j])<1e-4:
                    label_stiffness_dd_d_A_star = ttk.Label(frame_stiffness_dd_d_A_star, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_dd_d_A_star = ttk.Label(frame_stiffness_dd_d_A_star, text=formatted_value, padding=2)
                label_stiffness_dd_d_A_star.grid(row=i, column=j)
                value = D_star[i][j]
                formatted_value = f"{value:.2e}"
                if abs(D_star[i][j])<1e-4:
                    label_stiffness_dd_d_D_star = ttk.Label(frame_stiffness_dd_d_D_star, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_dd_d_D_star = ttk.Label(frame_stiffness_dd_d_D_star, text=formatted_value, padding=2)
                label_stiffness_dd_d_D_star.grid(row=i, column=j)
                value = B_star[i][j]
                formatted_value = f"{value:.2e}"
                if abs(B_star[i][j])<1e-4:
                    label_stiffness_dd_d_B_star = ttk.Label(frame_stiffness_dd_d_B_star, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_dd_d_B_star = ttk.Label(frame_stiffness_dd_d_B_star, text=formatted_value, padding=2)
                label_stiffness_dd_d_B_star.grid(row=i, column=j)
        
        # Difference with quad for equivalent [D*]
        for i in range(3):
            for j in range(3):
                value = abs(A_star[i][j] - A_star_quad[i][j] ) / A_star_quad[i][j]  * 100
                if (value < 10000) and (value > -10000):
                    formatted_value = "   {:.{prec}f}".format(value, prec=(2 if value>=10 else 3))
                    label_stiffness_dd_d_A_star_diff = ttk.Label(frame_stiffness_dd_d_A_star_diff, text=formatted_value, padding=2)
                elif (value > 10000):
                    formatted_value = "   > 1e5   "
                    label_stiffness_dd_d_A_star_diff = ttk.Label(frame_stiffness_dd_d_A_star_diff, text=formatted_value, padding=2, foreground='blue')
                else:
                    formatted_value = "   < -1e5  "
                    label_stiffness_dd_d_A_star_diff = ttk.Label(frame_stiffness_dd_d_A_star_diff, text=formatted_value, padding=2, foreground='blue')
                label_stiffness_dd_d_A_star_diff.grid(row=i, column=j)
                value = abs(D_star[i][j] - D_star_quad[i][j] ) / D_star_quad[i][j]  * 100
                if (value < 10000) and (value > -10000):
                    formatted_value = "   {:.{prec}f}".format(value, prec=(2 if value>=10 else 3))
                    label_stiffness_dd_d_D_star_diff = ttk.Label(frame_stiffness_dd_d_D_star_diff, text=formatted_value, padding=2)
                elif (value > 10000):
                    formatted_value = "   > 1e5   "
                    label_stiffness_dd_d_D_star_diff = ttk.Label(frame_stiffness_dd_d_D_star_diff, text=formatted_value, padding=2, foreground='blue')
                else:
                    formatted_value = "   < -1e5  "
                    label_stiffness_dd_d_D_star_diff = ttk.Label(frame_stiffness_dd_d_D_star_diff, text=formatted_value, padding=2, foreground='blue')
                label_stiffness_dd_d_D_star_diff.grid(row=i, column=j)
                value = abs(B_star[i][j] - B_star_quad[i][j] ) / B_star_quad[i][j] * 100
                if (value < 10000) and (value > -10000):
                    formatted_value = "   {:.{prec}f}".format(value, prec=(2 if value>=10 else 3))
                    label_stiffness_dd_d_B_star_diff = ttk.Label(frame_stiffness_dd_d_B_star_diff, text=formatted_value, padding=2)
                elif (value > 10000):
                    formatted_value = "   > 1e5   "
                    label_stiffness_dd_d_B_star_diff = ttk.Label(frame_stiffness_dd_d_B_star_diff, text=formatted_value, padding=2, foreground='blue')
                else:
                    formatted_value = "   < -1e5  "
                    label_stiffness_dd_d_B_star_diff = ttk.Label(frame_stiffness_dd_d_B_star_diff, text=formatted_value, padding=2, foreground='blue')
                label_stiffness_dd_d_B_star_diff.grid(row=i, column=j)
        
        # Homogenization conditions check for equivalent [D*]
        value_homogenization_AD_D = abs(A_star - D_star) / Tsai * 100
        value_homogenization_B_D = abs(B_star) / Tsai * 100
        for i in range(3):
            for j in range(3):
                formatted_value = "{:.{prec}f}".format(value_homogenization_AD_D[i][j], prec=(2 if value_homogenization_AD_D[i][j]>=10 else 3))
                if abs(value_homogenization_AD_D[i][j]) <= 2:
                    label_stiffness_dd_d_homo_ad = ttk.Label(frame_stiffness_dd_d_homo_ad, text=formatted_value, padding=2, foreground='green')
                else:
                    label_stiffness_dd_d_homo_ad = ttk.Label(frame_stiffness_dd_d_homo_ad, text=formatted_value, padding=2, background='#f5a99f')
                label_stiffness_dd_d_homo_ad.grid(row=i, column=j)

                prec=(2 if value_homogenization_B_D[i][j]>=10 else 3)
                formatted_value = "{:.{prec}f}".format(value_homogenization_B_D[i][j], prec=(2 if value_homogenization_B_D[i][j]>=10 else 3))
                if abs(value_homogenization_B_D[i][j]) <= 2:
                    label_stiffness_dd_d_homo_b = ttk.Label(frame_stiffness_dd_d_homo_b, text=formatted_value, padding=2, foreground='green')
                else:
                    label_stiffness_dd_d_homo_b = ttk.Label(frame_stiffness_dd_d_homo_b, text=formatted_value, padding=2, background='#f5a99f')
                label_stiffness_dd_d_homo_b.grid(row=i, column=j)

        ###############################################################################################
        ##################### Equivalent DD with a similar [A*] and [D*] ##############################
        ###############################################################################################
        if (entry_dd_mid_repeat.get() == "") or (int(entry_dd_mid_repeat.get()) < 1):
            messagebox.showinfo("Invalid input", "Please enter a valid number for the number of repeats of DD building block")

        theta_degrees_dd_mid = []
        if int(radio_var_mid.get()) == 1:
            for i in range(int(entry_dd_mid_repeat.get())):
                theta_degrees_dd_mid.append(phi_mid)
                theta_degrees_dd_mid.append(-psi_mid)
                theta_degrees_dd_mid.append(-phi_mid)
                theta_degrees_dd_mid.append(psi_mid)
        if int(radio_var_mid.get()) == 2:
            for i in range(int(entry_dd_mid_repeat.get())):
                theta_degrees_dd_mid.append(phi_mid)
                theta_degrees_dd_mid.append(psi_mid)
                theta_degrees_dd_mid.append(-phi_mid)
                theta_degrees_dd_mid.append(-psi_mid)
        if int(radio_var_mid.get()) == 3:
            for i in range(int(entry_dd_mid_repeat.get())):
                theta_degrees_dd_mid.append(phi_mid)
                theta_degrees_dd_mid.append(-psi_mid)
                theta_degrees_dd_mid.append(psi_mid)
                theta_degrees_dd_mid.append(-phi_mid)
        if int(radio_var_mid.get()) == 4:
            for i in range(int(entry_dd_mid_repeat.get())):
                theta_degrees_dd_mid.append(phi_mid)
                theta_degrees_dd_mid.append(-phi_mid)
                theta_degrees_dd_mid.append(psi_mid)
                theta_degrees_dd_mid.append(-psi_mid)
        A, B, D, A_star, B_star, D_star, va_1_star, va_2_star, vd_1_star, vd_2_star = stiffness_matrix(theta_degrees_input=theta_degrees_dd_mid, Ex=Ex, Ey=Ey, Es=Es, vx=vx, thickness=thickness)
        for i in range(3):
            for j in range(3):
                value = A_star[i][j]
                formatted_value = f"{value:.2e}"
                if abs(A_star[i][j])<1e-4:
                    label_stiffness_dd_mid_A_star = ttk.Label(frame_stiffness_dd_mid_A_star, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_dd_mid_A_star = ttk.Label(frame_stiffness_dd_mid_A_star, text=formatted_value, padding=2)
                label_stiffness_dd_mid_A_star.grid(row=i, column=j)
                value = D_star[i][j]
                formatted_value = f"{value:.2e}"
                if abs(D_star[i][j])<1e-4:
                    label_stiffness_dd_mid_D_star = ttk.Label(frame_stiffness_dd_mid_D_star, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_dd_mid_D_star = ttk.Label(frame_stiffness_dd_mid_D_star, text=formatted_value, padding=2)
                label_stiffness_dd_mid_D_star.grid(row=i, column=j)
                value = B_star[i][j]
                formatted_value = f"{value:.2e}"
                if abs(B_star[i][j])<1e-4:
                    label_stiffness_dd_mid_B_star = ttk.Label(frame_stiffness_dd_mid_B_star, text=formatted_value, padding=2, foreground='red')
                else:
                    label_stiffness_dd_mid_B_star = ttk.Label(frame_stiffness_dd_mid_B_star, text=formatted_value, padding=2)
                label_stiffness_dd_mid_B_star.grid(row=i, column=j)
        
        # Difference with quad for equivalent [A*] and [D*]
        for i in range(3):
            for j in range(3):
                value = abs(A_star[i][j] - A_star_quad[i][j] ) / A_star_quad[i][j]  * 100
                if (value < 10000) and (value > -10000):
                    formatted_value = "   {:.{prec}f}".format(value, prec=(2 if value>=10 else 3))
                    label_stiffness_dd_mid_A_star_diff = ttk.Label(frame_stiffness_dd_mid_A_star_diff, text=formatted_value, padding=2)
                elif (value > 10000):
                    formatted_value = "   > 1e5   "
                    label_stiffness_dd_mid_A_star_diff = ttk.Label(frame_stiffness_dd_mid_A_star_diff, text=formatted_value, padding=2, foreground='blue')
                else:
                    formatted_value = "   < -1e5  "
                    label_stiffness_dd_mid_A_star_diff = ttk.Label(frame_stiffness_dd_mid_A_star_diff, text=formatted_value, padding=2, foreground='blue')
                label_stiffness_dd_mid_A_star_diff.grid(row=i, column=j)
                value = abs(D_star[i][j] - D_star_quad[i][j] ) / D_star_quad[i][j]  * 100
                if (value < 10000) and (value > -10000):
                    formatted_value = "   {:.{prec}f}".format(value, prec=(2 if value>=10 else 3))
                    label_stiffness_dd_mid_D_star_diff = ttk.Label(frame_stiffness_dd_mid_D_star_diff, text=formatted_value, padding=2)
                elif (value > 10000):
                    formatted_value = "   > 1e5   "
                    label_stiffness_dd_mid_D_star_diff = ttk.Label(frame_stiffness_dd_mid_D_star_diff, text=formatted_value, padding=2, foreground='blue')
                else:
                    formatted_value = "   < -1e5  "
                    label_stiffness_dd_mid_D_star_diff = ttk.Label(frame_stiffness_dd_mid_D_star_diff, text=formatted_value, padding=2, foreground='blue')
                label_stiffness_dd_mid_D_star_diff.grid(row=i, column=j)
                value = abs(B_star[i][j] - B_star_quad[i][j] ) / B_star_quad[i][j] * 100
                if (value < 10000) and (value > -10000):
                    formatted_value = "   {:.{prec}f}".format(value, prec=(2 if value>=10 else 3))
                    label_stiffness_dd_mid_B_star_diff = ttk.Label(frame_stiffness_dd_mid_B_star_diff, text=formatted_value, padding=2)
                elif (value > 10000):
                    formatted_value = "   > 1e5   "
                    label_stiffness_dd_mid_B_star_diff = ttk.Label(frame_stiffness_dd_mid_B_star_diff, text=formatted_value, padding=2, foreground='blue')
                else:
                    formatted_value = "   < -1e5  "
                    label_stiffness_dd_mid_B_star_diff = ttk.Label(frame_stiffness_dd_mid_B_star_diff, text=formatted_value, padding=2, foreground='blue')
                label_stiffness_dd_mid_B_star_diff.grid(row=i, column=j)
        
        # Homogenization conditions check for equivalent [A*] and [D*]
        value_homogenization_AD_mid = abs(A_star - D_star) / Tsai * 100
        value_homogenization_B_mid = abs(B_star) / Tsai * 100
        for i in range(3):
            for j in range(3):
                formatted_value = "{:.{prec}f}".format(value_homogenization_AD_mid[i][j], prec=(2 if value_homogenization_AD_mid[i][j]>=10 else 3))
                if abs(value_homogenization_AD_mid[i][j]) <= 2:
                    label_stiffness_dd_mid_homo_ad = ttk.Label(frame_stiffness_dd_mid_homo_ad, text=formatted_value, padding=2, foreground='green')
                else:
                    label_stiffness_dd_mid_homo_ad = ttk.Label(frame_stiffness_dd_mid_homo_ad, text=formatted_value, padding=2, background='#f5a99f')
                label_stiffness_dd_mid_homo_ad.grid(row=i, column=j)

                prec=(2 if value_homogenization_B_mid[i][j]>=10 else 3)
                formatted_value = "{:.{prec}f}".format(value_homogenization_B_mid[i][j], prec=(2 if value_homogenization_B_mid[i][j]>=10 else 3))
                if abs(value_homogenization_B_mid[i][j]) <= 2:
                    label_stiffness_dd_mid_homo_b = ttk.Label(frame_stiffness_dd_mid_homo_b, text=formatted_value, padding=2, foreground='green')
                else:
                    label_stiffness_dd_mid_homo_b = ttk.Label(frame_stiffness_dd_mid_homo_b, text=formatted_value, padding=2, background='#f5a99f')
                label_stiffness_dd_mid_homo_b.grid(row=i, column=j)

        
        # for Plotting the DD range and indicating the QUAD and the equivalent DD 
        A, B, D, A_star, B_star, D_star, va_1_star, va_2_star, vd_1_star, vd_2_star = stiffness_matrix(theta_degrees_input=theta_degrees_quad, Ex=Ex, Ey=Ey, Es=Es, vx=vx, thickness=thickness)
        va1_quad = va_1_star
        va2_quad = va_2_star
        vd1_quad = vd_1_star
        vd2_quad = vd_2_star
        A, B, D, A_star, B_star, D_star, va_1_star, va_2_star, vd_1_star, vd_2_star = stiffness_matrix(theta_degrees_input=theta_degrees_dd_a, Ex=Ex, Ey=Ey, Es=Es, vx=vx, thickness=thickness)
        va1_dd_equi_a = va_1_star
        va2_dd_equi_a = va_2_star
        vd1_dd_equi_a = vd_1_star
        vd2_dd_equi_a = vd_2_star
        A, B, D, A_star, B_star, D_star, va_1_star, va_2_star, vd_1_star, vd_2_star = stiffness_matrix(theta_degrees_input=theta_degrees_dd_d, Ex=Ex, Ey=Ey, Es=Es, vx=vx, thickness=thickness)
        va1_dd_equi_d = va_1_star
        va2_dd_equi_d = va_2_star
        vd1_dd_equi_d = vd_1_star
        vd2_dd_equi_d = vd_2_star
        A, B, D, A_star, B_star, D_star, va_1_star, va_2_star, vd_1_star, vd_2_star = stiffness_matrix(theta_degrees_input=theta_degrees_dd_mid, Ex=Ex, Ey=Ey, Es=Es, vx=vx, thickness=thickness)
        va1_dd_equi_mid = va_1_star
        va2_dd_equi_mid = va_2_star
        vd1_dd_equi_mid = vd_1_star
        vd2_dd_equi_mid = vd_2_star

    if checkbox_material_variable.get() == True:
            Ex = float(entry_material_Ex_text.get())
            Ey = float(entry_material_Ey_text.get())
            Es = float(entry_material_Es_text.get())
            vx = float(entry_material_vx_text.get())
            thickness = float(entry_material_thickness_text.get())
    
    listbox_material.selection_set(index)

def write_output():
    file_output_name = filedialog.asksaveasfilename(
        defaultextension=".txt",
        filetypes=[("Text files", "*.txt"), ("All files", "*.*")] )
    # Check if the user selected a file
    if file_output_name:
        with open(file_output_name, 'w', encoding='utf-8') as file_output:
            file_output.write('************************************************************************* \n')
            file_output.write('********************************* Inputs ******************************** \n')
            file_output.write('************************************************************************* \n')
            if checkbox_material_variable.get() == False:
                selected_items = [listbox_material.get(i) for i in listbox_material.curselection()]
                file_output.write('Material: ' + str(selected_items) + '\n')
            else:
                file_output.write('Material: ' + 'User-defined' + '\n')
            file_output.write('Ex: ' + str(Ex) + ' GPa\n')
            file_output.write('Ey: ' + str(Ey) + ' GPa\n')
            file_output.write('Es: ' + str(Es) + ' GPa\n')
            file_output.write('vx: ' + str(vx) + '\n')
            file_output.write('Ply Thickness: ' + str(thickness) + ' mm\n')
            file_output.write("Tsai's modulus: " + str(round(Tsai/1e9, 2)) + ' GPa\n')
            theta_degrees_quad_output = str(theta_degrees_quad)
            theta_degrees_quad_output = theta_degrees_quad_output.replace(",","/")
            theta_degrees_quad_output = theta_degrees_quad_output.replace(" ","")
            theta_degrees_quad_output = round_numbers_in_string(theta_degrees_quad_output)
            file_output.write("Stacking Sequence: " + theta_degrees_quad_output + '\n\n')
            file_output.write('************************************************************************* \n')
            file_output.write('******************************** DD Angles ****************************** \n')
            file_output.write('************************************************************************* \n')
            file_output.write('************************* DD with an equal [A*] ************************* \n')
            file_output.write('DD angles for a laminate with an equal [A] matrix: \n')
            file_output.write('ψ: ' + str(round(psi_A,1)) + '°\n')
            file_output.write('φ: ' + str(round(phi_A,1)) + '°\n')
            if int(radio_var_a.get()) == 1:
                file_output.write('Stacking Sequence: Staggered 1 [φ/-ψ/-φ/ψ]\n')
            elif int(radio_var_a.get()) == 2:
                file_output.write('Stacking Sequence: Staggered 2 [φ/ψ/-φ/-ψ]\n')
            elif int(radio_var_a.get()) == 3:
                file_output.write('Stacking Sequence: Staggered 3 [φ/-ψ/ψ/-φ]\n')
            elif int(radio_var_a.get()) == 4:
                file_output.write('Stacking Sequence: Paired [φ/-φ/ψ/-ψ]\n')
            file_output.write('Number of Repeats: ' + str(int(entry_dd_a_repeat.get())) + '\n')
            if (np.all(value_homogenization_AD_A < 2) and np.all(value_homogenization_B_A < 2)):
                file_output.write('Homogenization condition: Homogenized\n')
            else:
                file_output.write('Homogenization condition: Not Homogenized\n')
            file_output.write("\n([A*] - [D*])/Tsai (%): \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_AD_A[0][0], value_homogenization_AD_A[0][1], value_homogenization_AD_A[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_AD_A[1][0], value_homogenization_AD_A[1][1], value_homogenization_AD_A[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_AD_A[2][0], value_homogenization_AD_A[2][1], value_homogenization_AD_A[2][2]))
            file_output.write("\n[B*]/Tsai (%): \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_B_A[0][0], value_homogenization_B_A[0][1], value_homogenization_B_A[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_B_A[1][0], value_homogenization_B_A[1][1], value_homogenization_B_A[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_B_A[2][0], value_homogenization_B_A[2][1], value_homogenization_B_A[2][2]))
                    
            file_output.write('\n************************* DD with an equal [D*] ************************* \n')
            file_output.write('DD angles for a laminate with an equal [D] matrix: \n')
            file_output.write('ψ: ' + str(round(psi_D,1)) + '°\n')
            file_output.write('φ: ' + str(round(phi_D,1)) + '°\n\n')
            if int(radio_var_d.get()) == 1:
                file_output.write('Stacking Sequence: Staggered 1 [φ/-ψ/-φ/ψ]\n')
            elif int(radio_var_d.get()) == 2:
                file_output.write('Stacking Sequence: Staggered 2 [φ/ψ/-φ/-ψ]\n')
            elif int(radio_var_d.get()) == 3:
                file_output.write('Stacking Sequence: Staggered 3 [φ/-ψ/ψ/-φ]\n')
            elif int(radio_var_d.get()) == 4:
                file_output.write('Stacking Sequence: Paired [φ/-φ/ψ/-ψ]\n')
            file_output.write('Number of Repeats: ' + str(int(entry_dd_d_repeat.get())) + '\n')
            if (np.all(value_homogenization_AD_D < 2) and np.all(value_homogenization_B_D < 2)):
                file_output.write('Homogenization condition: Homogenized\n')
            else:
                file_output.write('Homogenization condition: Not Homogenized\n')
            file_output.write("\n([A*] - [D*])/Tsai (%): \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_AD_D[0][0], value_homogenization_AD_D[0][1], value_homogenization_AD_D[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_AD_D[1][0], value_homogenization_AD_D[1][1], value_homogenization_AD_D[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_AD_D[2][0], value_homogenization_AD_D[2][1], value_homogenization_AD_D[2][2]))
            file_output.write("\n[B*]/Tsai (%): \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_B_D[0][0], value_homogenization_B_D[0][1], value_homogenization_B_D[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_B_D[1][0], value_homogenization_B_D[1][1], value_homogenization_B_D[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_B_D[2][0], value_homogenization_B_D[2][1], value_homogenization_B_D[2][2]))
                    
            file_output.write('\n******************* DD with a similar [A*] and [D*] ********************* \n')
            file_output.write('DD angles for a laminate with a similar [A] and [D] matrix: \n')
            file_output.write('In-plane/Flexural equivalency: ' + str(scrollbar_mid.get()) + '\n')
            file_output.write('ψ: ' + str(round(psi_mid,1)) + '°\n')
            file_output.write('φ: ' + str(round(phi_mid,1)) + '°\n\n')
            if int(radio_var_mid.get()) == 1:
                file_output.write('Stacking Sequence: Staggered 1 [φ/-ψ/-φ/ψ]\n')
            elif int(radio_var_mid.get()) == 2:
                file_output.write('Stacking Sequence: Staggered 2 [φ/ψ/-φ/-ψ]\n')
            elif int(radio_var_mid.get()) == 3:
                file_output.write('Stacking Sequence: Staggered 3 [φ/-ψ/ψ/-φ]\n')
            elif int(radio_var_mid.get()) == 4:
                file_output.write('Stacking Sequence: Paired [φ/-φ/ψ/-ψ]\n')
            file_output.write('Number of Repeats: ' + str(int(entry_dd_mid_repeat.get())) + '\n')
            if (np.all(value_homogenization_AD_mid < 2) and np.all(value_homogenization_B_mid < 2)):
                file_output.write('Homogenization condition: Homogenized\n')
            else:
                file_output.write('Homogenization condition: Not Homogenized\n')
            file_output.write("\n([A*] - [D*])/Tsai (%): \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_AD_mid[0][0], value_homogenization_AD_mid[0][1], value_homogenization_AD_mid[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_AD_mid[1][0], value_homogenization_AD_mid[1][1], value_homogenization_AD_mid[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_AD_mid[2][0], value_homogenization_AD_mid[2][1], value_homogenization_AD_mid[2][2]))
            file_output.write("\n[B*]/Tsai (%): \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_B_mid[0][0], value_homogenization_B_mid[0][1], value_homogenization_B_mid[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_B_mid[1][0], value_homogenization_B_mid[1][1], value_homogenization_B_mid[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(value_homogenization_B_mid[2][0], value_homogenization_B_mid[2][1], value_homogenization_B_mid[2][2]))
                    

            file_output.write('\n************************************************************************* \n')
            file_output.write('****************************** ABD Matrices ***************************** \n')
            file_output.write('************************************************************************* \n')

            A, B, D, A_star, B_star, D_star, va_1_star, va_2_star, vd_1_star, vd_2_star = stiffness_matrix(theta_degrees_input=theta_degrees_quad, Ex=Ex, Ey=Ey, Es=Es, vx=vx, thickness=thickness)
            file_output.write('****************************** Given Laminate *************************** \n')
            file_output.write("Stiffness matrices for the given laminate: \n")
            file_output.write("[A]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A[0][0], A[0][1], A[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A[1][0], A[1][1], A[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(A[2][0], A[2][1], A[2][2]))
            file_output.write("[B]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B[0][0], B[0][1], B[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B[1][0], B[1][1], B[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(B[2][0], B[2][1], B[2][2]))
            file_output.write("[D]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D[0][0], D[0][1], D[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D[1][0], D[1][1], D[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(D[2][0], D[2][1], D[2][2]))

            file_output.write("Normalized stiffness matrices for the given laminate: \n")
            file_output.write("[A*]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A_star[0][0], A_star[0][1], A_star[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A_star[1][0], A_star[1][1], A_star[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(A_star[2][0], A_star[2][1], A_star[2][2]))
            file_output.write("[B*]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B_star[0][0], B_star[0][1], B_star[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B_star[1][0], B_star[1][1], B_star[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(B_star[2][0], B_star[2][1], B_star[2][2]))
            file_output.write("[D*]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D_star[0][0], D_star[0][1], D_star[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D_star[1][0], D_star[1][1], D_star[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(D_star[2][0], D_star[2][1], D_star[2][2]))

            A, B, D, A_star, B_star, D_star, va_1_star, va_2_star, vd_1_star, vd_2_star = stiffness_matrix(theta_degrees_input=theta_degrees_dd_a, Ex=Ex, Ey=Ey, Es=Es, vx=vx, thickness=thickness)
            file_output.write('************************* DD with an equal [A*] ************************* \n')
            file_output.write("Stiffness matrices for the DD laminate with an equal [A*]: \n")
            file_output.write("[A]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A[0][0], A[0][1], A[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A[1][0], A[1][1], A[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(A[2][0], A[2][1], A[2][2]))
            file_output.write("[B]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B[0][0], B[0][1], B[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B[1][0], B[1][1], B[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(B[2][0], B[2][1], B[2][2]))
            file_output.write("[D]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D[0][0], D[0][1], D[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D[1][0], D[1][1], D[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(D[2][0], D[2][1], D[2][2]))

            file_output.write("Normalized stiffness matrices for the DD laminate with an equal [A*]: \n")
            file_output.write("[A*]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A_star[0][0], A_star[0][1], A_star[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A_star[1][0], A_star[1][1], A_star[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(A_star[2][0], A_star[2][1], A_star[2][2]))
            file_output.write("[B*]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B_star[0][0], B_star[0][1], B_star[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B_star[1][0], B_star[1][1], B_star[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(B_star[2][0], B_star[2][1], B_star[2][2]))
            file_output.write("[D*]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D_star[0][0], D_star[0][1], D_star[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D_star[1][0], D_star[1][1], D_star[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(D_star[2][0], D_star[2][1], D_star[2][2]))

            A, B, D, A_star, B_star, D_star, va_1_star, va_2_star, vd_1_star, vd_2_star = stiffness_matrix(theta_degrees_input=theta_degrees_dd_d, Ex=Ex, Ey=Ey, Es=Es, vx=vx, thickness=thickness)
            file_output.write('************************* DD with an equal [D*] ************************* \n')
            file_output.write("Stiffness matrices for the DD laminate with an equal [D*]: \n")
            file_output.write("[A]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A[0][0], A[0][1], A[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A[1][0], A[1][1], A[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(A[2][0], A[2][1], A[2][2]))
            file_output.write("[B]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B[0][0], B[0][1], B[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B[1][0], B[1][1], B[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(B[2][0], B[2][1], B[2][2]))
            file_output.write("[D]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D[0][0], D[0][1], D[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D[1][0], D[1][1], D[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(D[2][0], D[2][1], D[2][2]))

            file_output.write("Normalized stiffness matrices for the DD laminate with an equal [D*]: \n")
            file_output.write("[A*]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A_star[0][0], A_star[0][1], A_star[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A_star[1][0], A_star[1][1], A_star[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(A_star[2][0], A_star[2][1], A_star[2][2]))
            file_output.write("[B*]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B_star[0][0], B_star[0][1], B_star[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B_star[1][0], B_star[1][1], B_star[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(B_star[2][0], B_star[2][1], B_star[2][2]))
            file_output.write("[D*]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D_star[0][0], D_star[0][1], D_star[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D_star[1][0], D_star[1][1], D_star[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(D_star[2][0], D_star[2][1], D_star[2][2]))

            A, B, D, A_star, B_star, D_star, va_1_star, va_2_star, vd_1_star, vd_2_star = stiffness_matrix(theta_degrees_input=theta_degrees_dd_mid, Ex=Ex, Ey=Ey, Es=Es, vx=vx, thickness=thickness)
            file_output.write('******************* DD with a similar [A*] and [D*] ********************* \n')
            file_output.write("Stiffness matrices for the DD laminate with a similar [A*] and [D*]: \n")
            file_output.write("[A]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A[0][0], A[0][1], A[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A[1][0], A[1][1], A[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(A[2][0], A[2][1], A[2][2]))
            file_output.write("[B]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B[0][0], B[0][1], B[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B[1][0], B[1][1], B[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(B[2][0], B[2][1], B[2][2]))
            file_output.write("[D]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D[0][0], D[0][1], D[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D[1][0], D[1][1], D[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(D[2][0], D[2][1], D[2][2]))

            file_output.write("Normalized stiffness matrices for the DD laminate with a similar [A*] and [D*]: \n")
            file_output.write("[A*]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A_star[0][0], A_star[0][1], A_star[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(A_star[1][0], A_star[1][1], A_star[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(A_star[2][0], A_star[2][1], A_star[2][2]))
            file_output.write("[B*]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B_star[0][0], B_star[0][1], B_star[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(B_star[1][0], B_star[1][1], B_star[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(B_star[2][0], B_star[2][1], B_star[2][2]))
            file_output.write("[D*]: \n")
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D_star[0][0], D_star[0][1], D_star[0][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n".format(D_star[1][0], D_star[1][1], D_star[1][2]))
            file_output.write("{:.2e} {:.2e} {:.2e}\n\n".format(D_star[2][0], D_star[2][1], D_star[2][2]))
    file_output.close()
    messagebox.showinfo("Saved", f'The report file "{file_output_name}" has been saved.')

def toggle_material():
    global Ex, Ey, Es, vx, thickness
    if checkbox_material_variable.get() == True:
        listbox_material.config(state=tk.DISABLED)
        label_material.config(state=tk.DISABLED)
        label_material_thickness.config(state=tk.NORMAL)
        label_material_Ex.config(state=tk.NORMAL)
        label_material_Ey.config(state=tk.NORMAL)
        label_material_Es.config(state=tk.NORMAL)
        label_material_vx.config(state=tk.NORMAL)
        entry_material_thickness.config(state=tk.NORMAL)
        entry_material_Ex.config(state=tk.NORMAL)
        entry_material_Ey.config(state=tk.NORMAL)
        entry_material_Es.config(state=tk.NORMAL)
        entry_material_vx.config(state=tk.NORMAL)
        Ex = float(entry_material_Ex_text.get())
        Ey = float(entry_material_Ey_text.get())
        Es = float(entry_material_Es_text.get())
        vx = float(entry_material_vx_text.get())
        thickness = float(entry_material_thickness_text.get())
    else:
        listbox_material.config(state=tk.NORMAL)
        label_material.config(state=tk.NORMAL)
        label_material_thickness.config(state=tk.DISABLED)
        label_material_Ex.config(state=tk.DISABLED)
        label_material_Ey.config(state=tk.DISABLED)
        label_material_Es.config(state=tk.DISABLED)
        label_material_vx.config(state=tk.DISABLED)
        entry_material_thickness.config(state=tk.DISABLED)
        entry_material_Ex.config(state=tk.DISABLED)
        entry_material_Ey.config(state=tk.DISABLED)
        entry_material_Es.config(state=tk.DISABLED)
        entry_material_vx.config(state=tk.DISABLED)

def button_click(event=None):
    # Simulate a click on the desired button
    button_start.invoke()
##################################################################################################################
##################################################################################################################
#########################################################  GUI ###################################################
##################################################################################################################
##################################################################################################################
frame_app = tk.Frame(app, padx=5, pady=2)

############  frame_material ############
frame_material = tk.LabelFrame(frame_app, padx=5, pady=2, text="Material Selection")
frame_material.grid(row=0, column=0, padx=5, pady=2, sticky='news')
listbox_material = tk.Listbox(frame_material, height=3, width=15)
listbox_material.grid(row=0, column=0, sticky='news')
listbox_material.insert(0, "IM7/977-3")
listbox_material.insert(1, "T300/5208")
listbox_material.insert(2, "T800/Cytec")
listbox_material.insert(3, "IM7/8552")
listbox_material.insert(4, "AS4/H3501")
listbox_material.insert(5, "T700/2510")
listbox_material.select_set(0)
label_material = tk.Label(frame_material, text="Ex = \nEy = \nEs = \nvx= \nPly Thickness = \nTsai's Modulus = ", justify='left')
label_material.grid(row=0, column=1, columnspan=2, sticky=tk.N+tk.W)
listbox_material.bind("<<ListboxSelect>>", update_material_label)

label_material_thickness = tk.Label(frame_material, text=f"Thickness (mm)")
label_material_thickness.grid(row=3, column=0, padx=1, pady=1)
label_material_Ex = tk.Label(frame_material, text=f"Ex (GPa)")
label_material_Ex.grid(row=1, column=1, padx=1, pady=1)
label_material_Ey = tk.Label(frame_material, text=f"Ey (GPa)")
label_material_Ey.grid(row=1, column=2, padx=1, pady=1)
label_material_Es = tk.Label(frame_material, text=f"Es (GPa)")
label_material_Es.grid(row=3, column=1, padx=1, pady=1)
label_material_vx = tk.Label(frame_material, text=f"vx")
label_material_vx.grid(row=3, column=2, padx=1, pady=1)

entry_material_thickness_text = tk.StringVar()
entry_material_Ex_text = tk.StringVar()
entry_material_Ey_text = tk.StringVar()
entry_material_Es_text = tk.StringVar()
entry_material_vx_text = tk.StringVar()

entry_material_thickness_text.initialize(0.135)
entry_material_Ex_text.initialize(164)
entry_material_Ey_text.initialize(8.97)
entry_material_Es_text.initialize(5.67)
entry_material_vx_text.initialize(0.32)

entry_material_thickness = tk.Entry(frame_material, textvariable=entry_material_thickness_text, width=15)
entry_material_thickness.grid(row=4, column=0, padx=1, pady=1)
entry_material_Ex = tk.Entry(frame_material, textvariable=entry_material_Ex_text, width=15)
entry_material_Ex.grid(row=2, column=1, padx=1, pady=1) 
entry_material_Ey = tk.Entry(frame_material, textvariable=entry_material_Ey_text, width=15)
entry_material_Ey.grid(row=2, column=2, padx=1, pady=1) 
entry_material_Es = tk.Entry(frame_material, textvariable=entry_material_Es_text, width=15)
entry_material_Es.grid(row=4, column=1, padx=1, pady=1) 
entry_material_vx = tk.Entry(frame_material, textvariable=entry_material_vx_text, width=15)
entry_material_vx.grid(row=4, column=2, padx=1, pady=1) 

label_material_thickness.config(state=tk.DISABLED)
label_material_Ex.config(state=tk.DISABLED)
label_material_Ey.config(state=tk.DISABLED)
label_material_Es.config(state=tk.DISABLED)
label_material_vx.config(state=tk.DISABLED)
entry_material_thickness.config(state=tk.DISABLED)
entry_material_Ex.config(state=tk.DISABLED)
entry_material_Ey.config(state=tk.DISABLED)
entry_material_Es.config(state=tk.DISABLED)
entry_material_vx.config(state=tk.DISABLED)

checkbox_material = tk.Checkbutton(frame_material, text="User-defined", variable=checkbox_material_variable, command=toggle_material)
checkbox_material.grid(row=2, column=0, pady=1, sticky=tk.N+tk.W)

############ frame_stacking ###############
frame_stacking = tk.LabelFrame(frame_app, padx=5, pady=2, text="Quad Stacking Sequence")
frame_stacking.grid(row=1, column=0, rowspan=3, padx=5, pady=2, sticky='news')
frame_stacking.columnconfigure(0, weight=1)
label_stacking_notes = tk.Label(frame_stacking, text="(Please separate angles with comma\n or space)", justify='left')
label_stacking_notes.grid(row=0, column=0, sticky=tk.N+tk.W)
textbox_stacking = tk.Text(frame_stacking, height=2, width=12)
textbox_stacking.grid(row=1, column=0,  columnspan=2, sticky='news')
textbox_stacking.insert("1.0", "0 45 90 -45 0 45 90 -45 0 45 90 -45 0 45 90 -45")
checkbox_symmetry = tk.Checkbutton(frame_stacking, text="Symmetry", variable=checkbox_symmetry_variable)
checkbox_symmetry.grid(row=0, column=1, pady=1, sticky=tk.N+tk.W)

############ frame_quad ###############
frame_quad = tk.LabelFrame(frame_app, padx=5, pady=2, text="Quad Laminate Stiffness")
frame_quad.grid(row=0, column=1, columnspan=1, rowspan=2, padx=5, pady=2, sticky='news')

frame_stiffness_quad_A = ttk.Labelframe(frame_quad, borderwidth=2, relief="solid", text="               [A] (Pa.m)               ")
frame_stiffness_quad_A.grid(row=0, column=2, padx=0, pady=2, sticky='news')
label_stiffness_quad_A = ttk.Label(frame_stiffness_quad_A, text="", padding=2)
label_stiffness_quad_A.grid()

frame_stiffness_quad_D = ttk.Labelframe(frame_quad, borderwidth=2, relief="solid", text="               [D] (Pa.m³)               ")
frame_stiffness_quad_D.grid(row=0, column=3, padx=0, pady=2, sticky='news')
label_stiffness_quad_D = ttk.Label(frame_stiffness_quad_D, text="", padding=2)
label_stiffness_quad_D.grid()

frame_stiffness_quad_B = ttk.Labelframe(frame_quad, borderwidth=2, relief="solid", text="               [B] (Pa.m²)               ")
frame_stiffness_quad_B.grid(row=0, column=4, padx=0, pady=2, sticky='news')
label_stiffness_quad_B = ttk.Label(frame_stiffness_quad_B, text="", padding=2)
label_stiffness_quad_B.grid()

frame_stiffness_quad_A_star = ttk.Labelframe(frame_quad, borderwidth=2, relief="solid", text="               [A*] (Pa)               ")
frame_stiffness_quad_A_star.grid(row=1, column=2, padx=0, pady=2, sticky='news')
label_stiffness_quad_A_star = ttk.Label(frame_stiffness_quad_A_star, text="", padding=2)
label_stiffness_quad_A_star.grid()

frame_stiffness_quad_D_star = ttk.Labelframe(frame_quad, borderwidth=2, relief="solid", text="               [D*] (Pa)               ")
frame_stiffness_quad_D_star.grid(row=1, column=3, padx=0, pady=2, sticky='news')
label_stiffness_quad_D_star = ttk.Label(frame_stiffness_quad_D_star, text="", padding=2)
label_stiffness_quad_D_star.grid()

frame_stiffness_quad_B_star = ttk.Labelframe(frame_quad, borderwidth=2, relief="solid", text="               [B*] (Pa)               ")
frame_stiffness_quad_B_star.grid(row=1, column=4, padx=0, pady=2, sticky='news')
label_stiffness_quad_B_star = ttk.Label(frame_stiffness_quad_B_star, text="", padding=2)
label_stiffness_quad_B_star.grid()

label_stiffness_notes1 = tk.Label(frame_quad, text="Normalized stiffnesses:  [A*] = [A] / h   ,   [D*] = 12[D] / h³   ,   [B*] = 2[B] / h²")
label_stiffness_notes1. grid(row=2, column=2, columnspan=3, pady=1, sticky=tk.W+tk.S)

frame_plot = tk.Frame(frame_app, padx=5, pady=2)
frame_plot.grid(row=2, column=1, columnspan=1, padx=5, pady=2, sticky='news')

button_equi_inplane_plot = tk.Button(frame_plot, text="Plot in-plane DD range", font=('Aerial', 11),\
                                      bg='#f5e690', command=lambda: [quad_input(), update_ui(theta_degrees_input=theta_degrees_quad),\
                                                                      plot_parabola(va1_quad, va2_quad, va1_dd_equi_a, va2_dd_equi_a, \
                                                                                    va1_dd_equi_d, va2_dd_equi_d, va1_dd_equi_mid, va2_dd_equi_mid, True)])
button_equi_inplane_plot.grid(row=0, column=0, padx=5, pady=2, sticky=tk.S, ipadx=75)

button_equi_flexural_plot = tk.Button(frame_plot, text="Plot flexural DD range", font=('Aerial', 11), bg='#f5e690', command=lambda: [quad_input(), update_ui(theta_degrees_input=theta_degrees_quad),\
                                                                      plot_parabola(vd1_quad, vd2_quad, vd1_dd_equi_a, vd2_dd_equi_a, \
                                                                                    vd1_dd_equi_d, vd2_dd_equi_d, vd1_dd_equi_mid, vd2_dd_equi_mid, False)])
button_equi_flexural_plot.grid(row=0, column=1, padx=5, pady=2, sticky=tk.S, ipadx=75)

button_start = tk.Button(frame_plot, text="Start / Update", font=('Aerial', 11, 'bold'), bg='#aef099', command=lambda: [quad_input(), update_ui(theta_degrees_input=theta_degrees_quad)]) 
button_start.grid(row=1, column=0, columnspan=1, padx=5, pady=2, sticky='news', ipadx=75)

button_write = tk.Button(frame_plot, text="Save Results", font=('Aerial', 11), bg='#f5e690', command=lambda: [quad_input(), update_ui(theta_degrees_input=theta_degrees_quad), write_output()]) 
button_write.grid(row=1, column=1, columnspan=1, padx=5, pady=2, sticky='news', ipadx=75)

# Bind the <Return> key event to the button_click function
app.bind('<Return>', button_click)
##################################################################################################################
################################################### Tabs #########################################################
##################################################################################################################
notebook_equivalent_dd = ttk. Notebook(frame_app)
notebook_equivalent_dd.grid(row=4, column=0, columnspan=2)
tab_equivalentdd_A = ttk.Frame(notebook_equivalent_dd)
tab_equivalentdd_D = ttk.Frame(notebook_equivalent_dd)
tab_equivalentdd_mid = ttk.Frame(notebook_equivalent_dd)
notebook_equivalent_dd.add(tab_equivalentdd_A, text="Equivalent DD with similar [A*]")
notebook_equivalent_dd.add(tab_equivalentdd_D, text="Equivalent DD with similar [D*]")
notebook_equivalent_dd.add(tab_equivalentdd_mid, text="Equivalent DD with similar [A*] and [D*]")

################################################################################
############################# frame_equivalentdd_A ###########################
################################################################################
frame_equivalentdd_A = tk.Frame(tab_equivalentdd_A, padx=5, pady=2)
frame_equivalentdd_A.grid(row=0, column=0, columnspan=2, padx=5, pady=2, sticky='news')

radio_var_a = tk.IntVar()
radio_button1_a = tk.Radiobutton(frame_equivalentdd_A, text="Staggered 1\n[φ/-ψ/-φ/ψ]", variable=radio_var_a, value=1)
radio_button2_a = tk.Radiobutton(frame_equivalentdd_A, text="Staggered 2\n[φ/ψ/-φ/-ψ]", variable=radio_var_a, value=2)
radio_button3_a = tk.Radiobutton(frame_equivalentdd_A, text="Staggered 3\n[φ/-ψ/ψ/-φ]", variable=radio_var_a, value=3)
radio_button4_a = tk.Radiobutton(frame_equivalentdd_A, text="Paired\n[φ/-φ/ψ/-ψ]", variable=radio_var_a, value=4)
radio_button1_a.grid(row=0, column=2, sticky=tk.W)
radio_button2_a.grid(row=0, column=3, sticky=tk.W)
radio_button3_a.grid(row=0, column=4, sticky=tk.W)
radio_button4_a.grid(row=0, column=5, sticky=tk.W)
radio_var_a.set(1)

label_entry_dd_a_repeat = tk.Label(frame_equivalentdd_A, text="Number of repeats: ")
label_entry_dd_a_repeat.grid(row=0, column=6, sticky=tk.W)
entry_dd_a_repeat = tk.Entry(frame_equivalentdd_A, width=5)
entry_dd_a_repeat.grid(row=0, column=7, sticky=tk.W)
entry_dd_a_repeat.insert(0, "3")

label_homogenization_a = tk.Label(frame_equivalentdd_A, text="Homogenization Check\nAll should be green (<2%)", font=('Aerial', 11, 'bold'))
label_homogenization_a.grid(row=0, column=8, columnspan=4)

frame_stiffness_dd_a_A_star = ttk.Labelframe(frame_equivalentdd_A, borderwidth=2, relief="solid", text="                 [A*] (Pa)                 ")
frame_stiffness_dd_a_A_star.grid(row=1, column=2, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_a_A_star = ttk.Label(frame_stiffness_dd_a_A_star, text="", padding=2)
label_stiffness_dd_a_A_star.grid()

frame_stiffness_dd_a_D_star = ttk.Labelframe(frame_equivalentdd_A, borderwidth=2, relief="solid", text="                 [D*] (Pa)                 ")
frame_stiffness_dd_a_D_star.grid(row=1, column=4, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_a_D_star = ttk.Label(frame_stiffness_dd_a_D_star, text="", padding=2)
label_stiffness_dd_a_D_star.grid()

frame_stiffness_dd_a_B_star = ttk.Labelframe(frame_equivalentdd_A, borderwidth=2, relief="solid", text="                 [B*] (Pa)                 ")
frame_stiffness_dd_a_B_star.grid(row=1, column=6, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_a_B_star = ttk.Label(frame_stiffness_dd_a_B_star, text="", padding=2)
label_stiffness_dd_a_B_star.grid()

frame_stiffness_dd_a_homo_ad = ttk.Labelframe(frame_equivalentdd_A, borderwidth=2, relief="solid", text="([A*] - [D*]) / Tsai (%)")
frame_stiffness_dd_a_homo_ad.grid(row=1, column=8, rowspan=4, columnspan=2, padx=5, pady=2, sticky='news')
label_stiffness_dd_a_homo_ad = ttk.Label(frame_stiffness_dd_a_homo_ad, text="", padding=2)
label_stiffness_dd_a_homo_ad.grid()

frame_stiffness_dd_a_homo_b = ttk.Labelframe(frame_equivalentdd_A, borderwidth=5, relief="solid", text="  [B*] / Tsai (%)  ")
frame_stiffness_dd_a_homo_b.grid(row=1, column=10, rowspan=4, columnspan=2, padx=5, pady=2, sticky='news')
label_stiffness_dd_a_homo_b = ttk.Label(frame_stiffness_dd_a_homo_b, text="", padding=2)
label_stiffness_dd_a_homo_b.grid()

# Difference with the Quad
frame_stiffness_dd_a_A_star_diff = ttk.Labelframe(frame_equivalentdd_A, borderwidth=2, relief="solid", text="Difference with quad [A*] (%)")
frame_stiffness_dd_a_A_star_diff.grid(row=5, column=2, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_a_A_star_diff = ttk.Label(frame_stiffness_dd_a_A_star_diff, text="", padding=2)
label_stiffness_dd_a_A_star_diff.grid()

frame_stiffness_dd_a_D_star_diff = ttk.Labelframe(frame_equivalentdd_A, borderwidth=2, relief="solid", text="Difference with quad [D*] (%)")
frame_stiffness_dd_a_D_star_diff.grid(row=5, column=4, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_a_D_star_diff = ttk.Label(frame_stiffness_dd_a_D_star_diff, text="", padding=2)
label_stiffness_dd_a_D_star_diff.grid()

frame_stiffness_dd_a_B_star_diff = ttk.Labelframe(frame_equivalentdd_A, borderwidth=2, relief="solid", text="Difference with quad [B*] (%)")
frame_stiffness_dd_a_B_star_diff.grid(row=5, column=6, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_a_B_star_diff = ttk.Label(frame_stiffness_dd_a_B_star_diff, text="", padding=2)
label_stiffness_dd_a_B_star_diff.grid()

label_equivalentdd_equalA = tk.Label(frame_equivalentdd_A, text="ψ: \nφ:", font=('Aerial', 15, 'bold'), foreground='blue', bg='#f5e690')
label_equivalentdd_equalA.grid(row=5, column=8, rowspan=4, columnspan=4, padx=5, pady=2,  sticky='news')

################################################################################
############################ frame_equivalentdd_D ################################
################################################################################

frame_equivalentdd_D = tk.Label(tab_equivalentdd_D, padx=5, pady=2)
frame_equivalentdd_D.grid(row=0, column=0, columnspan=2, padx=5, pady=2, sticky='news')

radio_var_d = tk.IntVar()
radio_button1_d = tk.Radiobutton(frame_equivalentdd_D, text="Staggered 1\n[φ/-ψ/-φ/ψ]", variable=radio_var_d, value=1)
radio_button2_d = tk.Radiobutton(frame_equivalentdd_D, text="Staggered 2\n[φ/ψ/-φ/-ψ]", variable=radio_var_d, value=2)
radio_button3_d = tk.Radiobutton(frame_equivalentdd_D, text="Staggered 3\n[φ/-ψ/ψ/-φ]", variable=radio_var_d, value=3)
radio_button4_d = tk.Radiobutton(frame_equivalentdd_D, text="Paired\n[φ/-φ/ψ/-ψ]", variable=radio_var_d, value=4)
radio_button1_d.grid(row=0, column=2, sticky=tk.W)
radio_button2_d.grid(row=0, column=3, sticky=tk.W)
radio_button3_d.grid(row=0, column=4, sticky=tk.W)
radio_button4_d.grid(row=0, column=5, sticky=tk.W)
radio_var_d.set(1)

label_entry_dd_d_repeat = tk.Label(frame_equivalentdd_D, text="Number of repeats: ")
label_entry_dd_d_repeat.grid(row=0, column=6, sticky=tk.W)
entry_dd_d_repeat = tk.Entry(frame_equivalentdd_D, width=5)
entry_dd_d_repeat.grid(row=0, column=7, sticky=tk.W)
entry_dd_d_repeat.insert(0, "3")

label_homogenization_d = tk.Label(frame_equivalentdd_D, text="Homogenization Check\nAll should be green (<2%)", font=('Aerial', 11, 'bold'))
label_homogenization_d.grid(row=0, column=8, columnspan=4)

frame_stiffness_dd_d_A_star = ttk.Labelframe(frame_equivalentdd_D, borderwidth=2, relief="solid", text="                 [A*] (Pa)                 ")
frame_stiffness_dd_d_A_star.grid(row=1, column=2, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_d_A_star = ttk.Label(frame_stiffness_dd_d_A_star, text="", padding=2)
label_stiffness_dd_d_A_star.grid()

frame_stiffness_dd_d_D_star = ttk.Labelframe(frame_equivalentdd_D, borderwidth=2, relief="solid", text="                 [D*] (Pa)                 ")
frame_stiffness_dd_d_D_star.grid(row=1, column=4, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_d_D_star = ttk.Label(frame_stiffness_dd_d_D_star, text="", padding=2)
label_stiffness_dd_d_D_star.grid()

frame_stiffness_dd_d_B_star = ttk.Labelframe(frame_equivalentdd_D, borderwidth=2, relief="solid", text="                 [B*] (Pa)                 ")
frame_stiffness_dd_d_B_star.grid(row=1, column=6, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_d_B_star = ttk.Label(frame_stiffness_dd_d_B_star, text="", padding=2)
label_stiffness_dd_d_B_star.grid()

frame_stiffness_dd_d_homo_ad = ttk.Labelframe(frame_equivalentdd_D, borderwidth=2, relief="solid", text="([A*] - [D*]) / Tsai (%)")
frame_stiffness_dd_d_homo_ad.grid(row=1, column=8, rowspan=4, columnspan=2, padx=5, pady=2, sticky='news')
label_stiffness_dd_d_homo_ad = ttk.Label(frame_stiffness_dd_d_homo_ad, text="", padding=2)
label_stiffness_dd_d_homo_ad.grid()

frame_stiffness_dd_d_homo_b = ttk.Labelframe(frame_equivalentdd_D, borderwidth=5, relief="solid", text="  [B*] / Tsai (%)  ")
frame_stiffness_dd_d_homo_b.grid(row=1, column=10, rowspan=4, columnspan=2, padx=5, pady=2, sticky='news')
label_stiffness_dd_d_homo_b = ttk.Label(frame_stiffness_dd_d_homo_b, text="", padding=2)
label_stiffness_dd_d_homo_b.grid()

# Difference with the Quad
frame_stiffness_dd_d_A_star_diff = ttk.Labelframe(frame_equivalentdd_D, borderwidth=2, relief="solid", text="Difference with quad [A*] (%)")
frame_stiffness_dd_d_A_star_diff.grid(row=5, column=2, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_d_A_star_diff = ttk.Label(frame_stiffness_dd_d_A_star_diff, text="", padding=2)
label_stiffness_dd_d_A_star_diff.grid()

frame_stiffness_dd_d_D_star_diff = ttk.Labelframe(frame_equivalentdd_D, borderwidth=2, relief="solid", text="Difference with quad [D*] (%)")
frame_stiffness_dd_d_D_star_diff.grid(row=5, column=4, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_d_D_star_diff = ttk.Label(frame_stiffness_dd_d_D_star_diff, text="", padding=2)
label_stiffness_dd_d_D_star_diff.grid()

frame_stiffness_dd_d_B_star_diff = ttk.Labelframe(frame_equivalentdd_D, borderwidth=2, relief="solid", text="Difference with quad [B*] (%)")
frame_stiffness_dd_d_B_star_diff.grid(row=5, column=6, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_d_B_star_diff = ttk.Label(frame_stiffness_dd_d_B_star_diff, text="", padding=2)
label_stiffness_dd_d_B_star_diff.grid()

label_equivalentdd_equalD = tk.Label(frame_equivalentdd_D, text="ψ: \nφ:", font=('Aerial', 15, 'bold'), foreground='blue', bg='#bae8e3')
label_equivalentdd_equalD.grid(row=5, column=8, rowspan=4, columnspan=4, padx=5, pady=2,  sticky='news')

######################################################################################
############################ frame_equivalentdd_mid ################################
######################################################################################
frame_equivalentdd_mid = tk.Label(tab_equivalentdd_mid, padx=5, pady=2)
frame_equivalentdd_mid.grid(row=0, column=0, columnspan=2, padx=5, pady=2, sticky='news')

frame_equivalentdd_mid_scroll = ttk.Labelframe(frame_equivalentdd_mid, text="In-plane / Flexural equivalency selection")
frame_equivalentdd_mid_scroll.grid(row=0, column=4, columnspan=4, sticky='news')

label_equivalentdd_mid_note2 = tk.Label(frame_equivalentdd_mid_scroll, text="Equal\n[A*]", justify='right')
label_equivalentdd_mid_note2.grid(row=0, column=0, columnspan=1, sticky=tk.E)

scrollbar_mid = tk.Scale(frame_equivalentdd_mid_scroll, orient=tk.HORIZONTAL, from_=0, to=1, resolution=0.05, length=350)
scrollbar_mid.set(0.5)  # Set the initial value of the scrollbar
scrollbar_mid.grid(row=0, column=1, columnspan=4, sticky='news')

label_equivalentdd_mid_note3 = tk.Label(frame_equivalentdd_mid_scroll, text="Equal\n[D*]", justify='left')
label_equivalentdd_mid_note3.grid(row=0, column=5, columnspan=1, sticky=tk.W)

radio_var_mid = tk.IntVar()
radio_button1_mid = tk.Radiobutton(frame_equivalentdd_mid, text="Staggered 1\n[φ/-ψ/-φ/ψ]", variable=radio_var_mid, value=1)
radio_button2_mid = tk.Radiobutton(frame_equivalentdd_mid, text="Staggered 2\n[φ/ψ/-φ/-ψ]", variable=radio_var_mid, value=2)
radio_button3_mid = tk.Radiobutton(frame_equivalentdd_mid, text="Staggered 3\n[φ/-ψ/ψ/-φ]", variable=radio_var_mid, value=3)
radio_button4_mid = tk.Radiobutton(frame_equivalentdd_mid, text="Paired\n[φ/-φ/ψ/-ψ]", variable=radio_var_mid, value=4)
radio_button1_mid.grid(row=1, column=2, sticky=tk.W)
radio_button2_mid.grid(row=1, column=3, sticky=tk.W)
radio_button3_mid.grid(row=1, column=4, sticky=tk.W)
radio_button4_mid.grid(row=1, column=5, sticky=tk.W)
radio_var_mid.set(1)

label_entry_dd_mid_repeat = tk.Label(frame_equivalentdd_mid, text="Number of repeats: ")
label_entry_dd_mid_repeat.grid(row=1, column=6, sticky=tk.W)
entry_dd_mid_repeat = tk.Entry(frame_equivalentdd_mid, width=5)
entry_dd_mid_repeat.grid(row=1, column=7, sticky=tk.W)
entry_dd_mid_repeat.insert(0, "3")

label_homogenization_mid = tk.Label(frame_equivalentdd_mid, text="Homogenization Check\nAll should be green (<2%)", font=('Aerial', 11, 'bold'))
label_homogenization_mid.grid(row=1, column=8, columnspan=4)

frame_stiffness_dd_mid_A_star = ttk.Labelframe(frame_equivalentdd_mid, borderwidth=2, relief="solid", text="                 [A*] (Pa)                 ")
frame_stiffness_dd_mid_A_star.grid(row=2, column=2, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_mid_A_star = ttk.Label(frame_stiffness_dd_mid_A_star, text="", padding=2)
label_stiffness_dd_mid_A_star.grid()

frame_stiffness_dd_mid_D_star = ttk.Labelframe(frame_equivalentdd_mid, borderwidth=2, relief="solid", text="                 [D*] (Pa)                 ")
frame_stiffness_dd_mid_D_star.grid(row=2, column=4, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_mid_D_star = ttk.Label(frame_stiffness_dd_mid_D_star, text="", padding=2)
label_stiffness_dd_mid_D_star.grid()

frame_stiffness_dd_mid_B_star = ttk.Labelframe(frame_equivalentdd_mid, borderwidth=2, relief="solid", text="                 [B*] (Pa)                 ")
frame_stiffness_dd_mid_B_star.grid(row=2, column=6, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_mid_B_star = ttk.Label(frame_stiffness_dd_mid_B_star, text="", padding=2)
label_stiffness_dd_mid_B_star.grid()

frame_stiffness_dd_mid_homo_ad = ttk.Labelframe(frame_equivalentdd_mid, borderwidth=2, relief="solid", text="([A*] - [D*]) / Tsai (%)")
frame_stiffness_dd_mid_homo_ad.grid(row=2, column=8, rowspan=4, columnspan=2, padx=5, pady=2, sticky='news')
label_stiffness_dd_mid_homo_ad = ttk.Label(frame_stiffness_dd_mid_homo_ad, text="", padding=2)
label_stiffness_dd_mid_homo_ad.grid()

frame_stiffness_dd_mid_homo_b = ttk.Labelframe(frame_equivalentdd_mid, borderwidth=5, relief="solid", text="  [B*] / Tsai (%)  ")
frame_stiffness_dd_mid_homo_b.grid(row=2, column=10, rowspan=4, columnspan=2, padx=5, pady=2, sticky='news')
label_stiffness_dd_mid_homo_b = ttk.Label(frame_stiffness_dd_mid_homo_b, text="", padding=2)
label_stiffness_dd_mid_homo_b.grid()

# Difference with the Quad
frame_stiffness_dd_mid_A_star_diff = ttk.Labelframe(frame_equivalentdd_mid, borderwidth=2, relief="solid", text="Difference with quad [A*] (%)")
frame_stiffness_dd_mid_A_star_diff.grid(row=6, column=2, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_mid_A_star_diff = ttk.Label(frame_stiffness_dd_mid_A_star_diff, text="", padding=2)
label_stiffness_dd_mid_A_star_diff.grid()

frame_stiffness_dd_mid_D_star_diff = ttk.Labelframe(frame_equivalentdd_mid, borderwidth=2, relief="solid", text="Difference with quad [D*] (%)")
frame_stiffness_dd_mid_D_star_diff.grid(row=6, column=4, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_mid_D_star_diff = ttk.Label(frame_stiffness_dd_mid_D_star_diff, text="", padding=2)
label_stiffness_dd_mid_D_star_diff.grid()

frame_stiffness_dd_mid_B_star_diff = ttk.Labelframe(frame_equivalentdd_mid, borderwidth=2, relief="solid", text="Difference with quad [B*] (%)")
frame_stiffness_dd_mid_B_star_diff.grid(row=6, column=6, rowspan=4, columnspan=2, padx=0, pady=2, sticky='news')
label_stiffness_dd_mid_B_star_diff = ttk.Label(frame_stiffness_dd_mid_B_star_diff, text="", padding=2)
label_stiffness_dd_mid_B_star_diff.grid()

label_equivalentdd_equalmid = tk.Label(frame_equivalentdd_mid, text="ψ: \nφ:", font=('Aerial', 15, 'bold'), foreground='blue', bg='#bae8b0')
label_equivalentdd_equalmid.grid(row=6, column=8, rowspan=4, columnspan=4, padx=5, pady=2,  sticky='news')

label_copyright = ttk.Label(frame_app, text="©2024 Carleton University and National Research Council Canada                                            Version 1.5 (2024-06)                   Peyman.Shabani@carleton.ca, Lucy.Li@nrc-cnrc.gc.ca, Jeremy.Laliberte@carleton.ca", font=('Aerial', 8, 'italic'))
label_copyright.grid(row=5, column=0, columnspan=2,)

frame_app.pack()

app.mainloop()