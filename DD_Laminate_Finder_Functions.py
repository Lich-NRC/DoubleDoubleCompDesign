import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import re
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


def stiffness_matrix(theta_degrees_input, Ex, Ey, Es, vx, thickness):
    theta = np.deg2rad(theta_degrees_input)
    # Number of plies and thickness
    ply_number = len(theta)
    
    # Material paraemters
    Ex = Ex * 10**9
    Ey = Ey * 10**9
    Es = Es * 10**9
    vy = Ey * vx / Ex
    thickness = thickness * 10**-3

    # Calculate distance from neutral axis
    h = 0  # Total thickness
    zmin = np.zeros(ply_number)
    zmax = np.zeros(ply_number)
    for i in range(ply_number):
        zmin[i] = h
        h += thickness
        zmax[i] = h

    # Calculating ABD Matrices
    ply_vf = np.zeros(ply_number)  # Volume fraction of each layer
    va_1 = np.zeros(ply_number)
    va_2 = np.zeros(ply_number)
    va_3 = np.zeros(ply_number)
    va_4 = np.zeros(ply_number)

    q = np.array([[Ex / (1 - vx * vy), vy * Ex / (1 - vx * vy), 0],
                [vx * Ey / (1 - vx * vy), Ey / (1 - vx * vy), 0],
                [0, 0, Es]])

    u1 = (1 / 8) * (3 * q[0, 0] + 3 * q[1, 1] + 2 * q[0, 1] + 4 * q[2, 2])
    u2 = 0.5 * (q[0, 0] - q[1, 1])
    u3 = (1 / 8) * (q[0, 0] + q[1, 1] - 2 * q[0, 1] - 4 * q[2, 2])
    u4 = (1 / 8) * (q[0, 0] + q[1, 1] + 6 * q[0, 1] - 4 * q[2, 2])
    u5 = (1 / 8) * (q[0, 0] + q[1, 1] - 2 * q[0, 1] + 4 * q[2, 2])

    A11 = 0
    A12 = 0
    A16 = 0
    A22 = 0
    A26 = 0
    A66 = 0

    for i in range(ply_number):
        ply_vf[i] = thickness / h
        va_1[i] = ply_vf[i] * np.cos(2 * theta[i])
        va_2[i] = ply_vf[i] * np.cos(4 * theta[i])
        va_3[i] = ply_vf[i] * np.sin(2 * theta[i])
        va_4[i] = ply_vf[i] * np.sin(4 * theta[i])

        A11 = A11 + u1 * ply_vf[i] + u2 * va_1[i] + u3 * va_2[i]
        A22 = A22 + u1 * ply_vf[i] - u2 * va_1[i] + u3 * va_2[i]
        A12 = A12 + u4 * ply_vf[i] - u3 * va_2[i]
        A66 = A66 + u5 * ply_vf[i] - u3 * va_2[i]
        A16 = A16 + 0.5 * u2 * va_3[i] + u3 * va_4[i]
        A26 = A26 + 0.5 * u2 * va_3[i] - u3 * va_4[i]

    A_star = np.array([[A11, A12, A16],
                    [A12, A22, A26],
                    [A16, A26, A66]])
    A = A_star * h

    va_1_star = (A_star[0, 0] - A_star[1, 1]) / (2 * u2)
    va_2_star = (A_star[0, 0] + A_star[1, 1] - 2 * u1) / (2 * u3)
    va_3_star = (A_star[2, 0] + A_star[2, 1]) / u2
    va_4_star = (A_star[2, 0] - A_star[2, 1]) / (2 * u3)

    # Initialize variables for D matrix calculations
    vd_1 = np.zeros(ply_number)
    vd_2 = np.zeros(ply_number)
    vd_3 = np.zeros(ply_number)
    vd_4 = np.zeros(ply_number)
    D11 = 0
    D12 = 0
    D16 = 0
    D22 = 0
    D26 = 0
    D66 = 0

    for i in range(ply_number):
        ply_vf[i] = thickness / h
        Q11 = u1 + u2 * np.cos(2 * theta[i]) + u3 * np.cos(4 * theta[i])
        Q22 = u1 - u2 * np.cos(2 * theta[i]) + u3 * np.cos(4 * theta[i])
        Q12 = u4 - u3 * np.cos(4 * theta[i])
        Q66 = u5 - u3 * np.cos(4 * theta[i])
        Q16 = 0.5 * np.sin(2 * theta[i]) * u2 + np.sin(4 * theta[i]) * u3
        Q26 = 0.5 * np.sin(2 * theta[i]) * u2 - np.sin(4 * theta[i]) * u3

        fd_11 = lambda z: (z ** 2) * Q11
        D11 += quad(fd_11, zmin[i] - h / 2, zmax[i] - h / 2)[0]

        fd_22 = lambda z: (z ** 2) * Q22
        D22 += quad(fd_22, zmin[i] - h / 2, zmax[i] - h / 2)[0]

        fd_12 = lambda z: (z ** 2) * Q12
        D12 += quad(fd_12, zmin[i] - h / 2, zmax[i] - h / 2)[0]

        fd_66 = lambda z: (z ** 2) * Q66
        D66 += quad(fd_66, zmin[i] - h / 2, zmax[i] - h / 2)[0]

        fd_16 = lambda z: (z ** 2) * Q16
        D16 += quad(fd_16, zmin[i] - h / 2, zmax[i] - h / 2)[0]

        fd_26 = lambda z: (z ** 2) * Q26
        D26 += quad(fd_26, zmin[i] - h / 2, zmax[i] - h / 2)[0]

    D = np.array([[D11, D12, D16],
                [D12, D22, D26],
                [D16, D26, D66]])
    D_star = 12 / h ** 3 * D
    vd_1_star = (D_star[0, 0] - D_star[1, 1]) / (2 * u2)
    vd_2_star = (D_star[0, 0] + D_star[1, 1] - 2 * u1) / (2 * u3)
    vd_3_star = (D_star[2, 0] + D_star[2, 1]) / u2
    vd_4_star = (D_star[2, 0] - D_star[2, 1]) / (2 * u3)

    B11 = 0
    B12 = 0
    B16 = 0
    B22 = 0
    B26 = 0
    B66 = 0
    vb_1 = np.zeros(ply_number)
    vb_2 = np.zeros(ply_number)
    vb_3 = np.zeros(ply_number)
    vb_4 = np.zeros(ply_number)
    hi = 0

    for i in range(ply_number):
        hi += thickness
        Q11 = u1 + u2 * np.cos(2 * theta[i]) + u3 * np.cos(4 * theta[i])
        Q22 = u1 - u2 * np.cos(2 * theta[i]) + u3 * np.cos(4 * theta[i])
        Q12 = u4 - u3 * np.cos(4 * theta[i])
        Q66 = u5 - u3 * np.cos(4 * theta[i])
        Q16 = 0.5 * np.sin(2 * theta[i]) * u2 + np.sin(4 * theta[i]) * u3
        Q26 = 0.5 * np.sin(2 * theta[i]) * u2 - np.sin(4 * theta[i]) * u3

        fd_11 = lambda z: z * Q11
        B11 += quad(fd_11, zmin[i] - h / 2, zmax[i] - h / 2)[0]

        fd_22 = lambda z: z * Q22
        B22 += quad(fd_22, zmin[i] - h / 2, zmax[i] - h / 2)[0]

        fd_12 = lambda z: z * Q12
        B12 += quad(fd_12, zmin[i] - h / 2, zmax[i] - h / 2)[0]

        fd_66 = lambda z: z * Q66
        B66 += quad(fd_66, zmin[i] - h / 2, zmax[i] - h / 2)[0]

        fd_16 = lambda z: z * Q16
        B16 += quad(fd_16, zmin[i] - h / 2, zmax[i] - h / 2)[0]

        fd_26 = lambda z: z * Q26
        B26 += quad(fd_26, zmin[i] - h / 2, zmax[i] - h / 2)[0]

    B = np.array([[B11, B12, B16],
                [B12, B22, B26],
                [B16, B26, B66]])
    B_star = 2/h**2 * B
    
    return A, B, D, A_star, B_star, D_star, va_1_star, va_2_star, vd_1_star, vd_2_star    



# Function to calculate the parabola
def calculate_parabola(x_values, y_values):
    coefficients = np.polyfit(x_values, y_values, 2)
    a, b, c = coefficients
    return a, b, c

plt.rcParams['figure.figsize'] = (6, 6)

# Function to plot the parabola
def plot_parabola(v1_quad, v2_quad, v1_dd_a, v2_dd_a, v1_dd_d, v2_dd_d, v1_dd_mid, v2_dd_mid, inplane):
    plt.clf()
    #plt.figure(figsize=(6, 6))
    plt.grid(True)
    x_values = [0, 1, -1]
    y_values = [-1, 1, 1]
    a, b, c = calculate_parabola(x_values, y_values)
    x = np.linspace(-1, 1, 400)
    y = a * x ** 2 + b * x + c
    plt.plot(x, y, color='black')

    x_values = [-1, -0.5, 0]
    y_values = [1, 0, 1]
    a, b, c = calculate_parabola(x_values, y_values)
    x = np.linspace(-1, 0, 400)
    y = a * x ** 2 + b * x + c
    plt.plot(x, y, color='black')

    x_values = [0, 0.5, 1]
    y_values = [1, 0, 1]
    a, b, c = calculate_parabola(x_values, y_values)
    x = np.linspace(0, 1, 400)
    y = a * x ** 2 + b * x + c
    plt.plot(x, y, color='black')

    if inplane==True:
        plt.xlabel("V$^{*}$$_{A1}$", fontsize=14)
        plt.ylabel("V$^{*}$$_{A2}$", fontsize=14)
        plt.title('Double-Double design range (in-plane)')
    else:
        plt.xlabel("V$^{*}$$_{D1}$", fontsize=14)
        plt.ylabel("V$^{*}$$_{D2}$", fontsize=14)
        plt.title('Double-Double design range (flexural)')
    
    plt.subplots_adjust(left=0.15, right=0.9, bottom=0.3, top=0.9)

    plt.scatter([v1_quad], [v2_quad], color='#ff0000', marker='x',  s=50, linewidths=2, label='Given QUAD laminate')
    plt.scatter([v1_dd_a], [v2_dd_a], color='#fcba03', label='Equivalent DD with similar [A*]', marker='o', facecolors='none', s=150, linewidths=2)
    plt.scatter([v1_dd_d], [v2_dd_d], color='#037ffc', label='Equivalent DD with similar [D*]', marker='o', facecolors='none', s=150, linewidths=2)
    plt.scatter([v1_dd_mid], [v2_dd_mid], color='#12b800', label='Equivalent DD with similar [A*] and [D*]', marker='o', facecolors='none', s=150, linewidths=2)
    #plt.legend()
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=1)
    plt.show()

def round_numbers_in_string(input_string):
    # Extract numbers from the string using regex
    numbers = re.findall(r"-?\d+\.\d+", input_string)
    
    # Round numbers and format them
    rounded_numbers = []
    for number in numbers:
        rounded_number = round(float(number), 1)
        if rounded_number.is_integer():
            rounded_numbers.append(str(int(rounded_number)))
        else:
            rounded_numbers.append(str(rounded_number))
    
    # Reconstruct the string
    parts = re.split(r"(-?\d+\.\d+)", input_string)
    for i, part in enumerate(parts):
        if re.match(r"-?\d+\.\d+", part):
            parts[i] = rounded_numbers.pop(0)
    
    return ''.join(parts)
