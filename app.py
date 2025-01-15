import streamlit as st
import math
import pandas as pd
import plotly.express as px
from scipy.optimize import fsolve

# Set the title and subtitle with custom HTML and CSS for styling
st.markdown("""
    <style>
    .title {
        font-family: 'Roboto', sans-serif;
        text-align: center;
        font-size: 36px;
        margin-top: 20px;
    }
    .subtitle {
        font-family: 'Roboto', sans-serif;
        text-align: center;
        font-size: 24px;
        margin-top: 10px;
    }
    </style>
    <h1 class="title">Gas Compressibility Factor Calculator</h1>
    <h2 class="subtitle">Hall-Yarborough (1973)</h2>
    """, unsafe_allow_html=True)

# Streamlit app content
st.write("Calculate Z-factor using Hall-Yarborough method")

# Box for Gas Gravity (air)
SG = st.text_input("Gas Gravity (air):")
st.write("Note: Acceptable Range for Gas Gravity is 0.55 - 1.0")

# Validate Gas Gravity input
if SG:
    try:
        SG = float(SG)
        if not (0.55 <= SG <= 1.0):
            st.error("Gas Gravity must be between 0.55 and 1.0")
            SG = None
    except ValueError:
        st.error("Please enter a valid number for Gas Gravity")
        SG = None

# Text input for Pressure (psi)
P = st.text_input("Pressure (psi):")
T = st.text_input("Temperature (deg F):")

# Validate Pressure input
if P:
    try:
        P = float(P)
        if P < 14.7:
            st.error("Pressure must be equal to or above pressure at standard conditions (14.7 psi).")
            P = None
    except ValueError:
        st.error("Please enter a valid number for Pressure.")
        P = None

# Validate Temperature input
if T:
    try:
        T = float(T)
        if T < 60:
            st.error("Temperature must be equal to or above temperature at standard conditions (60 deg F).")
            T = None
    except ValueError:
        st.error("Please enter a valid number for Temperature.")
        T = None

# Text inputs for mole fractions
y_H2S = st.text_input("H2S content (mole frac):")
y_CO2 = st.text_input("CO2 content (mole frac):")
y_N2  = st.text_input("N2 content (mole frac):")

# Initialize variables to store the mole fractions
y_H2S_val, y_CO2_val, y_N2_val = 0, 0, 0

# Validate the inputs
if y_H2S:
    try:
        y_H2S_val = float(y_H2S)
        if not (0 <= y_H2S_val <= 1):
            st.error("H2S content must be between 0 and 1")
            y_H2S_val = 0
    except ValueError:
        st.error("Please enter a valid number for H2S content")
        y_H2S_val = 0

if y_CO2:
    try:
        y_CO2_val = float(y_CO2)
        if not (0 <= y_CO2_val <= 1):
            st.error("CO2 content must be between 0 and 1")
            y_CO2_val = 0
    except ValueError:
        st.error("Please enter a valid number for CO2 content")
        y_CO2_val = 0

if y_N2:
    try:
        y_N2_val = float(y_N2)
        if not (0 <= y_N2_val <= 1):
            st.error("N2 content must be between 0 and 1")
            y_N2_val = 0
    except ValueError:
        st.error("Please enter a valid number for N2 content")
        y_N2_val = 0

# Check if the sum of the mole fractions exceeds 1
total_mole_frac = y_H2S_val + y_CO2_val + y_N2_val
if total_mole_frac > 1:
    st.error("The sum of H2S, CO2, and N2 contents must not exceed 1. Currently, it is {:.2f}".format(total_mole_frac))
else:
    st.write("The sum of H2S, CO2, and N2 contents is {:.2f}".format(total_mole_frac))

# Hall-Yarborough calculations

# Sutton's correlation for pseudo-critical properties
Tpc = 169.2 + 349.5 * SG - 74.0 * SG**2
Ppc = 756.8 - 131.0 * SG - 3.6 * SG**2

# Wichert-Aziz correction for non-hydrocarbon gases
e = 120 * ((y_N2_val + y_CO2_val)**0.9 - (y_N2_val + y_CO2_val)**1.6) + 15 * (y_H2S_val**0.5 - y_H2S_val**4)
Tpc_corr = Tpc - e
Ppc_corr = Ppc * Tpc_corr / (Tpc + y_H2S_val * (1 - y_H2S_val) * (304.2 - Tpc_corr))

# Define the Hall-Yarborough functions
def fy(y, alpha, Ppr, t):
    return (-alpha * Ppr + 
            (y + y ** 2 + y ** 3 - y ** 4) / (1 - y) ** 3 - 
            (14.76 * t - 9.76 * t ** 2 + 4.58 * t ** 3) * y ** 2 + 
            (90.7 * t - 242.2 * t ** 2 + 42.4 * t ** 3) * y ** (2.18 + 2.82 * t))

def Zfac(Tpr, Ppr):
    t = 1 / Tpr
    alpha = 0.06125 * t * math.exp(-1.2 * (1 - t) ** 2)
    y_initial = 0.001  # Initial guess for y
    y = fsolve(fy, y_initial, args=(alpha, Ppr, t))[0]  # Solve for y using fsolve
    return alpha * Ppr / y

# Calculate Z-factor when the button is clicked and inputs are valid
if st.button("Calculate Z-factor"):
    # Check for missing inputs
    missing_inputs = []
    if SG is None:
        missing_inputs.append("Gas Gravity")
    if P is None:
        missing_inputs.append("Pressure")
    if T is None:
        missing_inputs.append("Temperature")
    if not y_H2S:
        missing_inputs.append("H2S content")
    if not y_CO2:
        missing_inputs.append("CO2 content")
    if not y_N2:
        missing_inputs.append("N2 content")

    if missing_inputs:
        st.warning(f"Please provide values for: {', '.join(missing_inputs)}")
    else:
        # Convert temperature from Fahrenheit to Rankine
        T_rankine = T + 459.67

        # Reduced properties
        Pr = P / Ppc_corr
        Tr = T_rankine / Tpc_corr

        z_factor = Zfac(Tr, Pr)
        st.write(f"The gas compressibility factor (Z-factor) is: {z_factor:.3f}")

        # Generate pressure values from 14.7 to (input pressure + 1000) with increments of 200
        pressures = [14.7] + [i for i in range(200, int(P + 1000) + 1, 200)]

        # Calculate Z-factors for each pressure value
        z_factors = [Zfac(T_rankine / Tpc_corr, p / Ppc_corr) for p in pressures]

        # Create a DataFrame for the table
        df = pd.DataFrame({'Pressure (psi)': pressures, 'Z-factor': z_factors})

        # Display the DataFrame as a table
        st.write(df)

        # Plot Z-factor vs. Pressure using Plotly
        fig = px.line(df, x='Pressure (psi)', y='Z-factor', title='Z-factor vs. Pressure')

        # Highlight the input pressure point
        highlight_point = df[df['Pressure (psi)'] == P]

        fig.add_scatter(
            x=highlight_point['Pressure (psi)'],
            y=highlight_point['Z-factor'],
            mode='markers+text',
            marker=dict(size=12, color='red'),
            text=[f"<b>{highlight_point['Z-factor'].values[0]:.4f}</b>"],
            textposition='top center',
            name='Input Pressure'
        )

        st.plotly_chart(fig)