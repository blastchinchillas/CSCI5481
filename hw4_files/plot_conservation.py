import numpy as np
import matplotlib.pyplot as plt

# Load conservation rates from the text file generated in step (1)
conservation_rates_file = "solution-problem-1.txt"
with open(conservation_rates_file, "r") as file:
    conservation_rates = [float(line.strip()) for line in file]

# Define the window size for smoothing
window_size = 50

# Perform sliding window average for smoothing
smoothed_conservation_rates = np.convolve(conservation_rates, np.ones(window_size)/window_size, mode='valid')

# Generate x-axis values (positions in the gapped alignment)
positions = np.arange(1, len(smoothed_conservation_rates) + 1)

# Plot the smoothed variability against position
plt.plot(positions, smoothed_conservation_rates, label='Smoothed Variability', color='blue')
plt.xlabel('Position in Gapped Alignment')
plt.ylabel('Smoothed Variability')
plt.title('Variability in Gapped Alignment with Smoothing')
plt.legend()
plt.grid(True)

# Save the figure to a PDF file
output_figure_path = "solution-problem-2.pdf"
plt.savefig(output_figure_path)

plt.show()
