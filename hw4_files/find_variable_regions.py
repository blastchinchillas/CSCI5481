import numpy as np

# Load conservation rates from the text file generated in step (1)
conservation_rates_file = "solution-problem-1.txt"
with open(conservation_rates_file, "r") as file:
    conservation_rates = np.array([float(line.strip()) for line in file])

# Set a threshold for identifying variable regions
threshold = 0.5  # Adjust the threshold as needed

# Identify variable regions based on the threshold
variable_regions = np.where(conservation_rates < threshold)[0]

# Initialize variables for tracking variable regions
start_coordinates = [variable_regions[0]]
end_coordinates = []

# Identify start and end coordinates of variable regions
for i in range(len(variable_regions) - 1):
    if variable_regions[i] + 1 != variable_regions[i + 1]:
        end_coordinates.append(variable_regions[i])
        start_coordinates.append(variable_regions[i + 1])    

# Add the last coordinate
end_coordinates.append(variable_regions[-1])

# Write start and end coordinates to a tab-delimited text file
output_file_path = "solution-problem-3.txt"
with open(output_file_path, "w") as file:
    for start, end in zip(start_coordinates, end_coordinates):
        # only write the region longer than 9
        if int(end) - int(start) >= 9 :
            file.write(f"{start + 1}\t{end + 1}\n")  # Adding 1 to convert from 0-based to 1-based coordinates

print(f"Variable regions saved to {output_file_path}")

