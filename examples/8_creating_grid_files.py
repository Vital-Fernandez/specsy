import numpy as np
import h5netcdf

# Create an example dictionary of numpy arrays
data_dict = {
    'array1': np.random.rand(10, 5),  # Random array of shape (10, 5)
    'array2': np.random.rand(10, 5)  # Another random array of shape (10, 5)
}

# Specify the path for the HDF5 file
filename = 'example.h5'

# Create a new HDF5 file
with h5netcdf.File(filename, 'w') as f:
    # Create dimensions
    f.dimensions['m'] = data_dict['array1'].shape[0]
    f.dimensions['n'] = data_dict['array1'].shape[1]

    # Iterate over the dictionary and create a variable for each array
    for key, values in data_dict.items():
        # Create a variable with the same name as the key
        var = f.create_variable(key, ('m', 'n'), data=values)
        # Set the variable attributes if needed
        var.attrs['description'] = f'Data from {key}'

with h5netcdf.File(filename, 'r') as f:
    # Create a new dictionary to store the loaded data
    loaded_data_dict = {}

    # Iterate over variables in the file
    for var_name in f.variables:
        # Add the data to the dictionary
        loaded_data_dict[var_name] = f.variables[var_name][...]

# Print the loaded data dictionary to check if it matches the original
print(f'array1: {np.all(loaded_data_dict["array1"] == data_dict["array1"])}')
print(f'array2: {np.all(loaded_data_dict["array2"] == data_dict["array2"])}')

# for key, value in loaded_data_dict.items():
#     print(f'{key}: {value}')