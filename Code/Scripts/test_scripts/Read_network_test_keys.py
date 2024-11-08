# Test network see duplicate keys
from netCDF4 import Dataset

# Open the NetCDF file
filepath = r"\\dc02\Project\HL-P24050\05_Analysis\02_Model\HHSK_08nov\dflowfm\network.nc"
dataset = Dataset(filepath, 'r')

# Print the file structure
print("Dimensions:")
for dim in dataset.dimensions.values():
    print(dim)

print("\nVariables:")
for var in dataset.variables.values():
    print(var)

print("\nGlobal Attributes:")
for attr in dataset.ncattrs():
    print(f"{attr}: {dataset.getncattr(attr)}")

# Check for the Conventions attribute
conventions = dataset.getncattr('Conventions') if 'Conventions' in dataset.ncattrs() else "No Conventions attribute found"
print(f"\nConventions: {conventions}")

if 'links' in dataset.variables:
    link_id_var = dataset.variables['links']  # Access the link_id variable
    link_ids = link_id_var[:].tolist()  # Convert it to a list for easier handling

    print("\nLink IDs:")
    print(link_ids)
else:
    print("Link ID variable not found in the dataset.")

# Close the dataset
dataset.close()
