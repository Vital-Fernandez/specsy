from astropy.table import Table
from specsy.models.tools_pySB99 import find_ionizing_files, ionizing_table

root_folder = '/home/vital/Astrodata/pystarburst99/compiled_spectra_v2'
file_list, mets, timesteps_file = find_ionizing_files(root_folder)

output_file = '/home/vital/Astrodata/pystarburst99/compiled_spectra_v2/pySB99_ionization_Qtable_v2.txt'
nion_table = ionizing_table(file_list, mets, timesteps_file, output_file)