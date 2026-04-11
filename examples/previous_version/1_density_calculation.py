import lime
import specsy
import pyneb

pyneb.lab

logs_file = f'my/logs/file/address.txt'
lines_log = lime.load_lines_log(logs_file)

lines = ['S2_6716A', 'S2_6731A']

flux_dict = specsy.flux_distribution(lines_log, flux_type='intg')

nSII_array = specsy.truncated_SII_density_dist(SII_lines=lines, flux_dict=flux_dict)

