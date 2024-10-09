import openmc
import matplotlib.pyplot as plt
import math

# cross_sections_path = r'/home/davide/openmc_models/CROSS_SECTIONS/endfb-viii.0-hdf5/cross_sections.xml'
# openmc.config['cross_sections'] = cross_sections_path

fw = openmc.Material(name='first_wall')  #first wall, TUNGSTEN
fw.add_nuclide('W182', 26.5e-2)
fw.add_nuclide('W180', 0.12e-2)
fw.add_nuclide('W183', 14.31e-2)
fw.add_nuclide('W184', 30.64e-2)
fw.add_nuclide('W186', 28.43e-2)
fw.set_density('g/cm3', 19.250)
fw.temperature = 900.0  # temperature in Kelvin

inconel = openmc.Material(name='structural_material')  #structural layers, INCONEL 18
inconel.add_element('Al', 0.52e-2, 'wo')
inconel.add_element('C', 0.021e-2, 'wo')
inconel.add_element('Co', 0.11e-2, 'wo')
inconel.add_element('Cr', 19.06e-2, 'wo')
inconel.add_element('Cu', 0.02e-2, 'wo')
inconel.add_element('Fe', 18.15e-2, 'wo')
inconel.add_element('Mo', 3.04e-2, 'wo')
inconel.add_element('Ti', 0.93e-2, 'wo')
inconel.add_element('Nb', 5.08e-2, 'wo')
inconel.add_element('Ni', 53.0e-2, 'wo')
inconel.set_density('g/cm3', 8.19)
inconel.temperature = 900.0

flibe = openmc.Material(name='molten_salt')  #breeding blanket (Li-6 enrichment has been kept constant to 90 %), FLiBe
flibe.add_element('F', 4.)
flibe.add_element('Be', 1.)
flibe.add_nuclide('Li6', 1.8)
flibe.add_nuclide('Li7', 0.2)
flibe.set_density('g/cm3', 1.94)
flibe.temperature = 900.0

nm = openmc.Material(name='neutron_multiplier')  #neutron multiplier, BERILLIUM
nm.add_element('Be', 1.)
nm.set_density('g/cm3', 1.848)
nm.temperature = 900.0

void = openmc.Material(name='chamber')  #hydrogen to simulate vacuum
void.add_element('H', 1.0)
void.set_density('g/cm3', 0.00000001)
void.temperature = 900.0

# Collect the materials together and export to XML
vv_material = inconel
materials = openmc.Materials([fw, vv_material, flibe, nm, void])
materials.export_to_xml()

FW_surface = 2458369.11 # cm2
R_major = 330
r_minor= FW_surface / (4 *  math.pi ** 2 * R_major)

r_str1 = r_minor+0.1
r_flibe1 = r_str1+1
r_be = r_flibe1+2
r_str2 = r_be+1
r_flibe2 = r_str2+3
r_str3 = r_flibe2+100
r_outer = r_str3+3

fw_inner = openmc.ZTorus(x0=0.0, y0=0.0, z0=0.0, a=R_major, b=r_minor, c=r_minor)
str1_inner = openmc.ZTorus(x0=0.0, y0=0.0, z0=0.0, a=R_major, b=r_str1, c=r_str1)
flibe1_inner = openmc.ZTorus(x0=0.0, y0=0.0, z0=0.0, a=R_major, b=r_flibe1, c=r_flibe1)
be_inner = openmc.ZTorus(x0=0.0, y0=0.0, z0=0.0, a=R_major, b=r_be, c=r_be)
str2_inner = openmc.ZTorus(x0=0.0, y0=0.0, z0=0.0, a=R_major, b=r_str2, c=r_str2)
flibe2_iner = openmc.ZTorus(x0=0.0, y0=0.0, z0=0.0, a=R_major, b=r_flibe2, c=r_flibe2)
str3_inner = openmc.ZTorus(x0=0.0, y0=0.0, z0=0.0, a=R_major, b=r_str3, c=r_str3)
str3_outer = openmc.ZTorus(x0=0.0, y0=0.0, z0=0.0, a=R_major, b=r_outer, c=r_outer, boundary_type="vacuum")

inner_region = -fw_inner
fw_region = +fw_inner & -str1_inner
str1_region = +str1_inner & -flibe1_inner
flibe1_region = +flibe1_inner & -be_inner
be_region = +be_inner & -str2_inner
str2_region = +str2_inner & -flibe2_iner
flibe2_region = +flibe2_iner & -str3_inner
str3_region = +str3_inner & -str3_outer

inner_cell = openmc.Cell(1, name='inner', region=inner_region, fill=void)
firstwall_cell = openmc.Cell(2, region=fw_region, fill=fw)
str1_cell = openmc.Cell(3, region=str1_region, fill=inconel)
flibe1_cell = openmc.Cell(4, region=flibe1_region, fill=flibe)
be_cell = openmc.Cell(5, region=be_region, fill=nm)
str2_cell = openmc.Cell(6, region=str2_region, fill=inconel)
flibe2_cell = openmc.Cell(7, region=flibe2_region, fill=flibe)
str3_cell= openmc.Cell(8, region=str3_region, fill=inconel)

universe = openmc.Universe(cells=[inner_cell, firstwall_cell, str1_cell, flibe1_cell, be_cell, str2_cell, flibe2_cell, str3_cell])

geometry = openmc.Geometry()
geometry.root_universe = universe
geometry.export_to_xml()

my_source = openmc.Source()
radius = openmc.stats.Discrete([330], [1])
z_values = openmc.stats.Discrete([0], [1])
angle = openmc.stats.Uniform(a=0, b=math.radians(360))
my_source.space = openmc.stats.CylindricalIndependent(r=radius, phi=angle, z=z_values, origin=(0.0, 0.0, 0.0))
my_source.angle = openmc.stats.Isotropic()
my_source.energy = openmc.stats.muir(e0=14080000.0, m_rat=5.0, kt=20000.0)

settings = openmc.Settings()
settings.batches = 240
settings.particles = 4200000
settings.inactive = 0
settings.run_mode = 'fixed source'
settings.source = my_source
settings.photon_transport = True  # This line is required to switch on photons tracking
settings.export_to_xml()

# adds a tally to record neutron spectra on STR1
energy_bins = openmc.mgxs.GROUP_STRUCTURES['VITAMIN-J-175']
energy_filter = openmc.EnergyFilter(energy_bins)

fw_filter = openmc.CellFilter(firstwall_cell)
str1_filter = openmc.CellFilter(str1_cell)
flibe1_filter = openmc.CellFilter(flibe1_cell)
nm_filter = openmc.CellFilter(be_cell)
str2_filter = openmc.CellFilter(str2_cell)
flibe2_filter = openmc.CellFilter(flibe2_cell)
str3_filter = openmc.CellFilter(str3_cell)

neutron_particle_filter = openmc.ParticleFilter(['neutron'])
photon_particle_filter = openmc.ParticleFilter(['photon'])
electron_particle_filter = openmc.ParticleFilter(['electron'])
positron_particle_filter = openmc.ParticleFilter(['positron'])

# makes a tally to distinguish tritium production from Li6 and Li7
tbr_tally_channel = openmc.Tally(name='TBR channel')
tbr_tally_channel.filters = [flibe1_filter]
tbr_tally_channel.scores = ['(n,Xt)']
tbr_tally_channel.nuclides = [openmc.Nuclide('Li6'), openmc.Nuclide('Li7')]

tbr_tally_tank = openmc.Tally(name='TBR tank')
tbr_tally_tank.filters = [flibe2_filter]
tbr_tally_tank.scores = ['(n,Xt)']
tbr_tally_tank.nuclides = [openmc.Nuclide('Li6'), openmc.Nuclide('Li7')]

# makes a tally for neutron flux on FW
cell_spectra_tally1 = openmc.Tally(name='cell_flux_tally_fw')
cell_spectra_tally1.scores = ['flux']
cell_spectra_tally1.filters = [fw_filter, neutron_particle_filter, electron_particle_filter, positron_particle_filter, photon_particle_filter]

# makes a tally for neutron flux on STR1 (the unit of measurement will then be changed)
cell_spectra_tally2 = openmc.Tally(name='cell_flux_tally_str1')
cell_spectra_tally2.scores = ['flux']
cell_spectra_tally2.filters = [str1_filter, neutron_particle_filter]

# Heating FW
heating_n_fw_tally = openmc.Tally(name="heating fw n")
heating_n_fw_tally.scores = ["heating"]
heating_n_fw_tally.filters = [fw_filter, neutron_particle_filter]

heating_n_fw_tally1 = openmc.Tally(name="heating fw electron")
heating_n_fw_tally1.scores = ["heating"]
heating_n_fw_tally1.filters = [fw_filter, electron_particle_filter]

heating_n_fw_tally2 = openmc.Tally(name="heating fw photon")
heating_n_fw_tally2.scores = ["heating"]
heating_n_fw_tally2.filters = [fw_filter, photon_particle_filter]

heating_n_fw_tally3 = openmc.Tally(name="heating fw positron")
heating_n_fw_tally3.scores = ["heating"]
heating_n_fw_tally3.filters = [fw_filter, positron_particle_filter]

heating_p_fw_tally = openmc.Tally(name="heating fw tot n")
heating_p_fw_tally.scores = ["heating-local"]
heating_p_fw_tally.filters = [fw_filter, neutron_particle_filter]

heating_p_fw_tally1 = openmc.Tally(name="heating fw tot electron")
heating_p_fw_tally1.scores = ["heating-local"]
heating_p_fw_tally1.filters = [fw_filter, electron_particle_filter]

heating_p_fw_tally2 = openmc.Tally(name="heating fw tot photon")
heating_p_fw_tally2.scores = ["heating-local"]
heating_p_fw_tally2.filters = [fw_filter, photon_particle_filter]

heating_p_fw_tally3 = openmc.Tally(name="heating fw tot positron")
heating_p_fw_tally3.scores = ["heating-local"]
heating_p_fw_tally3.filters = [fw_filter, positron_particle_filter]

# Heating STR1
heating_n_str1_tally = openmc.Tally(name="heating str1 n")
heating_n_str1_tally.scores = ["heating"]
heating_n_str1_tally.filters = [str1_filter, neutron_particle_filter]

heating_n_str1_tally1 = openmc.Tally(name="heating str1 electron")
heating_n_str1_tally1.scores = ["heating"]
heating_n_str1_tally1.filters = [str1_filter, electron_particle_filter]

heating_n_str1_tally2 = openmc.Tally(name="heating str1 photon")
heating_n_str1_tally2.scores = ["heating"]
heating_n_str1_tally2.filters = [str1_filter, photon_particle_filter]

heating_n_str1_tally3 = openmc.Tally(name="heating str1 positron")
heating_n_str1_tally3.scores = ["heating"]
heating_n_str1_tally3.filters = [str1_filter, positron_particle_filter]

heating_p_str1_tally = openmc.Tally(name="heating str1 tot n")
heating_p_str1_tally.scores = ["heating-local"]
heating_p_str1_tally.filters = [str1_filter, neutron_particle_filter]

heating_p_str1_tally1 = openmc.Tally(name="heating str1 tot electron")
heating_p_str1_tally1.scores = ["heating-local"]
heating_p_str1_tally1.filters = [str1_filter, electron_particle_filter]

heating_p_str1_tally2 = openmc.Tally(name="heating str1 tot photon")
heating_p_str1_tally2.scores = ["heating-local"]
heating_p_str1_tally2.filters = [str1_filter, photon_particle_filter]

heating_p_str1_tally3 = openmc.Tally(name="heating str1 tot positron")
heating_p_str1_tally3.scores = ["heating-local"]
heating_p_str1_tally3.filters = [str1_filter, positron_particle_filter]

# Heating FLiBe1
heating_n_flibe1_tally = openmc.Tally(name="heating flibe1 n")
heating_n_flibe1_tally.scores = ["heating"]
heating_n_flibe1_tally.filters = [flibe1_filter, neutron_particle_filter]

heating_n_flibe1_tally1 = openmc.Tally(name="heating flibe1 electron")
heating_n_flibe1_tally1.scores = ["heating"]
heating_n_flibe1_tally1.filters = [flibe1_filter, electron_particle_filter]

heating_n_flibe1_tally2 = openmc.Tally(name="heating flibe1 photon")
heating_n_flibe1_tally2.scores = ["heating"]
heating_n_flibe1_tally2.filters = [flibe1_filter, photon_particle_filter]

heating_n_flibe1_tally3 = openmc.Tally(name="heating flibe1 positron")
heating_n_flibe1_tally3.scores = ["heating"]
heating_n_flibe1_tally3.filters = [flibe1_filter, positron_particle_filter]

heating_p_flibe1_tally = openmc.Tally(name="heating flibe1 tot n")
heating_p_flibe1_tally.scores = ["heating-local"]
heating_p_flibe1_tally.filters = [flibe1_filter, neutron_particle_filter]

heating_p_flibe1_tally1 = openmc.Tally(name="heating flibe1 tot electron")
heating_p_flibe1_tally1.scores = ["heating-local"]
heating_p_flibe1_tally1.filters = [flibe1_filter, electron_particle_filter]

heating_p_flibe1_tally2 = openmc.Tally(name="heating flibe1 tot photon")
heating_p_flibe1_tally2.scores = ["heating-local"]
heating_p_flibe1_tally2.filters = [flibe1_filter, photon_particle_filter]

heating_p_flibe1_tally3 = openmc.Tally(name="heating flibe1 tot positron")
heating_p_flibe1_tally3.scores = ["heating-local"]
heating_p_flibe1_tally3.filters = [flibe1_filter, positron_particle_filter]

# Heating Neutron Multiplier
heating_n_nm_tally = openmc.Tally(name="heating nm n")
heating_n_nm_tally.scores = ["heating"]
heating_n_nm_tally.filters = [nm_filter, neutron_particle_filter]

heating_n_nm_tally1 = openmc.Tally(name="heating nm electron")
heating_n_nm_tally1.scores = ["heating"]
heating_n_nm_tally1.filters = [nm_filter, electron_particle_filter]

heating_n_nm_tally2 = openmc.Tally(name="heating nm photon")
heating_n_nm_tally2.scores = ["heating"]
heating_n_nm_tally2.filters = [nm_filter, photon_particle_filter]

heating_n_nm_tally3 = openmc.Tally(name="heating nm positron")
heating_n_nm_tally3.scores = ["heating"]
heating_n_nm_tally3.filters = [nm_filter, positron_particle_filter]

heating_p_nm_tally = openmc.Tally(name="heating nm tot n")
heating_p_nm_tally.scores = ["heating-local"]
heating_p_nm_tally.filters = [nm_filter, neutron_particle_filter]

heating_p_nm_tally1 = openmc.Tally(name="heating nm tot electron")
heating_p_nm_tally1.scores = ["heating-local"]
heating_p_nm_tally1.filters = [nm_filter, electron_particle_filter]

heating_p_nm_tally2 = openmc.Tally(name="heating nm tot photon")
heating_p_nm_tally2.scores = ["heating-local"]
heating_p_nm_tally2.filters = [nm_filter, photon_particle_filter]

heating_p_nm_tally3 = openmc.Tally(name="heating nm tot positron")
heating_p_nm_tally3.scores = ["heating-local"]
heating_p_nm_tally3.filters = [nm_filter, positron_particle_filter]

# Heating STR2
heating_n_str2_tally = openmc.Tally(name="heating str2 n")
heating_n_str2_tally.scores = ["heating"]
heating_n_str2_tally.filters = [str2_filter, neutron_particle_filter]

heating_n_str2_tally1 = openmc.Tally(name="heating str2 electron")
heating_n_str2_tally1.scores = ["heating"]
heating_n_str2_tally1.filters = [str2_filter, electron_particle_filter]

heating_n_str2_tally2 = openmc.Tally(name="heating str2 photon")
heating_n_str2_tally2.scores = ["heating"]
heating_n_str2_tally2.filters = [str2_filter, photon_particle_filter]

heating_n_str2_tally3 = openmc.Tally(name="heating str2 positron")
heating_n_str2_tally3.scores = ["heating"]
heating_n_str2_tally3.filters = [str2_filter, positron_particle_filter]

heating_p_str2_tally = openmc.Tally(name="heating str2 tot n")
heating_p_str2_tally.scores = ["heating-local"]
heating_p_str2_tally.filters = [str2_filter, neutron_particle_filter]

heating_p_str2_tally1 = openmc.Tally(name="heating str2 tot electron")
heating_p_str2_tally1.scores = ["heating-local"]
heating_p_str2_tally1.filters = [str2_filter, electron_particle_filter]

heating_p_str2_tally2 = openmc.Tally(name="heating str2 tot photon")
heating_p_str2_tally2.scores = ["heating-local"]
heating_p_str2_tally2.filters = [str2_filter, photon_particle_filter]

heating_p_str2_tally3 = openmc.Tally(name="heating str2 tot positron")
heating_p_str2_tally3.scores = ["heating-local"]
heating_p_str2_tally3.filters = [str2_filter, positron_particle_filter]

# Heating FLiBe2
heating_n_flibe2_tally = openmc.Tally(name="heating flibe2 n")
heating_n_flibe2_tally.scores = ["heating"]
heating_n_flibe2_tally.filters = [flibe2_filter, neutron_particle_filter]

heating_n_flibe2_tally1 = openmc.Tally(name="heating flibe2 electron")
heating_n_flibe2_tally1.scores = ["heating"]
heating_n_flibe2_tally1.filters = [flibe2_filter, electron_particle_filter]

heating_n_flibe2_tally2 = openmc.Tally(name="heating flibe2 photon")
heating_n_flibe2_tally2.scores = ["heating"]
heating_n_flibe2_tally2.filters = [flibe2_filter, photon_particle_filter]

heating_n_flibe2_tally3 = openmc.Tally(name="heating flibe2 positron")
heating_n_flibe2_tally3.scores = ["heating"]
heating_n_flibe2_tally3.filters = [flibe2_filter, positron_particle_filter]

heating_p_flibe2_tally = openmc.Tally(name="heating flibe2 tot n")
heating_p_flibe2_tally.scores = ["heating-local"]
heating_p_flibe2_tally.filters = [flibe2_filter, neutron_particle_filter]

heating_p_flibe2_tally1 = openmc.Tally(name="heating flibe2 tot electron")
heating_p_flibe2_tally1.scores = ["heating-local"]
heating_p_flibe2_tally1.filters = [flibe2_filter, electron_particle_filter]

heating_p_flibe2_tally2 = openmc.Tally(name="heating flibe2 tot photon")
heating_p_flibe2_tally2.scores = ["heating-local"]
heating_p_flibe2_tally2.filters = [flibe2_filter, photon_particle_filter]

heating_p_flibe2_tally3 = openmc.Tally(name="heating flibe2 tot positron")
heating_p_flibe2_tally3.scores = ["heating-local"]
heating_p_flibe2_tally3.filters = [flibe2_filter, positron_particle_filter]

# Heating STR3
heating_n_str3_tally = openmc.Tally(name="heating str3 n")
heating_n_str3_tally.scores = ["heating"]
heating_n_str3_tally.filters = [str3_filter, neutron_particle_filter]

heating_n_str3_tally1 = openmc.Tally(name="heating str3 electron")
heating_n_str3_tally1.scores = ["heating"]
heating_n_str3_tally1.filters = [str3_filter, electron_particle_filter]

heating_n_str3_tally2 = openmc.Tally(name="heating str3 photon")
heating_n_str3_tally2.scores = ["heating"]
heating_n_str3_tally2.filters = [str3_filter, photon_particle_filter]

heating_n_str3_tally3 = openmc.Tally(name="heating str3 positron")
heating_n_str3_tally3.scores = ["heating"]
heating_n_str3_tally3.filters = [str3_filter, positron_particle_filter]

heating_p_str3_tally = openmc.Tally(name="heating str3 tot n")
heating_p_str3_tally.scores = ["heating-local"]
heating_p_str3_tally.filters = [str3_filter, neutron_particle_filter]

heating_p_str3_tally1 = openmc.Tally(name="heating str3 tot electron")
heating_p_str3_tally1.scores = ["heating-local"]
heating_p_str3_tally1.filters = [str3_filter, electron_particle_filter]

heating_p_str3_tally2 = openmc.Tally(name="heating str3 tot photon")
heating_p_str3_tally2.scores = ["heating-local"]
heating_p_str3_tally2.filters = [str3_filter, photon_particle_filter]

heating_p_str3_tally3 = openmc.Tally(name="heating str3 tot positron")
heating_p_str3_tally3.scores = ["heating-local"]
heating_p_str3_tally3.filters = [str3_filter, positron_particle_filter]

# NEUTRON SPECTRA
# FW
neutronspectra_fw_cell_tally = openmc.Tally(name='neutron_spectra_fw')
neutronspectra_fw_cell_tally.filters = [fw_filter, neutron_particle_filter, energy_filter]
neutronspectra_fw_cell_tally.scores = ['flux']

# STR1
neutronspectra_str1_cell_tally = openmc.Tally(name='neutron_spectra_str1')
neutronspectra_str1_cell_tally.filters = [str1_filter, neutron_particle_filter, energy_filter]
neutronspectra_str1_cell_tally.scores = ['flux']

# FLiBe1
neutronspectra_flibe1_cell_tally = openmc.Tally(name='neutron_spectra_flibe1')
neutronspectra_flibe1_cell_tally.filters = [flibe1_filter, neutron_particle_filter, energy_filter]
neutronspectra_flibe1_cell_tally.scores = ['flux']

# NM
neutronspectra_nm_cell_tally = openmc.Tally(name='neutron_spectra_nm')
neutronspectra_nm_cell_tally.filters = [nm_filter, neutron_particle_filter, energy_filter]
neutronspectra_nm_cell_tally.scores = ['flux']

# STR2
neutronspectra_str2_cell_tally = openmc.Tally(name='neutron_spectra_str2')
neutronspectra_str2_cell_tally.filters = [str2_filter, neutron_particle_filter, energy_filter]
neutronspectra_str2_cell_tally.scores = ['flux']

# FLiBe2
neutronspectra_flibe2_cell_tally = openmc.Tally(name='neutron_spectra_flibe2')
neutronspectra_flibe2_cell_tally.filters = [flibe2_filter, neutron_particle_filter, energy_filter]
neutronspectra_flibe2_cell_tally.scores = ['flux']

# STR3
neutronspectra_str3_cell_tally = openmc.Tally(name='neutron_spectra_str3')
neutronspectra_str3_cell_tally.filters = [str3_filter, neutron_particle_filter, energy_filter]
neutronspectra_str3_cell_tally.scores = ['flux']

# groups the tallies
tallies = openmc.Tallies([tbr_tally_channel, tbr_tally_tank, cell_spectra_tally1, cell_spectra_tally2,
                          heating_n_fw_tally, heating_n_fw_tally1, heating_n_fw_tally2, heating_n_fw_tally3,
                          heating_p_fw_tally, heating_p_fw_tally1, heating_p_fw_tally2, heating_p_fw_tally3,
                          heating_n_str1_tally, heating_n_str1_tally1, heating_n_str1_tally2, heating_n_str1_tally3,
                          heating_p_str1_tally, heating_p_str1_tally1, heating_p_str1_tally2, heating_p_str1_tally3,
                          heating_n_flibe1_tally, heating_n_flibe1_tally1, heating_n_flibe1_tally2, heating_n_flibe1_tally3,
                          heating_p_flibe1_tally, heating_p_flibe1_tally1, heating_p_flibe1_tally2, heating_p_flibe1_tally3,
                          heating_n_nm_tally, heating_n_nm_tally1, heating_n_nm_tally2, heating_n_nm_tally3,
                          heating_p_nm_tally, heating_p_nm_tally1, heating_p_nm_tally2, heating_p_nm_tally3,
                          heating_n_str2_tally, heating_n_str2_tally1, heating_n_str2_tally2, heating_n_str2_tally3,
                          heating_p_str2_tally, heating_p_str2_tally1, heating_p_str2_tally2, heating_p_str2_tally3,
                          heating_n_flibe2_tally, heating_n_flibe2_tally1, heating_n_flibe2_tally2, heating_n_flibe2_tally3,
                          heating_p_flibe2_tally, heating_p_flibe2_tally1, heating_p_flibe2_tally2, heating_p_flibe2_tally3,
                          heating_n_str3_tally, heating_n_str3_tally1, heating_n_str3_tally2, heating_n_str3_tally3,
                          heating_p_str3_tally, heating_p_str3_tally1, heating_p_str3_tally2, heating_p_str3_tally3,
                          neutronspectra_fw_cell_tally, neutronspectra_str1_cell_tally, neutronspectra_flibe1_cell_tally,
                          neutronspectra_nm_cell_tally, neutronspectra_str2_cell_tally, neutronspectra_flibe2_cell_tally,
                          neutronspectra_str3_cell_tally])
tallies.export_to_xml()

# builds the openmc model
my_model = openmc.Model(
    materials=materials, geometry=geometry, settings=settings, tallies=tallies
)